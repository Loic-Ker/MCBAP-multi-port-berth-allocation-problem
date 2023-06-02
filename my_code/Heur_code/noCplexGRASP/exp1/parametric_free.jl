using UnPack 
using Statistics
import XLSX
using CSV, Tables
using DataFrames
using Random
using StatsBase
include("../../../MBAP_INST.jl")
include("../get_iterations.jl")
include("../MBAP_SOL.jl")
include("../toolsMatrixTimes.jl")
include("../check_solution.jl")
include("../utilInit.jl")
#include("../localCPLEX.jl")
include("../localSearch.jl")

mutable struct NewVisit
    n::Int64    # ship number
    c::Int64    # visit number
    b::Int64  #postion
    t::Int64    # port
    cost::Float64    
    distance::Float64
    time::Float64   
    constrained::Bool
    store::ToStoreVisit 
end


function findConstrainedPos(inst::Instance, sol::Sol, port::Int, t1::Int, t2::Int, n::Int, c::Int, l::Int)
    # find ship visits during t1-t2
    pos = Tuple{Int, Int}[]
    for n_ in 1:inst.N
        for (c_,p) in enumerate(inst.Pi[n_])
            if p == port
                if sol.visits[n_][c_].planned
                    @unpack b,t = sol.visits[n_][c_]
                    hand = ceil(Int, inst.h[n_][c_][b])
                    if doOverlap(t1, t2, t, t + hand)
                        len = ceil(Int, inst.shipsIn[n_].l/inst.qli)
                        # check 12 possible positions (16-4)
                        # this positions may still be infeasible
                        X = [b - l, b, b + len - l, b + len]
                        # Y = [t - h, t, t + hand - h, t + hand]
                        for x in b-l:b+len
                            if 1 <= x && x + l <= inst.Bp[p]
                                h = ceil(Int, inst.h[n][c][x])
                                for y in [t - h, t + hand]
                                    if t1 <= y && y + h <= t2
                                        
                                        push!(pos, (x, y))
                                    end
                                end
                            end
                        end
                        for x in [b - l, b - l-1, b + len, b+len+1]
                            if 1 <= x && x + l <= inst.Bp[p]
                                h = ceil(Int, inst.h[n][c][x])
                                for y in t-h:t+hand
                                    if t1 <= y && y + h <= t2
                                        push!(pos, (x, y))
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    # also consider external ships and quay limits
    for extS in inst.shipsOut[port]
        @unpack time, berth, hand, length = extS
        if doOverlap(t1, t2, time, time + hand)
            # check 12 possible positions (16-4)
            # this positions may still be infeasible
            X = [berth - l, berth, berth + length - l, berth + length]
            # Y = [t - h, t, t + hand - h, t + hand]
            for (ix, x) in enumerate(X)
                if 1 <= x && x + l <= inst.Bp[port]
                    h = ceil(Int, inst.h[n][c][x])
                    if ix in [2,3]
                        if t1 <= time - h && time <= t2
                            push!(pos, (x, time - h))
                        end
                        if t1 <= time + hand && time + hand + h <= t2
                            push!(pos, (x, time + hand))
                        end
                    end
                    if ix in [1,4]
                        if t1 <= time && time + h <= t2
                            push!(pos, (x, time))
                        end
                        if t1 <= time + hand - h && time + hand <= t2
                            push!(pos, (x, time + hand - h))
                        end
                    end
                end
            end
        end
    end
    #for t in t1:t2
    #    h0 = ceil(Int, inst.h[n][c][1])
    #    hB = ceil(Int, inst.h[n][c][inst.Bp[port]-l-1])
    #    if t + h0 <= t2
    #        push!(pos, (1, t))
    #        push!(pos, (inst.Bp[port]-l-1, t))
    #    end
    #end
    return pos
end

# compute feasible positions (b,t) in different ways :  
## Order by best costs
## Get random ones
## Get the ones with less feasible solutions
## Get the ones closer to other ports visit
## Ordered by port visit number et the end
function SelectNewVisitPerShip(inst::Instance, sol::Sol, paramchosen::ChosenParameters, paramfixed::FixedParameters, n, c, p, printdetails) 
    @unpack N, h, qli, shipsIn, Pi, Bp, delta, dist, T, maxC, maxT, Pc = inst
    @unpack M = sol
    @unpack OneBoatCost, OneBoatDistance, OneBoatTime, AllBoatsCost, AllBoatsCount, AllBoatsTime, LocalSearchRandom, LocalSearchBoat = paramfixed
    l = ceil(Int, shipsIn[n].l/qli)
    new_visits = Vector{NewVisit}()
    new_visits_constrained = Vector{NewVisit}()
    max_dist=0
    max_cost=0
    max_time=0
    min_dist=1000000000
    min_cost=1000000000
    min_time=1000000000
    max_dist_const=0
    max_cost_const=0
    max_time_const=0
    min_dist_const=1000000000
    min_cost_const=1000000000
    min_time_const=1000000000
    #find_feasible = false
    #sch = sol.visits[n]

    t1=sol.visits[n][c].minT
    t2=sol.visits[n][c].maxT
    
    bestPos=shipsIn[n].Bi[c]/qli
    ## Check all the possible solutions
    for b in 1:Bp[p]-l
        hand = ceil(Int, h[n][c][b])
        for t in t1:t2
            if t + hand <= t2
                if M[n][c][b,t+1]
                    this_cost, delay_cost, waiting_cost, penalty, handling_cost, fuel_cost, feas = computeCostPosSol(inst, n, c, b, t,sol)
                    time=penalty
                    if feas
                        distance = abs(bestPos-b)/l
                        if distance>=max_dist
                            max_dist=distance
                        end
                        if this_cost>=max_cost
                            max_cost=this_cost
                        end
                        if time>=max_time
                            max_time=time
                        end
                        if distance<=min_dist
                            min_dist=distance
                        end
                        if this_cost<=min_cost
                            min_cost=this_cost
                        end
                        if time<=min_time
                            min_time=time
                        end
                        push!(new_visits, NewVisit(n,c,b,t,this_cost,distance,time, false, ToStoreVisit(SplitCosts(ceil(Int, this_cost), ceil(Int,delay_cost), ceil(Int,waiting_cost), ceil(Int,penalty), ceil(Int,handling_cost), ceil(Int,fuel_cost)),0,"","")))
                    end
                end
            end
        end
    end  
    
    ## Check the constrained ones
    #pos = findConstrainedPos(inst, sol, p, t1, t2, n, c, l)
    #for (b,t) in pos
    #    hand = ceil(Int, h[n][c][b])
    #    if t + hand <= t2
    #        if M[n][c][b,t+1]
    ##            this_cost, delay_cost, waiting_cost, penalty, handling_cost, fuel_cost, feas = computeCostPosSol(inst, n, c, b, t,sol)
    #            time=penalty
    #            if feas
    #                distance = abs(bestPos-b)/l
    #                if distance>=max_dist_const
    #                    max_dist=distance
    #                end
    #                if this_cost>=max_cost_const
    #                    max_cost_const=this_cost
    #                end
    #                if time>=max_time_const
    #                    max_time_const=time
    #                end
    #                if distance<=min_dist_const
    #                    min_dist_const=distance
    #                end
    #                if this_cost<=min_cost_const
    #                    min_cost_const=this_cost
    #                end
    #                if time<=min_time_const
    #                    min_time_const=time
    #                end
    #                push!(new_visits_constrained, NewVisit(n,c,b,t,this_cost,distance,time, true, ToStoreVisit(SplitCosts(ceil(Int, this_cost), ceil(Int,delay_cost), ceil(Int,waiting_cost), ceil(Int,penalty), ceil(Int,handling_cost), ceil(Int,fuel_cost)),0,"","")))
    #            end
    #        end
    #    end
    #end
    pos_chosen=NewVisit(0,0,0,0,0.0,0.0,0.0, false, ToStoreVisit(SplitCosts(0,0,0,0,0,0),0,"",""))
    pos_chosen_constrained=NewVisit(0,0,0,0,0.0,0.0,0.0, false, ToStoreVisit(SplitCosts(0,0,0,0,0,0),0,"",""))
    tactic=paramchosen.TacticOneBoat
    if length(new_visits)>0
        if tactic=="cost"
            alpha_value= OneBoatCost
            list_available_shippos = Vector{NewVisit}()
            for visit in new_visits
                if visit.constrained==false
                    if visit.cost<=min_cost+alpha_value*(max_cost-min_cost)
                        visit.store.tacticBoat="cost"
                        push!(list_available_shippos, visit)
                    end
                end
            end
            pos_chosen = list_available_shippos[rand(1:length(list_available_shippos))]
        end
        if tactic=="dist"
            alpha_value= OneBoatDistance
            list_available_shippos = Vector{NewVisit}()
            for visit in new_visits
                if visit.constrained==false
                    if visit.distance<=min_dist+alpha_value*(max_dist-min_dist)
                        visit.store.tacticBoat="dist"
                        push!(list_available_shippos, visit)
                    end                   
                end
            end
            pos_chosen = list_available_shippos[rand(1:length(list_available_shippos))]
        end
        if tactic=="time"
            alpha_value= OneBoatTime
            list_available_shippos = Vector{NewVisit}()
            for visit in new_visits
                if visit.constrained==false
                    if visit.time<=min_time+alpha_value*(max_time-min_time)
                        visit.store.tacticBoat="time"
                        push!(list_available_shippos, visit)
                    end                   
                end
            end
            pos_chosen = list_available_shippos[rand(1:length(list_available_shippos))]
        end
    end
    
    #if length(new_visits_constrained)>0
    #    if tactic=="cost"
    #        alpha_value= OneBoatCost
    #        list_available_shippos = Vector{NewVisit}()
    #        for visit in new_visits_constrained
    #            if visit.constrained==true
    #                if visit.cost<=min_cost_const+alpha_value*(max_cost_const-min_cost_const)
    #                    visit.store.tacticBoat="cost"
    #                    push!(list_available_shippos, visit)
    #                end
    #            end
    #        end
    #        pos_chosen_constrained = list_available_shippos[rand(1:length(list_available_shippos))]
    #    end
    #    if tactic=="dist"
    #        alpha_value= OneBoatDistance
    #        list_available_shippos = Vector{NewVisit}()
    #        for visit in new_visits_constrained
    #            if visit.constrained==true
    #                if visit.distance<=min_dist_const+alpha_value*(max_dist_const-min_dist_const)
    #                    visit.store.tacticBoat="dist"
    #                    push!(list_available_shippos, visit)
    #                end                   
    #            end
    #        end
    #        pos_chosen_constrained = list_available_shippos[rand(1:length(list_available_shippos))]
    #    end
    #    if tactic=="time"
    #        alpha_value= OneBoatTime
    #        list_available_shippos = Vector{NewVisit}()
    #        for visit in new_visits_constrained
    #            if visit.constrained==true
    #                if visit.time<=min_time_const+alpha_value*(max_time_const-min_time_const)
    #                    visit.store.tacticBoat="time"
    #                    push!(list_available_shippos, visit)
    #                end                   
    #            end
    #        end
    #        pos_chosen_constrained = list_available_shippos[rand(1:length(list_available_shippos))]
    #    end
    #end
    
    if length(new_visits)+length(new_visits_constrained)>0
        return pos_chosen, pos_chosen_constrained, true, length(new_visits)+length(new_visits_constrained)
    else
        return pos_chosen, pos_chosen_constrained, false, 0
    end
end


# compute feasible positions (b,t) in different ways :  
## Order by best costs
## Get random ones
## Get the ones with less feasible solutions
## Get the ones closer to other ports visit
## Ordered by port visit number et the end
function SelectNewVisitPerShipReconstruct(inst::Instance, sol::Sol, paramchosen::ChosenParameters, paramfixed::FixedParameters, n, c, p, printdetails) 
    @unpack N, h, qli, shipsIn, Pi, Bp, delta, dist, T, maxC, maxT, Pc = inst
    @unpack M = sol
    @unpack OneBoatCost, OneBoatDistance, OneBoatTime, AllBoatsCost, AllBoatsCount, AllBoatsTime, LocalSearchRandom, LocalSearchBoat = paramfixed
    l = ceil(Int, shipsIn[n].l/qli)
    new_visits = Vector{NewVisit}()
    new_visits_constrained = Vector{NewVisit}()
    max_dist=0
    max_cost=0
    max_time=0
    min_dist=1000000000
    min_cost=1000000000
    min_time=1000000000
    #max_dist_const=0
    #max_cost_const=0
    #min_dist_const=1000000000
    #min_cost_const=1000000000
    #find_feasible = false
    #sch = sol.visits[n]

    t1=sol.visits[n][c].minT
    t2=sol.visits[n][c].maxT
    
    bestPos=shipsIn[n].Bi[c]/qli
    ## Check all the possible solutions

    #print('\n')
    #print("The visit $n $c :")
    minT=20*maxT
    for b in 1:Bp[p]-l
        t=t1
        find_feasible=true
        while find_feasible
            if M[n][c][b,t+1]
                if t<minT
                    minT=t
                end
                #print('\n')
                #print("The maximum given time :")
                #print('\n')
                #print(t2)
                #print("The time found :")
                #print('\n')
                #print(t)
                this_cost, delay_cost, waiting_cost, penalty, handling_cost, fuel_cost, feas = computeCostPosSol(inst, n, c, b, t,sol)
                time=penalty
                find_feasible = false
                if feas
                    distance = abs(bestPos-b)/l
                    if distance>=max_dist
                        max_dist=distance
                    end
                    if this_cost>=max_cost
                        max_cost=this_cost
                    end
                    if time>=max_time
                        max_time=time
                    end
                    if distance<=min_dist
                        min_dist=distance
                    end
                    if this_cost<=min_cost
                        min_cost=this_cost
                    end
                    if time<=min_time
                        min_time=time
                    end
                    push!(new_visits, NewVisit(n,c,b,t,this_cost,distance,time, false, ToStoreVisit(SplitCosts(ceil(Int, this_cost), ceil(Int,delay_cost), ceil(Int,waiting_cost), ceil(Int,penalty), ceil(Int,handling_cost), ceil(Int,fuel_cost)),0,"","")))
                end
            end
            t=t+1
            if t>minT
                find_feasible=false
            end
        end
    end  
    pos_chosen=NewVisit(0,0,0,0,0.0,0.0,0.0, false, ToStoreVisit(SplitCosts(0,0,0,0,0,0),0,"",""))
    pos_chosen_constrained=NewVisit(0,0,0,0,0.0,0.0,0.0, false, ToStoreVisit(SplitCosts(0,0,0,0,0,0),0,"",""))
    tactic=paramchosen.TacticOneBoat
    if length(new_visits)>0
        if tactic=="cost"
            alpha_value= OneBoatCost
            list_available_shippos = Vector{NewVisit}()
            for visit in new_visits
                if visit.constrained==false
                    if visit.cost<=min_cost+alpha_value*(max_cost-min_cost)
                        visit.store.tacticBoat="cost"
                        push!(list_available_shippos, visit)
                    end
                end
            end
            pos_chosen = list_available_shippos[rand(1:length(list_available_shippos))]
        end
        if tactic=="dist"
            alpha_value= OneBoatDistance
            list_available_shippos = Vector{NewVisit}()
            for visit in new_visits
                if visit.constrained==false
                    if visit.distance<=min_dist+alpha_value*(max_dist-min_dist)
                        visit.store.tacticBoat="dist"
                        push!(list_available_shippos, visit)
                    end                   
                end
            end
            pos_chosen = list_available_shippos[rand(1:length(list_available_shippos))]
        end
        if tactic=="time"
            alpha_value= OneBoatTime
            list_available_shippos = Vector{NewVisit}()
            for visit in new_visits
                if visit.constrained==false
                    if visit.time<=min_time+alpha_value*(max_time-min_time)
                        visit.store.tacticBoat="time"
                        push!(list_available_shippos, visit)
                    end                   
                end
            end
            pos_chosen = list_available_shippos[rand(1:length(list_available_shippos))]
        end
    end
    if length(new_visits)+length(new_visits_constrained)>0
        return pos_chosen, pos_chosen_constrained, true, length(new_visits)+length(new_visits_constrained)
    else
        return pos_chosen, pos_chosen_constrained, false, 0
    end
end

# compute feasible positions (b,t) in different ways :  
## Order by best costs
## Get random ones
## Get the ones with less feasible solutions
## Get the ones closer to other ports visit
## Ordered by port visit number et the end
function SelectNewVisitAllShips(inst::Instance, sol::Sol, paramchosen::ChosenParameters,paramfixed::FixedParameters) #, placed::Vector{Tuple{Int64,Int64}})
    #TODO: t1 cannot be greater than t2!
    # adapt to given planned previous visits for each ship
    @unpack N, h, qli, shipsIn, Pi, Bp, delta, dist, T, maxC, maxT = inst
    @unpack M = sol
    @unpack OneBoatCost, OneBoatDistance, OneBoatTime, AllBoatsCost, AllBoatsCount, AllBoatsTime, LocalSearchRandom, LocalSearchBoat = paramfixed
    # C = maximum(length.(Pi))
    new_visits = Vector{NewVisit}()
    new_visits_constrained = Vector{NewVisit}()
    max_count=0
    max_cost=0
    max_time = 0
    min_count=1000000000
    min_cost=1000000000
    min_time=1000000000
    max_count_const=0
    max_cost_const=0
    max_time_const = 0
    min_count_const=1000000000
    min_cost_const=1000000000
    min_time_const=1000000000
    count_available=Dict()
    for n in 1:N
        count_available[n]=Dict()
        for (c,p) in enumerate(Pi[n])
            sch = sol.visits[n]
            if sch[c].planned == false
                visit, visit_constrained, feasible, count_visit=SelectNewVisitPerShip(inst, sol, paramchosen, paramfixed, n, c, p, false)
                if feasible
                    count_available[n][c]=count_visit
                    if count_visit>max_count
                        max_count=count_visit
                    end
                    if count_visit<min_count
                        min_count=count_visit
                    end
                    if visit.n!=0
                        this_cost = visit.cost
                        this_time = visit.time
                        if this_cost>max_cost
                            max_cost=this_cost
                        end
                        if this_cost<min_cost
                            min_cost=this_cost
                        end
                        if this_time>max_time
                            max_time=this_time
                        end
                        if this_time<min_time
                            min_time=this_time
                        end
                        push!(new_visits,visit)
                    end
                    #if visit_constrained.n!=0
                    #    this_cost = visit_constrained.cost
                    #    if this_cost>=max_cost_const
                    #        max_cost_const=this_cost
                    #    end
                    #    if this_cost<=min_cost_const
                    #        min_cost_const=this_cost
                    #    end
                    #    if this_time>max_time_const
                    #        max_time_const=this_time
                    #    end
                    #    if this_time<min_time_const
                    #        min_time_const=this_time
                    #    end
                    #    push!(new_visits_constrained,visit_constrained)
                    #end
                end
            end
        end
    end
    pos_chosen=[]
    if length(new_visits)>0
        tactic=paramchosen.TacticAllBoats
        if tactic=="cost"
            alpha_value = AllBoatsCost
            list_available_shippos = Vector{NewVisit}()
            for visit in new_visits
                if visit.cost<=min_cost+alpha_value*(max_cost-min_cost)
                    visit.store.tacticAll="cost"
                    push!(list_available_shippos, visit)                  
                end
            end
            pos_chosen = list_available_shippos
        end

        if tactic=="count"
            alpha_value = AllBoatsCount
            list_available_shippos = Vector{NewVisit}()
            for visit in new_visits
                if count_available[visit.n][visit.c]<=min_count+alpha_value*(max_count-min_count)
                    visit.store.tacticAll="count"
                    push!(list_available_shippos, visit)                  
                end
            end
            pos_chosen = list_available_shippos
        end

        if tactic=="time"
            alpha_value = AllBoatsTime
            list_available_shippos = Vector{NewVisit}()
            for visit in new_visits
                if visit.time<=min_time+alpha_value*(max_time-min_time)
                    visit.store.tacticAll="time"
                    push!(list_available_shippos, visit)                  
                end
            end
            pos_chosen = list_available_shippos
        end
    end

    pos_chosen_constrained=[]
    #if length(new_visits_constrained)>0
    #    tactic=paramchosen.TacticAllBoats
    #    if tactic=="cost"
    #        alpha_value = AllBoatsCost
    #        list_available_shippos = Vector{NewVisit}()
    #        for visit in new_visits_constrained
    #            if visit.cost<=min_cost_const+alpha_value*(max_cost_const-min_cost_const)
    #                visit.store.tacticAll="cost"
    #                push!(list_available_shippos, visit)                  
    #            end
    #        end
    #        pos_chosen_constrained = list_available_shippos
    #    end

    #    if tactic=="count"
    #        alpha_value = AllBoatsCount
    #        list_available_shippos = Vector{NewVisit}()
    #        for visit in new_visits_constrained
    #            if count_available[visit.n][visit.c]<=min_count_const+alpha_value*(max_count_const-min_count_const)
    #                visit.store.tacticAll="count"
    #                push!(list_available_shippos, visit)                  
    #            end
    #        end
    #        pos_chosen_constrained = list_available_shippos
    #    end

    #    if tactic=="time"
    #        alpha_value = AllBoatsTime
    #        list_available_shippos = Vector{NewVisit}()
    #        for visit in new_visits_constrained
    #            if visit.time<=min_time_const+alpha_value*(max_time_const-min_time_const)
    #                visit.store.tacticAll="time"
    #                push!(list_available_shippos, visit)                  
    #            end
    #        end
    #        pos_chosen_constrained = list_available_shippos
    #    end
    #end

    if length(pos_chosen)+length(pos_chosen_constrained)>0
        #alpha_value = Alpha.RateConstrained[paramchosen.IndexRateConstrained]
        #epsilon=0.0000001
        #if length(pos_chosen)>1
        #    pos_chosen = sample(pos_chosen, ceil(Int,alpha_value*length(pos_chosen)+epsilon))
        #end
        if length(pos_chosen_constrained)>1
            return pos_chosen_constrained[rand(1:length(pos_chosen_constrained))], true
        else
            return pos_chosen[rand(1:length(pos_chosen))], true
        end
    else
        return NewVisit(0,0,0,0,0.0,0.0,0.0, false, ToStoreVisit(SplitCosts(0,0,0,0,0,0),0,"","")), false
    end
end


function greedyrandomizedconstruction(inst::Instance, paramchosen::ChosenParameters, paramfixed::FixedParameters, max_time)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT = inst
    sol = initializeSol(inst)
    new_visit, feasible = SelectNewVisitAllShips(inst, sol, paramchosen, paramfixed)
    continue_=true
    start = time_ns()
    elapsed = round((time_ns()-start)/1e9,digits=3)
    count_when=0
    while continue_ && elapsed<max_time
        @unpack n,c,b,t,cost,distance,constrained,store= new_visit
        count_when+=1
        store.when=count_when
        l = ceil(Int, shipsIn[n].l/qli)
        hand = ceil(Int, h[n][c][b])
        @unpack M, visits = sol
        sol.visits = updateTimesAfterVisit(inst, visits, n, c, b, t)
        sol.M = updateMpositions(inst, n, c, b, t, l, hand, M)
        @unpack M, visits = sol
        sol.visits[n][c].p = Pi[n][c]
        sol.visits[n][c].b = b
        sol.visits[n][c].t = t
        sol.visits[n][c].planned = true
        sol.visits[n][c].store = store
                
        continue_=false
        for boat_vis in sol.visits
            for vis in boat_vis
                if vis.planned == false
                    continue_ = true
                end
            end
        end


        if continue_
            new_visit, feasible = SelectNewVisitAllShips(inst, sol, paramchosen, paramfixed)
            if feasible==false
                #print('\n')
                #print("Continue but not feasible")
                #print('\n')
                #print(round((time_ns()-start)/1e9,digits=3))
                sol.failed=1
                nb_to_remove_visits=0
                nb_to_remove_ships=0
                #print('\n')
                #print("###################")
                #print('\n')
                #print("We did not find a at first try solution")
                #print('\n')
                #print(sol.visits)
                #print('\n')
                for n in 1:N
                    this_ship_remove = false
                    for (c,p) in enumerate(inst.Pi[n])
                        if sol.visits[n][c].planned == false
                            sol.visits[n][c].failed=1
                            this_ship_remove = true
                            if c<length(inst.Pi[n])
                                for i in c:length(inst.Pi[n])
                                    nb_to_remove_visits+=1
                                    sol.visits[n][i].p = -1
                                    sol.visits[n][i].b = -1
                                    sol.visits[n][i].t = -1
                                    sol.visits[n][i].planned = false
                                    sol.visits[n][i].minT = max(shipsIn[n].sT[i], T[n,i,1])
                                    sol.visits[n][i].maxT = maxT
                                end
                            end
                        end
                    end
                    if this_ship_remove
                        nb_to_remove_ships+=1
                    end
                end
                #print('\n')
                #print("Number of visits removed :")
                #print('\n')
                #print(nb_to_remove_visits)
                #print('\n')
                #print("Number of ships removed :")
                #print('\n')
                #print(nb_to_remove_ships)
                sol.M = generateOccupiedMx(inst, sol.visits)
                boat_order = 1:N
                random_boat_order = shuffle(boat_order)
                for n in random_boat_order
                    for (c,p) in enumerate(inst.Pi[n])
                        if sol.visits[n][c].planned == false
                            for i in c:length(inst.Pi[n])
                                new_visit, new_visit_constrained, feasible=SelectNewVisitPerShipReconstruct(inst, sol, paramchosen, paramfixed, n, i, Pi[n][i], false)
                                if feasible
                                    n_ = new_visit.n
                                    c_ = new_visit.c
                                    b_ = new_visit.b
                                    t_ = new_visit.t
                                    l_ = ceil(Int, shipsIn[n_].l/qli)
                                    store = new_visit.store
                                    store.when=count_when
                                    store.tacticAll = "reconstruct"
                                    hand = ceil(Int, h[n_][c_][b_])
                                    @unpack M, visits = sol
                                    sol.visits =updateTimesAfterVisit(inst, visits, n_, c_, b_, t_)
                                    sol.M = updateMpositions(inst, n_, c_, b_, t_, l_, hand, M)
                                    sol.visits[n_][c_].p = Pi[n][c_]
                                    sol.visits[n_][c_].b = b_
                                    sol.visits[n_][c_].t = t_
                                    sol.visits[n_][c_].planned = true
                                    sol.visits[n_][c_].store = store
                                else
                                    elapsed = round((time_ns()-start)/1e9,digits=3)
                                    #print('\n')
                                    #print("Huston we have a problem : $elapsed")
                                    #print('\n')
                                    #print("The problematic visit : $n, $i")
                                    #print('\n')
                                    #print("###################")
                                    #print('\n')
                                    #print(sol.visits[n])
                                    #print('\n')
                                    return initializeSol(inst)
                                end
                            end
                        end
                    end
                end
                continue_=false
                for boat_vis in sol.visits
                    for vis in boat_vis
                        if vis.planned == false
                            continue_ = true
                        end
                    end
                end
                if continue_==false
                    elapsed = round((time_ns()-start)/1e9,digits=3)
                    #print('\n')
                    #print("The time it took : $elapsed")
                end
            end
        else
            elapsed = round((time_ns()-start)/1e9,digits=3)
            #print('\n')
            #print("The time it took : $elapsed")
        end
    end

    #print("###################")
    #print('\n')
    #print("Number of time we had to stop the heuristic")
    #print('\n')
    #print(not_find_total)
    #print('\n')
    return sol #, solution
end


function ChooseParam(allparam::AllParameters, type1, type2, type3)
    @unpack Proba, AverageCost, Nbexp, Q = allparam
    if type1=="all"
        tacticoneship = sample(["cost","dist","time"], Weights(Proba.TacticOneBoat))
    end
    if type1=="cost"
        tacticoneship="cost"
    end
    if type1=="dist"
        tacticoneship="dist"
    end
    if type1=="time"
        tacticoneship="time"
    end

    if type2=="all"
        tacticallship = sample(["cost","count","time"], Weights(Proba.TacticAllBoats))
    end
    if type2=="cost"
        tacticallship="cost"
    end
    if type2=="count"
        tacticallship="count"
    end
    if type2=="time"
        tacticallship="time"
    end

    if type3=="all"
        tacticlocal = sample(["random","boat"], Weights(Proba.TacticLocalSearch))
    end
    if type3=="random"
        tacticlocal="random"
    end
    if type3=="boat"
        tacticlocal="boat"
    end
    return ChosenParameters(tacticoneship,tacticallship,tacticlocal)
end


function UpdateParameters(paramchosen::ChosenParameters, allparam::AllParameters, cost, new_cost)
    newparam=deepcopy(allparam)

    index=0
    if paramchosen.TacticOneBoat=="cost"
        index=1
    end
    if paramchosen.TacticOneBoat=="dist"
        index=2
    end
    if paramchosen.TacticOneBoat=="time"
        index=3
    end
    nbexp = allparam.Nbexp.TacticOneBoat[index]
    newparam.AverageCost.TacticOneBoat[index] = (allparam.AverageCost.TacticOneBoat[1]*nbexp+new_cost)/(nbexp+1)
    newparam.Nbexp.TacticOneBoat[index]=nbexp+1
    for i in 1:3
        newav = deepcopy(newparam.AverageCost.TacticOneBoat[i])
        newparam.Q.TacticOneBoat[i]=cost/newav
    end
    for i in 1:3
        newq=deepcopy(newparam.Q.TacticOneBoat[i])
        newparam.Proba.TacticOneBoat[i] = newq/sum(newparam.Q.TacticOneBoat)
    end

    index=0
    if paramchosen.TacticAllBoats=="cost"
        index=1
    end
    if paramchosen.TacticAllBoats=="count"
        index=2
    end
    if paramchosen.TacticAllBoats=="time"
        index=3
    end
    nbexp = allparam.Nbexp.TacticAllBoats[index]
    newparam.AverageCost.TacticAllBoats[index] = (allparam.AverageCost.TacticAllBoats[1]*nbexp+new_cost)/(nbexp+1)
    newparam.Nbexp.TacticAllBoats[index]=nbexp+1
    for i in 1:3
        newav = deepcopy(newparam.AverageCost.TacticAllBoats[i])
        newparam.Q.TacticAllBoats[i]=cost/newav
    end
    for i in 1:3
        newq=deepcopy(newparam.Q.TacticAllBoats[i])
        newparam.Proba.TacticAllBoats[i] = newq/sum(newparam.Q.TacticAllBoats)
    end

    index=0
    if paramchosen.TacticLocalSearch=="random"
        index=1
    end
    if paramchosen.TacticLocalSearch=="boat"
        index=2
    end
    nbexp = allparam.Nbexp.TacticLocalSearch[index]
    newparam.AverageCost.TacticLocalSearch[index] = (allparam.AverageCost.TacticLocalSearch[1]*nbexp+new_cost)/(nbexp+1)
    newparam.Nbexp.TacticLocalSearch[index]=nbexp+1
    for i in 1:2
        newav = deepcopy(newparam.AverageCost.TacticLocalSearch[i])
        newparam.Q.TacticLocalSearch[i]=cost/newav
    end
    for i in 1:2
        newq=deepcopy(newparam.Q.TacticLocalSearch[i])
        newparam.Proba.TacticLocalSearch[i] = newq/sum(newparam.Q.TacticLocalSearch)
    end
    return newparam
end


function GRASP_reactive(seed,N,Nout,qli, type1, type2, type3, paramfixed, time_local, max_time_heur, max_time, expname, location)
    inst = readInstFromFile(location*"MCBAP-multi-port-berth-allocation-problem/data_small/CP2_Inst_$seed"*"_$N"*"_$Nout"*"_$qli"*".txt")
    dictionnary_csv =
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp = inst
    cost=1000000000
    worst_cost=1000000000
    sol=initializeSol(inst)
    allparam = initializeParam(inst)
    start = time_ns()
    elapsed = round((time_ns()-start)/1e9,digits=3)
    start_iter=time_ns()
    elapsed_iter = round((time_ns()-start_iter)/1e9,digits=3)
    nb_iter=0
    d_alliter_before=Dict()
    d_alliter_after=Dict()
    penalty_cost_list = Vector{Tuple}()
    min_cost_heur=1000000000
    while elapsed<max_time
        paramchosen = ChooseParam(allparam, type1, type2, type3)
        start_heur = time_ns()
        new_sol = greedyrandomizedconstruction(inst, paramchosen, paramfixed, max_time_heur)

        elapsed_heur = round((time_ns()-start_heur)/1e9,digits=3)
        feasible = true
        ## No conflicts with within one boat schedule :
        for n in 1:N
            # The times :
            for (c,p) in enumerate(inst.Pi[n])
                if new_sol.visits[n][c].planned == false
                    feasible = false              
                end
            end
        end
        if feasible && checkSolutionFeasability(inst, new_sol)
            new_cost_heur, delay_cost_heur, waiting_cost_heur, penalty_cost_heur, handling_cost_heur, fuel_cost_heur = checkSolutionCost(inst, new_sol)
            
            #print('\n')
            #print("Before local search")
            #print('\n')
            #print(new_sol.visits)
            d_before=prepareSolIter(seed,N,Nout,qli,nb_iter,inst, new_sol, new_cost_heur, paramfixed, expname)
            d_alliter_before[nb_iter]=d_before
            start_local=time_ns()

            new_sol, new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost = pushTime(inst, new_sol, time_local)
            if new_cost<min_cost_heur+0.1*min_cost_heur
                new_sol, new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost = local_search(inst, deepcopy(new_sol), ceil(Int,cost), paramchosen, paramfixed, time_local)
                if new_cost<min_cost_heur
                    min_cost_heur=new_cost
                end
            end

            #new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost = new_cost_heur, delay_cost_heur, waiting_cost_heur, penalty_cost_heur, handling_cost_heur, fuel_cost_heur

            elapsed_local = round((time_ns()-start_local)/1e9,digits=3)
            #print('\n')
            #print("After local search")
            #print('\n')
            #print(new_sol.visits)
            #print('\n')
            #print("Cost at the end of local search")
            #print('\n')
            #print(new_cost)
            #print('\n')
            #print("Old cost")
            #print('\n')
            #print(new_cost_heur)
            new_sol.store.costHeur=SplitCosts(ceil(Int,new_cost_heur), ceil(Int,delay_cost_heur), ceil(Int,waiting_cost_heur), ceil(Int,penalty_cost_heur), ceil(Int,handling_cost_heur), ceil(Int,fuel_cost_heur))
            new_sol.store.costLocal=SplitCosts(ceil(Int,new_cost), ceil(Int,delay_cost), ceil(Int,waiting_cost), ceil(Int,penalty_cost), ceil(Int, handling_cost), ceil(Int,fuel_cost))
            new_sol.store.timeHeur=elapsed_heur
            new_sol.store.timeLocalSearch=elapsed_local
            new_sol.store.parameters=allparam.Proba
            elapsed_iter = round((time_ns()-start_iter)/1e9,digits=3)
            penalty_cost_list = Vector{Tuple}()
            for n in 1:N
                for (c,p) in enumerate(Pi[n])
                    t = new_sol.visits[n][c].t
                    b = new_sol.visits[n][c].b
                    this_cost, delay_cost, waiting_cost, penalty, handling_cost, fuel_cost, feas = computeCostPosSol(inst, n, c, b, t,new_sol)
                    push!(penalty_cost_list, (n,c,penalty))
                    new_sol.visits[n][c].store.cost =  SplitCosts(ceil(Int, this_cost), ceil(Int,delay_cost), ceil(Int,waiting_cost), ceil(Int,penalty), ceil(Int,handling_cost), ceil(Int,fuel_cost))
                end
            end

            d_after = prepareSolIter(seed,N,Nout,qli,nb_iter,inst, new_sol, new_cost, paramfixed, expname)
            
            if cost==1000000000
                worst_cost=deepcopy(new_cost)
            end
            if new_cost<cost
                #print('\n')
                #print("################# lol")
                d_after["better"]=1
                cost=deepcopy(new_cost)
                sol=deepcopy(new_sol)
            end
            if new_cost>=worst_cost
                worst_cost=deepcopy(new_cost)
            end

            d_alliter_after[nb_iter]=d_after
            nb_iter+=1

            allparam = UpdateParameters(paramchosen, allparam, cost, new_cost)    
            #allparam = UpdateParametersExp(paramchosen, allparam, N, cost, new_cost)                
        else
            allparam = UpdateParameters(paramchosen, allparam, cost, worst_cost) 
            #allparam = UpdateParametersExp(paramchosen, allparam, N, cost, worst_cost) 
        end
        elapsed = round((time_ns()-start)/1e9,digits=3)
    end
    #print('\n')
    #print("Cost at the end of heur")
    #print('\n')
    #print(cost)
    ## No conflicts with within one boat schedule :
    feasible=true
    for n in 1:N
        # The times :
        for (c,p) in enumerate(inst.Pi[n])
            if sol.visits[n][c].planned == false
                feasible = false              
            end
        end
    end
    if feasible && checkSolutionFeasability(inst, sol)
        print('\n')
        print("The solution is feasible")
        for iter in keys(d_alliter_after)
            CSV.write(location*"results_jobs/benchmarks_HEUR/noCplexGRASP/$expname"*"/iterations_before_local/sol_$seed"*"_$N"*"_$Nout"*"_$qli"*"/iter_$iter"*"_beforelocal"*".csv", d_alliter_before[iter])
            CSV.write(location*"results_jobs/benchmarks_HEUR/noCplexGRASP/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli"*"/iter_$iter"*".csv", d_alliter_after[iter])
        end
        print('\n')
        print("The number of iterations is :")
        print('\n')
        print(nb_iter)
        nb_iter+=1
        d = prepareSolIter(seed,N,Nout,qli,nb_iter,inst, sol, cost, paramfixed, expname)
        CSV.write(location*"results_jobs/benchmarks_HEUR/noCplexGRASP/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli"*"/iter_$nb_iter"*".csv", d)
        return sol, cost, allparam
    else
        cost=1000000000
        sol=initializeSol(inst)
        allparam = initializeParam(inst)
        return sol, cost, allparam
    end
end


function testallfunction()
    location = "D:/DTU-Courses/DTU-Thesis/berth_allocation/"
    #location="/zhome/c3/6/164957/code_git/"


    # The experience name :
    expname="test"

    # The tactic types :
    type1="time" 
    type2="cost" 
    type3="random"

    # The alpha parameters for each tactic :
    prop_oneboatcost = 0.1
    prop_oneboatdist = 0.2
    prop_oneboattime = 0.1
    prop_allboatcost = 0.2
    prop_allboatcount = 0.2
    prop_allboattime = 0.5

    # The number of boat to remove for the local search :
    alphaboat=5
    alpharandom=15

    # All the parameters :
    paramfixed=FixedParameters(prop_oneboatcost,prop_oneboatdist,prop_oneboattime,prop_allboatcost,prop_allboatcount,prop_allboattime,alpharandom,alphaboat)

    # Maximum time for the local search :
    time_local=5

    # Maximum time for the heuristic :
    max_time_heur=30

    # maximum time for the experiment :
    max_time=300


    #print("Start")
    #print('\n')
    for N in 11:11
        for qli in [10]
            for Nout in 3:3
                for seed in 1:1
                    #print("The instance : $seed"*"_$N"*"_$Nout"*"_$qli")
                    inst = readInstFromFile(location*"MCBAP-multi-port-berth-allocation-problem/data_small/CP2_Inst_$seed"*"_$N"*"_$Nout"*"_$qli"*".txt")
                    sol, cost, allparam = GRASP_reactive(seed,N,Nout,qli, type1, type2, type3, paramfixed, time_local, max_time_heur, max_time, expname, location)
                    feasible=true
                    for n in 1:N
                        # The times :
                        for (c,p) in enumerate(inst.Pi[n])
                            if sol.visits[n][c].planned == false
                                feasible = false              
                            end
                        end
                    end
                    if feasible && checkSolutionFeasability(inst, sol)
                        print('\n')
                        print("The probas :")
                        print('\n')
                        print(allparam.Proba)
                        print('\n')
                        print("The solution :")
                        print('\n')
                        print(sol.visits)
                        print('\n')
                        print("And the cost is ")
                        print('\n')
                        print(cost)
                    end
                end
            end
        end
    end
end

