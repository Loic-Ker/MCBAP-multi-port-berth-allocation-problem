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
include("../localCPLEX.jl")

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
    #        if sol.M[n][c][b,t+1]
    #            this_cost, delay_cost, waiting_cost, penalty, handling_cost, fuel_cost, feas = computeCostPosSol(inst, n, c, b, t ,sol)
    #            if feas
    #                distance = abs(bestPos-b)/l
    #                if distance>=max_dist_const
    #                    max_dist_const=distance
    #                end
    #                if this_cost>=max_cost_const
    #                    max_cost_const=this_cost
    #                end
    #                if distance<=min_dist_const
    #                    min_dist_const=distance
    #                end
    #                if this_cost<=min_cost_const
    #                    min_cost_const=this_cost
    #                end
    #                push!(new_visits_constrained, NewVisit(n,c,b,t,this_cost,distance, true, ToStoreVisit(SplitCosts(ceil(Int, this_cost), ceil(Int,delay_cost), ceil(Int,waiting_cost), ceil(Int,penalty), ceil(Int,handling_cost), ceil(Int,fuel_cost)),0,"","")))
    #            end
    #        end
    #    end
    #end
    if printdetails
        print("lol")
        print('\n')
        print(new_visits)
        print('\n')
        print(new_visits_constrained)
        print('\n')
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
    
    #if length(new_visits_constrained)>0
    #    if tactic=="cost"
    #        alpha_value= Alpha.CostTypeShip[paramchosen.IndexTypeShip]
    #        list_available_shippos = Vector{NewVisit}()
    #        for visit in new_visits_constrained
    #            if visit.constrained
    #                if visit.cost<=min_cost_const+alpha_value*(max_cost_const-min_cost_const)
    #                    visit.store.tacticBoat="cost"
    #                    push!(list_available_shippos, visit)
    #                end                    
    #            end
    #        end
    #        pos_chosen_constrained = list_available_shippos[rand(1:length(list_available_shippos))]
    #    end
    #    if tactic=="dist"
    #        alpha_value= Alpha.DistanceTypeShip[paramchosen.IndexTypeShip]
    #        list_available_shippos = Vector{NewVisit}()
    #        for visit in new_visits_constrained
    #            if visit.constrained
    #                if visit.distance<=min_dist_const+alpha_value*(max_dist_const-min_dist_const)
    #                    visit.store.tacticBoat="dist"
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
    #max_count_const=0
    #max_cost_const=0
    #min_count_const=1000000000
    #min_cost_const=1000000000
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

    #pos_chosen_constrained=[]
    #if length(new_visits_constrained)>0
    #    tactic=paramchosen.TacticAllShip
    #    if tactic=="cost"
    #        alpha_value = Alpha.CostAllShip[paramchosen.IndexAllShip]
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
    #        alpha_value = Alpha.CountAllShip[paramchosen.IndexAllShip]
    #        list_available_shippos = Vector{NewVisit}()
    #        for visit in new_visits_constrained
    #            if count_available[visit.n][visit.c]<=min_count_const+alpha_value*(max_count_const-min_count_const)
    #                visit.store.tacticAll="count"
    #                push!(list_available_shippos, visit)                  
    #            end
    #        end
    #        pos_chosen_constrained = list_available_shippos
    #    end
    #end

    if length(pos_chosen)>0#+length(pos_chosen_constrained)>0
        #alpha_value = Alpha.RateConstrained[paramchosen.IndexRateConstrained]
        #epsilon=0.0000001
        #if length(pos_chosen)>1
        #    pos_chosen = sample(pos_chosen, ceil(Int,alpha_value*length(pos_chosen)+epsilon))
        #end
        #if length(pos_chosen_constrained)>1
        #    pos_chosen_constrained = sample(pos_chosen_constrained, ceil(Int,(1-alpha_value)*length(pos_chosen_constrained)+epsilon))
        #end
        #all_pos=vcat(pos_chosen,pos_chosen_constrained)
        #return all_pos[rand(1:length(all_pos))], true
        return pos_chosen[rand(1:length(pos_chosen))], true
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
        if b==0
            print('\n')
            print(new_visit)
        end
        l = ceil(Int, shipsIn[n].l/qli)
        hand = ceil(Int, h[n][c][b])
        @unpack M, visits = sol
        sol.visits = updateTimesAfterVisit(inst, visits, n, c, b, t)
        sol.M = updateMpositions(inst, n, c, b, t, l, hand, M)
        @unpack M, visits = sol
        sol.M = updateMtimes(inst, visits, n, c, M) 
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
                sol.failed=1
                print('\n')
                print("###################")
                print('\n')
                print("We did not find a at first try solution")
                print('\n')
                print(sol.visits)
                print('\n')
                for n in 1:N
                    for (c,p) in enumerate(inst.Pi[n])
                        if sol.visits[n][c].planned == false
                            sol.visits[n][c].failed=1
                            if c<length(inst.Pi[n])
                                for i in c:length(inst.Pi[n])
                                    sol.visits[n][i].p = -1
                                    sol.visits[n][i].b = -1
                                    sol.visits[n][i].t = -1
                                    sol.visits[n][i].planned = false
                                    sol.visits[n][i].minT = max(shipsIn[n].sT[i], T[n,i,1])
                                    sol.visits[n][i].maxT = min(5*T[n,i,2],maxT)
                                end
                                sol.M = generateOccupiedMx(inst, sol.visits)
                                sol = updateTimesInitialization(inst,sol)
                                for i in c:length(inst.Pi[n])
                                    new_visit, new_visit_constrained, feasible=SelectNewVisitPerShip(inst, sol, paramchosen, paramfixed, n, i, Pi[n][i], false)
                                    all_visits=Vector{NewVisit}()
                                    if new_visit.n !=0
                                        push!(all_visits, new_visit)
                                    end
                                    if new_visit_constrained.n !=0
                                        push!(all_visits, new_visit_constrained)
                                    end                                
                                    if feasible
                                        new_visit=all_visits[rand(1:length(all_visits))]
                                        count_when+=1
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
                                        @unpack M, visits = sol
                                        sol.M = updateMtimes(inst, visits, n_, c_, M) 
                                        sol.visits[n_][c_].p = Pi[n][c_]
                                        sol.visits[n_][c_].b = b_
                                        sol.visits[n_][c_].t = t_
                                        sol.visits[n_][c_].planned = true
                                        sol.visits[n_][c_].store = store
                                    else
                                        return initializeSol(inst)
                                    end
                                end
                            else
                                return initializeSol(inst)
                            end
                        end
                    end
                end   
                continue_=false
                elapsed = round((time_ns()-start)/1e9,digits=3)
                print('\n')
                print("The time it took : $elapsed")
            end
        else
            elapsed = round((time_ns()-start)/1e9,digits=3)
            print('\n')
            print("The time it took : $elapsed")
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
            print('\n')
            print("Before local search")
            print('\n')
            print(new_sol.visits)

            d=prepareSolIter(seed,N,Nout,qli,nb_iter,inst, new_sol, new_cost_heur, paramfixed, expname)
            CSV.write(location*"results_jobs/benchmarks_HEUR/parametricGRASP/$expname"*"/iterations_before_local/sol_$seed"*"_$N"*"_$Nout"*"_$qli"*"/iter_$nb_iter"*"_beforelocal"*".csv", d)
            
            start_local=time_ns()
            new_sol, new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost = local_search(inst, deepcopy(new_sol), ceil(Int,cost), paramchosen, paramfixed, time_local)
            elapsed_local = round((time_ns()-start_local)/1e9,digits=3)
            print('\n')
            print("After local search")
            print('\n')
            print(new_sol.visits)
            print('\n')
            print("Cost at the end of local search")
            print('\n')
            print(new_cost)
            print('\n')
            print("Old cost")
            print('\n')
            print(new_cost_heur)
            new_sol.store.costHeur=SplitCosts(ceil(Int,new_cost_heur), ceil(Int,delay_cost_heur), ceil(Int,waiting_cost_heur), ceil(Int,penalty_cost_heur), ceil(Int,handling_cost_heur), ceil(Int,fuel_cost_heur))
            new_sol.store.costLocal=SplitCosts(ceil(Int,new_cost), ceil(Int,delay_cost), ceil(Int,waiting_cost), ceil(Int,penalty_cost), ceil(Int, handling_cost), ceil(Int,fuel_cost))
            new_sol.store.timeHeur=elapsed_heur
            new_sol.store.timeLocalSearch=elapsed_local
            new_sol.store.parameters=allparam.Proba
            elapsed_iter = round((time_ns()-start_iter)/1e9,digits=3)

            d = prepareSolIter(seed,N,Nout,qli,nb_iter,inst, new_sol, cost, paramfixed, expname)
            #CSV.write(location*"results_jobs/benchmarks_HEUR/parametricGRASP/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli"*"/iter_$nb_iter"*".csv", d)
            nb_iter+=1
            
            if elapsed_iter>max_time/50
                start_iter=time_ns()
                elapsed_iter = round((time_ns()-start_iter)/1e9,digits=3)
            end
            if cost==1000000000
                worst_cost=deepcopy(new_cost)
            end
            if new_cost<cost
                print('\n')
                print("################# lol")
                cost=deepcopy(new_cost)
                sol=deepcopy(new_sol)
            end
            if new_cost>=worst_cost
                worst_cost=deepcopy(new_cost)
            end
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
        d = prepareSolIter(seed,N,Nout,qli,nb_iter,inst, sol, cost, paramfixed, expname)
        nb_iter+=1
        CSV.write(location*"results_jobs/benchmarks_HEUR/parametricGRASP/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli"*"/iter_$nb_iter"*".csv", d)
        return sol, cost, allparam
    else
        cost=1000000000
        sol=initializeSol(inst)
        allparam = initializeParam(inst)
        return sol, cost, allparam
    end
end

#inst = readInstFromFile("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_1_10_5_80.txt")
#GRASP_reactive(inst::Instance, type1, type2, adjustproba, alphaboat, alpharandom, time_local, max_time_heur, max_time, expname)
#CSV.write("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/reactiveGRASP/sols/HEUR_LOCAL_exp1_sol_$seed"*"_$N"*"_$Nout"*"_$qli"*".csv", d)
#expname="exp2"
#seed=1
#N=10
#Nout=5
#qli=80
#if isdir("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/parametricGRASP/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli")==false
#    mkdir("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/parametricGRASP/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli")
#end
#sol, cost, allparam = GRASP_reactive(seed,N,Nout,qli,"both","both", AdjustProba(3,3,3,3,3,3), 5, 18, 20, 25, 50, "exp2")
#print('\n')
#print("############")
#print('\n')
#print(cost)
#print('\n')
#print(sol.visits)
#print('\n')
#print(sol.store)
#print(cost)
#print('\n')
#print("Proba tactic type :")
#print('\n')
#print(allparam.Proba.ChooseTacticLocalSearch)
#print('\n')
#print("Proba distance :")
#print('\n')
#print(allparam.Proba.DistanceTypeShip)
#print('\n')
#print("Proba cost type :")
#print('\n')
#print(allparam.Proba.CostTypeShip)
#print('\n')
#print("Proba tactic all :")
#print('\n')
#print(allparam.Proba.ChooseTacticAllShip)
#print('\n')
#print("Proba count all :")
#print('\n')
#print(allparam.Proba.CountAllShip)
#print('\n')
#print("Proba cost all :")
#print('\n')
#print(allparam.Proba.CostAllShip)
#print("Rate of constrained all :")
#print('\n')
#print(allparam.Proba.RateConstrained)
