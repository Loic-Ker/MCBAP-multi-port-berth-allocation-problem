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
include("../localsearch2.jl")

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


################ Look at one visit :
# Find constrained positions
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



# compute feasible positions (b,t) for one visit : 
function SelectNewVisitPerShip(inst::Instance, sol::Sol, paramchosen::ChosenParameters, allparam::AllParameters, paramfixed::FixedParameters, n, c, p, printdetails) 
    @unpack N, h, qli, shipsIn, Pi, Bp, delta, dist, T, maxC, maxT, Pc = inst
    @unpack M = sol
    @unpack Alpha, Proba, AverageCost, Nbexp, Q = allparam
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
        find_it=false
        t=t1-1
        while t<t2 && find_it==false
            t=t+1
            if t + hand <= t2
                if M[n][c][b,t+1]
                    this_cost, delay_cost, waiting_cost, penalty, handling_cost, fuel_cost, feas = computeCostPosSol(inst, n, c, b, t,sol)
                    time=deepcopy(t)
                    #time=penalty
                    ## look at the the time window or only the penalty
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
                        if (max_time-min_time)/min_time>paramfixed.WindowSize
                            find_it=true
                        end
                        push!(new_visits, NewVisit(n,c,b,t,this_cost,distance,time, false, ToStoreVisit(SplitCosts(ceil(Int, this_cost), ceil(Int,delay_cost), ceil(Int,waiting_cost), ceil(Int,penalty), ceil(Int,handling_cost), ceil(Int,fuel_cost)),0,"","")))
                    end
                end
            end
        end
    end  
    if paramfixed.lookforconstrained==true
        ## Check the constrained ones
        pos = findConstrainedPos(inst, sol, p, t1, t2, n, c, l)
        for (b,t) in pos
            hand = ceil(Int, h[n][c][b])
            if t + hand <= t2
                if M[n][c][b,t+1]
                    this_cost, delay_cost, waiting_cost, penalty, handling_cost, fuel_cost, feas = computeCostPosSol(inst, n, c, b, t,sol)
                    time=penalty
                    if feas
                        distance = abs(bestPos-b)/l
                        if distance>=max_dist_const
                            max_dist_const=distance
                        end
                        if this_cost>=max_cost_const
                            max_cost_const=this_cost
                        end
                        if time>=max_time_const
                            max_time_const=time
                        end
                        if distance<=min_dist_const
                            min_dist_const=distance
                        end
                        if this_cost<=min_cost_const
                            min_cost_const=this_cost
                        end
                        if time<=min_time_const
                            min_time_const=time
                        end
                        push!(new_visits_constrained, NewVisit(n,c,b,t,this_cost,distance,time, true, ToStoreVisit(SplitCosts(ceil(Int, this_cost), ceil(Int,delay_cost), ceil(Int,waiting_cost), ceil(Int,penalty), ceil(Int,handling_cost), ceil(Int,fuel_cost)),0,"","")))
                    end
                end
            end
        end
    end
    pos_chosen=NewVisit(0,0,0,0,0.0,0.0,0.0, false, ToStoreVisit(SplitCosts(0,0,0,0,0,0),0,"",""))
    pos_chosen_constrained=NewVisit(0,0,0,0,0.0,0.0,0.0, false, ToStoreVisit(SplitCosts(0,0,0,0,0,0),0,"",""))
    tactic=paramchosen.TacticOneBoat
    if length(new_visits)>0
        if tactic=="cost"
            alpha_value= Alpha.CostOneShip[paramchosen.IndexOneShip]
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
            alpha_value= Alpha.DistOneShip[paramchosen.IndexOneShip]
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
            alpha_value= Alpha.TimeOneShip[paramchosen.IndexOneShip]
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
    
    if paramfixed.lookforconstrained==true
        if length(new_visits_constrained)>0
            if tactic=="cost"
                alpha_value= Alpha.CostOneShip[paramchosen.IndexOneShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits_constrained
                    if visit.constrained==true
                        if visit.cost<=min_cost_const+alpha_value*(max_cost_const-min_cost_const)
                            visit.store.tacticBoat="cost"
                            push!(list_available_shippos, visit)
                        end
                    end
                end
                pos_chosen_constrained = list_available_shippos[rand(1:length(list_available_shippos))]
            end
            if tactic=="dist"
                alpha_value= Alpha.DistOneShip[paramchosen.IndexOneShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits_constrained
                    if visit.constrained==true
                        if visit.distance<=min_dist_const+alpha_value*(max_dist_const-min_dist_const)
                            visit.store.tacticBoat="dist"
                            push!(list_available_shippos, visit)
                        end                   
                    end
                end
                pos_chosen_constrained = list_available_shippos[rand(1:length(list_available_shippos))]
            end
            if tactic=="time"
                alpha_value= Alpha.TimeOneShip[paramchosen.IndexOneShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits_constrained
                    if visit.constrained==true
                        if visit.time<=min_time_const+alpha_value*(max_time_const-min_time_const)
                            visit.store.tacticBoat="time"
                            push!(list_available_shippos, visit)
                        end                   
                    end
                end
                pos_chosen_constrained = list_available_shippos[rand(1:length(list_available_shippos))]
            end
        end
    end
    
    if length(new_visits)+length(new_visits_constrained)>0
        return pos_chosen, pos_chosen_constrained, true
    else
        return pos_chosen, pos_chosen_constrained, false
    end
end

function SelectNewVisitAllShips(inst::Instance, sol::Sol, paramchosen::ChosenParameters,allparam::AllParameters, paramfixed::FixedParameters) #, placed::Vector{Tuple{Int64,Int64}})
    #TODO: t1 cannot be greater than t2!
    # adapt to given planned previous visits for each ship
    @unpack N, h, qli, shipsIn, Pi, Bp, delta, dist, T, maxC, maxT = inst
    @unpack M = sol
    @unpack Alpha, Proba, AverageCost, Nbexp, Q = allparam
    # C = maximum(length.(Pi))
    new_visits = Vector{NewVisit}()
    new_visits_constrained = Vector{NewVisit}()
    max_cost=0
    max_time = 0
    max_dist=0
    min_cost=1000000000
    min_time=1000000000
    min_dist=1000000000
    max_cost_const=0
    max_time_const = 0
    max_dist_const=0
    min_cost_const=1000000000
    min_time_const=1000000000
    min_dist_const=1000000000
    for n in 1:N
        for (c,p) in enumerate(Pi[n])
            sch = sol.visits[n]
            if sch[c].planned == false
                visit, visit_constrained, feasible=SelectNewVisitPerShip(inst, sol, paramchosen, allparam, paramfixed, n, c, p, false)
                if feasible
                    if visit.n!=0
                        this_cost = visit.cost
                        this_time = visit.time
                        this_dist = visit.distance
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
                        if this_dist>=max_dist
                            max_dist=this_dist
                        end
                        if this_dist<=min_dist
                            min_dist_const=this_dist
                        end
                        push!(new_visits,visit)
                    end
                    if visit_constrained.n!=0
                        this_cost = visit_constrained.cost
                        this_time = visit_constrained.time
                        this_dist = visit_constrained.distance
                        if this_cost>=max_cost_const
                            max_cost_const=this_cost
                        end
                        if this_cost<=min_cost_const
                            min_cost_const=this_cost
                        end
                        if this_time>max_time_const
                            max_time_const=this_time
                        end
                        if this_time<min_time_const
                            min_time_const=this_time
                        end
                        if this_dist>=max_dist_const
                            max_dist_const=this_dist
                        end
                        if this_dist<=min_dist_const
                            min_dist_const=this_dist
                        end
                        push!(new_visits_constrained,visit_constrained)
                    end
                end
            end
        end
    end

    pos_chosen=[]
    if paramchosen.Reversed=="no"
        if length(new_visits)>0
            tactic=paramchosen.TacticAllBoats
            if tactic=="cost"
                alpha_value = Alpha.CostAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits
                    if visit.cost<=min_cost+alpha_value*(max_cost-min_cost)
                        visit.store.tacticAll="cost"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen = list_available_shippos
            end

            if tactic=="dist"
                alpha_value = Alpha.DistAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits
                    if visit.distance<=min_dist+alpha_value*(max_dist-min_dist)
                        visit.store.tacticAll="dist"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen = list_available_shippos
            end

            if tactic=="time"
                alpha_value = Alpha.TimeAllShip[paramchosen.IndexAllShip]
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
        if length(new_visits_constrained)>0
            tactic=paramchosen.TacticAllBoats
            if tactic=="cost"
                alpha_value = Alpha.CostAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits_constrained
                    if visit.cost<=min_cost_const+alpha_value*(max_cost_const-min_cost_const)
                        visit.store.tacticAll="cost"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen_constrained = list_available_shippos
            end

            if tactic=="dist"
                alpha_value = Alpha.DistAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits_constrained
                    if visit.distance<=min_dist_const+alpha_value*(max_dist_const-min_dist_const)
                        visit.store.tacticAll="dist"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen_constrained = list_available_shippos
            end

            if tactic=="time"
                alpha_value = Alpha.TimeAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits_constrained
                    if visit.time<=min_time_const+alpha_value*(max_time_const-min_time_const)
                        visit.store.tacticAll="time"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen_constrained = list_available_shippos
            end
        end
    else
        if length(new_visits)>0
            tactic=paramchosen.ReversedTacticAllBoats
            if tactic=="cost"
                alpha_value = Alpha.ReversedCostAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits
                    if visit.cost>=min_cost+alpha_value*(max_cost-min_cost)
                        visit.store.tacticAll="cost"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen = list_available_shippos
            end

            if tactic=="dist"
                alpha_value = Alpha.ReversedDistAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits
                    if visit.distance>=min_dist+alpha_value*(max_dist-min_dist)
                        visit.store.tacticAll="dist"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen = list_available_shippos
            end

            if tactic=="time"
                alpha_value = Alpha.ReversedTimeAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits
                    if visit.time>=min_time+alpha_value*(max_time-min_time)
                        visit.store.tacticAll="time"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen = list_available_shippos
            end
        end

        pos_chosen_constrained=[]
        if length(new_visits_constrained)>0
            tactic=paramchosen.ReversedTacticAllBoats
            if tactic=="cost"
                alpha_value = Alpha.ReversedCostAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits_constrained
                    if visit.cost>=min_cost_const+alpha_value*(max_cost_const-min_cost_const)
                        visit.store.tacticAll="cost"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen_constrained = list_available_shippos
            end

            if tactic=="dist"
                alpha_value = Alpha.ReversedDistAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits_constrained
                    if visit.distance>=min_dist_const+alpha_value*(max_dist_const-min_dist_const)
                        visit.store.tacticAll="dist"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen_constrained = list_available_shippos
            end

            if tactic=="time"
                alpha_value = Alpha.ReversedTimeAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits_constrained
                    if visit.time>=min_time_const+alpha_value*(max_time_const-min_time_const)
                        visit.store.tacticAll="time"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen_constrained = list_available_shippos
            end
        end
    end

    if length(pos_chosen)+length(pos_chosen_constrained)>0
        if length(pos_chosen_constrained)>=1
            random_value = rand()
            if random_value<Alpha.RateConstrained[paramchosen.IndexRateConstrained]
                return pos_chosen_constrained[rand(1:length(pos_chosen_constrained))], true
            else
                return pos_chosen[rand(1:length(pos_chosen))], true
            end
        else
            return pos_chosen[rand(1:length(pos_chosen))], true
        end
    else
        return NewVisit(0,0,0,0,0.0,0.0,0.0, false, ToStoreVisit(SplitCosts(0,0,0,0,0,0),0,"","")), false
    end
end


# compute feasible positions (b,t) in different ways :  
## Order by best costs
## Get random ones
## Get the ones with less feasible solutions
## Get the ones closer to other ports visit
## Ordered by port visit number et the end
function SelectNewVisitOrderedShips(inst::Instance, sol::Sol, paramchosen::ChosenParameters,allparam::AllParameters, paramfixed::FixedParameters, list_visits::Vector{Tuple}) #, placed::Vector{Tuple{Int64,Int64}})
    #TODO: t1 cannot be greater than t2!
    # adapt to given planned previous visits for each ship
    @unpack N, h, qli, shipsIn, Pi, Bp, delta, dist, T, maxC, maxT = inst
    @unpack M = sol
    @unpack Alpha, Proba, AverageCost, Nbexp, Q = allparam
    # C = maximum(length.(Pi))
    new_visits = Vector{NewVisit}()
    new_visits_constrained = Vector{NewVisit}()
    max_cost=0
    max_time = 0
    max_dist=0
    min_cost=1000000000
    min_time=1000000000
    min_dist=1000000000
    max_cost_const=0
    max_time_const = 0
    max_dist_const=0
    min_cost_const=1000000000
    min_time_const=1000000000
    min_dist_const=1000000000
    for (n,c) in list_visits
        sch = sol.visits[n]
        if sch[c].planned == false
            visit, visit_constrained, feasible=SelectNewVisitPerShip(inst, sol, paramchosen, allparam, paramfixed, n, c, Pi[n][c], false)
            if feasible
                if visit.n!=0
                    this_cost = visit.cost
                    this_time = visit.time
                    this_dist = visit.distance
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
                    if this_dist>max_dist
                        max_dist=this_dist
                    end
                    if this_dist<min_dist
                        min_dist=this_dist
                    end
                    push!(new_visits,visit)
                end
                if visit_constrained.n!=0
                    this_cost = visit_constrained.cost
                    this_time = visit_constrained.time
                    this_dist = visit_constrained.distance
                    if this_cost>=max_cost_const
                        max_cost_const=this_cost
                    end
                    if this_cost<=min_cost_const
                        min_cost_const=this_cost
                    end
                    if this_time>max_time_const
                        max_time_const=this_time
                    end
                    if this_time<min_time_const
                        min_time_const=this_time
                    end
                    if this_dist>max_dist_const
                        max_dist_const=this_dist
                    end
                    if this_dist<min_dist_const
                        min_dist_const=this_dist
                    end
                    push!(new_visits_constrained,visit_constrained)
                end
            end
        end
    end
    pos_chosen=[]
    if paramchosen.Reversed=="no"
        if length(new_visits)>0
            tactic=paramchosen.TacticAllBoats
            if tactic=="cost"
                alpha_value = Alpha.CostAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits
                    if visit.cost<=min_cost+alpha_value*(max_cost-min_cost)
                        visit.store.tacticAll="cost"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen = list_available_shippos
            end

            if tactic=="dist"
                alpha_value = Alpha.DistAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits
                    if visit.distance<=min_dist+alpha_value*(max_dist-min_dist)
                        visit.store.tacticAll="dist"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen = list_available_shippos
            end

            if tactic=="time"
                alpha_value = Alpha.TimeAllShip[paramchosen.IndexAllShip]
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
        if length(new_visits_constrained)>0
            tactic=paramchosen.TacticAllBoats
            if tactic=="cost"
                alpha_value = Alpha.CostAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits_constrained
                    if visit.cost<=min_cost_const+alpha_value*(max_cost_const-min_cost_const)
                        visit.store.tacticAll="cost"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen_constrained = list_available_shippos
            end

            if tactic=="dist"
                alpha_value = Alpha.DistAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits_constrained
                    if visit.distance<=min_dist_const+alpha_value*(max_dist_const-min_dist_const)
                        visit.store.tacticAll="dist"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen_constrained = list_available_shippos
            end

            if tactic=="time"
                alpha_value = Alpha.TimeAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits_constrained
                    if visit.time<=min_time_const+alpha_value*(max_time_const-min_time_const)
                        visit.store.tacticAll="time"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen_constrained = list_available_shippos
            end
        end
    else
        if length(new_visits)>0
            tactic=paramchosen.ReversedTacticAllBoats
            if tactic=="cost"
                alpha_value = Alpha.ReversedCostAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits
                    if visit.cost>=min_cost+alpha_value*(max_cost-min_cost)
                        visit.store.tacticAll="cost"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen = list_available_shippos
            end

            if tactic=="dist"
                alpha_value = Alpha.ReversedDistAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits
                    if visit.distance>=min_dist+alpha_value*(max_dist-min_dist)
                        visit.store.tacticAll="dist"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen = list_available_shippos
            end

            if tactic=="time"
                alpha_value = Alpha.ReversedTimeAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits
                    if visit.time>=min_time+alpha_value*(max_time-min_time)
                        visit.store.tacticAll="time"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen = list_available_shippos
            end
        end

        pos_chosen_constrained=[]
        if length(new_visits_constrained)>0
            tactic=paramchosen.ReversedTacticAllBoats
            if tactic=="cost"
                alpha_value = Alpha.ReversedCostAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits_constrained
                    if visit.cost>=min_cost_const+alpha_value*(max_cost_const-min_cost_const)
                        visit.store.tacticAll="cost"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen_constrained = list_available_shippos
            end

            if tactic=="dist"
                alpha_value = Alpha.ReversedDistAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits_constrained
                    if visit.distance>=min_dist_const+alpha_value*(max_dist_const-min_dist_const)
                        visit.store.tacticAll="dist"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen_constrained = list_available_shippos
            end

            if tactic=="time"
                alpha_value = Alpha.ReversedTimeAllShip[paramchosen.IndexAllShip]
                list_available_shippos = Vector{NewVisit}()
                for visit in new_visits_constrained
                    if visit.time>=min_time_const+alpha_value*(max_time_const-min_time_const)
                        visit.store.tacticAll="time"
                        push!(list_available_shippos, visit)                  
                    end
                end
                pos_chosen_constrained = list_available_shippos
            end
        end
    end

    if length(pos_chosen)+length(pos_chosen_constrained)>0
        if length(pos_chosen_constrained)>=1
            random_value = rand()
            if random_value<Alpha.RateConstrained[paramchosen.IndexRateConstrained]
                return pos_chosen_constrained[rand(1:length(pos_chosen_constrained))], true
            else
                return pos_chosen[rand(1:length(pos_chosen))], true
            end
        else
            return pos_chosen[rand(1:length(pos_chosen))], true
        end
    else
        return NewVisit(0,0,0,0,0.0,0.0,0.0, false, ToStoreVisit(SplitCosts(0,0,0,0,0,0),0,"","")), false
    end
end


################### Types of constructions
function greedyrandomizedconstruction(inst::Instance, paramchosen::ChosenParameters, allparam::AllParameters, paramfixed::FixedParameters, max_time)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT = inst
    sol = initializeSol(inst, allparam)
    continue_=true
    start = time_ns()
    elapsed = round((time_ns()-start)/1e9,digits=3)
    count_when=0
    list_to_be_visited = Vector{Tuple}()
    for n in 1:N
        for (c,p) in enumerate(Pi[n])
            push!(list_to_be_visited, (n,c))
        end
    end
    list_to_be_visited_first = Vector{Tuple}()
    for n in 1:N
        push!(list_to_be_visited_first, (n,1))
    end
    new_visit, feasible = SelectNewVisitOrderedShips(inst, sol, paramchosen, allparam, paramfixed, list_to_be_visited_first)
    while continue_ && elapsed<max_time
        @unpack n,c,b,t,cost,distance,constrained,store= new_visit
        deleteat!(list_to_be_visited, findall(x->x==(n,c),list_to_be_visited))
        deleteat!(list_to_be_visited_first, findall(x->x==(n,c),list_to_be_visited_first))
        if c<length(Pi[n])
            push!(list_to_be_visited_first, (n,c+1))
        end
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
            new_visit, feasible = SelectNewVisitOrderedShips(inst, sol, paramchosen, allparam, paramfixed, list_to_be_visited_first)
            if feasible==false
                feasible=true
                new_visit, feasible = SelectNewVisitAllShips(inst, sol, paramchosen, allparam, paramfixed)
            end
            elapsed = round((time_ns()-start)/1e9,digits=3)
            if feasible==false
                sol.failed=1
                nb_to_remove_visits=0
                nb_to_remove_ships=0
                #print('\n')
                #print("###################")
                #print('\n')
                #print("We did not find a at first try solution")
                return initializeSol(inst, allparam)
            end
        else
            elapsed = round((time_ns()-start)/1e9,digits=3)
            #print('\n')
            #print("The time it took : $elapsed")
        end
    end
    return sol
end

function makeSwapsFromOrder(inst::Instance,sol::Sol, when_list, nb_swap)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT = inst
    @unpack M, visits = sol
    for i in 1:nb_swap
        index_first_el_to_swap = rand(1:length(when_list))
        first_el_to_swap = deepcopy(when_list[index_first_el_to_swap])
        if first_el_to_swap[2]>1
            previous_visit_index = findfirst( x -> x == (first_el_to_swap[1],first_el_to_swap[2]-1), when_list )
            if previous_visit_index<length(when_list)
                index_second_el_to_swap = rand(previous_visit_index+1:length(when_list))
                second_el_to_swap = when_list[index_second_el_to_swap]
                if second_el_to_swap[2]>1
                    second_previous_visit_index = findfirst( x -> x == (second_el_to_swap[1],second_el_to_swap[2]-1), when_list )
                    if second_previous_visit_index<index_first_el_to_swap
                        when_list[index_first_el_to_swap] = deepcopy(second_el_to_swap)
                        when_list[index_second_el_to_swap] = deepcopy(first_el_to_swap)
                    end
                end
            end
        else
            index_second_el_to_swap = rand(1:length(when_list))
            second_el_to_swap = when_list[index_second_el_to_swap]
            if second_el_to_swap[2]>1
                second_previous_visit_index = findfirst( x -> x == (second_el_to_swap[1],second_el_to_swap[2]-1), when_list )
                print(second_previous_visit_index)
                if second_previous_visit_index<index_first_el_to_swap
                    when_list[index_first_el_to_swap] = deepcopy(second_el_to_swap)
                    when_list[index_second_el_to_swap] = deepcopy(first_el_to_swap)
                end
            end
        end
    end
    return when_list
end


function greedyorderedconstruction(inst::Instance, paramchosen::ChosenParameters, allparam::AllParameters, paramfixed::FixedParameters, when_list, max_time)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT = inst
    sol = initializeSol(inst, allparam)
    continue_=true
    start = time_ns()
    elapsed = round((time_ns()-start)/1e9,digits=3)
    count_when=0
    n=when_list[1][1]
    c=when_list[1][2]   
    p=Pi[n][c]
    new_visit, new_visit_constrained, feasible=SelectNewVisitPerShip(inst, sol, paramchosen, allparam, paramfixed, n, c, p, false)
    when_list=when_list[min(2,length(when_list)):length(when_list)]
    start_heur = time_ns()
    while continue_ && elapsed<max_time
        @unpack n,c,b,t,cost,distance,constrained,store= new_visit
        count_when=count_when+1
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
            n=when_list[1][1]
            c=when_list[1][2]   
            p=Pi[n][c]
            if sol.visits[n][c].planned == false
                new_visit, new_visit_constrained, feasible=SelectNewVisitPerShip(inst, sol, paramchosen, allparam, paramfixed, n, c, p, false)
                if feasible
                    if length(when_list)>1
                        when_list=when_list[2:length(when_list)]
                    else
                        when_list=[]
                    end
                end
                if feasible==false
                    feasible=true
                    new_visit, feasible = SelectNewVisitAllShips(inst, sol, paramchosen, allparam, paramfixed)
                end
            else
                if length(when_list)>1
                    when_list=when_list[2:length(when_list)]
                else
                    when_list=[]
                end
            end
            elapsed = round((time_ns()-start)/1e9,digits=3)
            if feasible==false
                sol.failed=1
                nb_to_remove_visits=0
                nb_to_remove_ships=0
                #print('\n')
                #print("###################")
                #print('\n')
                #print("We did not find a at first try solution")
                return initializeSol(inst, allparam), false
            end
        else
            #print('\n')
            #print("The time it took : $elapsed")
        end
    end
    #print('\n')
    #print("Time heur")
    #print('\n')
    #elapsed = round((time_ns()-start_heur)/1e9,digits=3)
    #print(elapsed_heur)
    return sol, true
end

function greedyremoverandomconstruction(inst::Instance, this_sol::Sol,  paramchosen::ChosenParameters, allparam::AllParameters, paramfixed::FixedParameters, frompath, max_time)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT = inst
    @unpack Alpha, Proba, AverageCost, Nbexp, Q = allparam
    sol=deepcopy(this_sol)
    continue_=true
    start = time_ns()
    elapsed = round((time_ns()-start)/1e9,digits=3)

    boats_to_remove = Vector{}()
    for n in 1:N
        push!(boats_to_remove, n)
    end
    if frompath
        boats_to_remove = sample(boats_to_remove, ceil(Int,paramfixed.RemovePathRelinking*N); replace=false)
    else
        boats_to_remove = sample(boats_to_remove, ceil(Int, Alpha.PropToRemove[paramchosen.IndexPropToRemove]*N); replace=false)
    end
    for n in boats_to_remove
        for (c,p) in enumerate(Pi[n])
            sol.visits[n][c].p = -1
            sol.visits[n][c].b = -1
            sol.visits[n][c].t = -1
            sol.visits[n][c].planned = false
            sol.visits[n][c].minT = max(shipsIn[n].sT[c], T[n,c,1])
            sol.visits[n][c].maxT = maxT
        end
    end

    count_when=0
    list_to_be_visited = Vector{Tuple}()
    for n in 1:N
        for (c,p) in enumerate(Pi[n])
            if sol.visits[n][c].planned == false
                push!(list_to_be_visited, (n,c))
            end
        end
    end
    list_to_be_visited_first = Vector{Tuple}()
    for n in 1:N
        this_boat_first = true
        for (c,p) in enumerate(Pi[n])
            if sol.visits[n][c].planned == false && this_boat_first
                push!(list_to_be_visited_first, (n,c))
                this_boat_first = false
            end
        end
    end
    sol.M=generateOccupiedMx(inst, sol.visits)
    new_visit, feasible = SelectNewVisitOrderedShips(inst, sol, paramchosen, allparam, paramfixed, list_to_be_visited_first)
    while continue_ && elapsed<max_time
        @unpack n,c,b,t,cost,distance,constrained,store= new_visit
        deleteat!(list_to_be_visited, findall(x->x==(n,c),list_to_be_visited))
        deleteat!(list_to_be_visited_first, findall(x->x==(n,c),list_to_be_visited_first))
        if c<length(Pi[n])
            push!(list_to_be_visited_first, (n,c+1))
        end
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
            new_visit, feasible = SelectNewVisitOrderedShips(inst, sol, paramchosen, allparam, paramfixed, list_to_be_visited_first)
            if feasible==false
                feasible=true
                new_visit, feasible = SelectNewVisitAllShips(inst, sol, paramchosen, allparam, paramfixed)
            end
            elapsed = round((time_ns()-start)/1e9,digits=3)
            if feasible==false
                sol.failed=1
                nb_to_remove_visits=0
                nb_to_remove_ships=0
                #print('\n')
                #print("###################")
                #print('\n')
                #print("We did not find a at first try solution")
                return initializeSol(inst, allparam)
            end
        else
            elapsed = round((time_ns()-start)/1e9,digits=3)
            #print('\n')
            #print("The time it took : $elapsed")
        end
    end
    return sol
end



#### Complete GRASP
function GRASP_reactive(seed,N,Nout,qli, type1, type2, type3, paramfixed, temperature, time_local, max_time_heur, max_time, expname, location)
    inst = readInstFromFile(location*"MCBAP-multi-port-berth-allocation-problem/Large/CP2_Inst_$seed"*"_$N"*"_$Nout"*"_$qli"*".txt")
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp = inst
    cost=1000000000
    worst_cost=10000000
    allparam = initializeParam(paramfixed)
    sol=initializeSol(inst, allparam)
    best_sol=initializeSol(inst, allparam)
    start = time_ns()
    elapsed = round((time_ns()-start)/1e9,digits=3)
    start_iter=time_ns()
    elapsed_iter = round((time_ns()-start_iter)/1e9,digits=3)
    nb_iter=0
    d_alliter_before=Dict()
    d_alliter_after=Dict()
    when_list = Vector{Tuple}()
    when_dict=Dict()
    count_better_heur_solution = 0
    min_cost_heur=1000000000
    best_cost=1000000000
    best_reconstruct_cost=1000000000
    best_relinking_cost=1000000000
    nb_iter_reconstruct=0
    reconstruct_no_improve=0
    greedy_no_improve=0
    relinking_no_improve=0
    usedcplex=0
    from_reconstruct=false
    from_reconstruct_iter=0
    nb_iter_restart_params = 0
    first_relinking = true
    elite_pool = Vector{}()
    elite_pool_heur = Vector{}()
    distance_btw_sols_min=0
    for n in 1:N
        when_dict[n]=Dict()
        for c in 1:length(inst.Pi[n])
            when_dict[n][c]=0
            
        end
    end
    for n in 1:N
        for c in 1:length(inst.Pi[n])
            push!(when_list,(n,c,1))            
        end
    end
    proba_temperature=1
    first_good_solution= false

    not_first_time=false
    focus_on_remove=false
    min_distance_elite=10000000

    while elapsed<max_time
        nb_iter_restart_params = nb_iter_restart_params + 1
        paramchosen = ChooseAllParam(allparam, paramfixed)
        start_heur = time_ns()
        random_value = rand()

        if ((random_value <= proba_temperature || nb_iter_reconstruct < paramfixed.Until) && (greedy_no_improve<paramfixed.GreedymaxNoImprove)) #|| first_good_solution == false
            proba_temperature = proba_temperature*temperature
            print('\n')
            print("First step")
            print('\n')
            print("GREEDY")
            print('\n')
            print("The temperature is : $proba_temperature")
            print('\n')
            print("=======================================")
            new_sol = greedyrandomizedconstruction(inst, paramchosen, allparam, paramfixed, max_time_heur)
            feasible=true
            for n in 1:N
                # The times :
                for (c,p) in enumerate(inst.Pi[n])
                    if new_sol.visits[n][c].planned == false
                        feasible = false              
                    end
                end
            end
            if feasible
                new_cost, delay_cost_heur, waiting_cost_heur, penalty_cost_heur, handling_cost_heur, fuel_cost_heur = checkSolutionCost(inst, new_sol)
            end
        else    
            print('\n')
            print("RECONSTRUCTION")
            print('\n')
            print("The temperature is : $proba_temperature")
            print('\n')
            print("#############################")
            proba_temperature = 1
            if reconstruct_no_improve>paramfixed.FocusRemoveUntil
                focus_on_remove=true
            end

            if focus_on_remove==false
                new_sol = greedyremoverandomconstruction(inst, sol, paramchosen, allparam, paramfixed, false, max_time_heur)
                feasible=true
                for n in 1:N
                    # The times :
                    for (c,p) in enumerate(inst.Pi[n])
                        if new_sol.visits[n][c].planned == false
                            feasible = false              
                        end
                    end
                end
                if feasible
                    new_cost, delay_cost_heur, waiting_cost_heur, penalty_cost_heur, handling_cost_heur, fuel_cost_heur = checkSolutionCost(inst, new_sol)
                    allparam = UpdateAfterReconstructParameters(paramchosen, allparam, cost, new_cost, paramfixed.lookforconstrained)
                end
                
            else
                foundbetter=false
                for i in 1:paramfixed.NbFocusRemove
                    new_sol = greedyremoverandomconstruction(inst, sol, paramchosen, allparam, paramfixed, false, max_time_heur)
                    feasible=true
                    for n in 1:N
                        # The times :
                        for (c,p) in enumerate(inst.Pi[n])
                            if new_sol.visits[n][c].planned == false
                                feasible = false              
                            end
                        end
                    end
                    if feasible
                        new_cost, delay_cost_heur, waiting_cost_heur, penalty_cost_heur, handling_cost_heur, fuel_cost_heur = checkSolutionCost(inst, new_sol)
                        allparam = UpdateAfterReconstructParameters(paramchosen, allparam, cost, new_cost, paramfixed.lookforconstrained) 
                    end
                    if new_cost<sol.store.costLocal.all
                        foundbetter=true
                        sol = deepcopy(new_sol)
                        cost=deepcopy(new_cost)
                    end
                end
                if foundbetter
                    new_sol = deepcopy(sol)
                    new_cost=deepcopy(cost)
                end
            end
            #new_sol=deepcopy(sol)
            from_reconstruct=true
            from_reconstruct_iter=1
        end
        proba_temperature = proba_temperature*temperature
        elapsed_heur = round((time_ns()-start_heur)/1e9,digits=3)
        feasible = true
        print('\n')
        print("Time heur")
        print('\n')
        print(elapsed_heur)
        ## No conflicts with within one boat schedule :
        for n in 1:N
            # The times :
            for (c,p) in enumerate(inst.Pi[n])
                if new_sol.visits[n][c].planned == false
                    feasible = false              
                end
            end
        end
        #print('\n')
        #print("Second step")
        #print('\n')
        #print("###########################")
        if feasible && checkSolutionFeasability(inst, new_sol)
            new_cost_heur, delay_cost_heur, waiting_cost_heur, penalty_cost_heur, handling_cost_heur, fuel_cost_heur = checkSolutionCost(inst, new_sol)
            
            #print('\n')
            #print("Before local search")
            #print('\n')
            #print(new_sol.visits)
            #d_before=prepareSolIter(seed,N,Nout,qli,nb_iter,inst, new_sol, new_cost_heur, allparam, paramfixed, expname)
            #d_alliter_before[nb_iter]=d_before
            start_local=time_ns()

            print('\n')
            print("New cost heur")
            print('\n')
            print(new_cost_heur)
            print('\n')
            print("Min cost heur")
            print('\n')
            print(min_cost_heur)
            #print('\n')
            #print("Does it come from a reconstruct?")
            #print('\n')
            #print(from_reconstruct)
            print('\n')
            print("Does it come from pathrelinking?")
            print('\n')
            print(new_sol.relinking)
            
            if from_reconstruct==false
                #print('\n')
                #print("The new cost without reconstruct ;")
                #print('\n')
                #print(new_cost_heur)
                if new_cost_heur<min_cost_heur
                    min_cost_heur=deepcopy(new_cost_heur)
                    best_sol_heur=deepcopy(new_sol)
                end
                if new_cost_heur>=min_cost_heur && nb_iter_reconstruct>paramfixed.Until
                    greedy_no_improve=greedy_no_improve+1
                else
                    greedy_no_improve=0
                end
            end

            if new_cost_heur<min_cost_heur+paramfixed.windowLocalSearch*min_cost_heur
                #new_sol, new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost = pushTime(inst, new_sol, paramfixed, time_local)
                new_sol, new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost = manualLocalSearch(inst, new_sol, new_cost_heur, delay_cost_heur, waiting_cost_heur, penalty_cost_heur, handling_cost_heur, fuel_cost_heur, allparam, paramfixed, paramchosen, time_local, reconstruct_no_improve, cost)
                new_sol.usedLocalSearch=1
                #new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost = new_cost_heur, delay_cost_heur, waiting_cost_heur, penalty_cost_heur, handling_cost_heur, fuel_cost_heur
            else
                new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost = new_cost_heur, delay_cost_heur, waiting_cost_heur, penalty_cost_heur, handling_cost_heur, fuel_cost_heur
            end
            
            if paramfixed.pathRelinking=="yes" && new_sol.usedLocalSearch==1
                print('\n')
                print("PATH RELINKING")
                print('\n')
                print("The temperature is : $proba_temperature")
                print('\n')
                print("#############################")
                new_sol.relinking==1
                if length(elite_pool)>0
                    sol_elite = deepcopy(rand(elite_pool)[1])
                else
                    sol_elite = deepcopy(sol)
                end
                distance_btw_sols = DistanceSols(inst, new_sol, sol_elite)
                if first_relinking
                    distance_btw_sols_min=deepcopy(distance_btw_sols)
                    first_relinking=false
                end
                start_time_relinking = time_ns()
                elapsed_relinking = round((time_ns()-start_time_relinking)/1e9,digits=3)
                new_sol_path=deepcopy(new_sol)
                delta_distance = deepcopy(distance_btw_sols)
                while elapsed_relinking<paramfixed.MaxTimeRelinking && distance_btw_sols>distance_btw_sols_min
                    new_sol_path_loop = greedyremoverandomconstruction(inst, new_sol_path, paramchosen, allparam, paramfixed, true, max_time_heur)
                    feasible=true
                    for n in 1:N
                        # The times :
                        for (c,p) in enumerate(inst.Pi[n])
                            if new_sol.visits[n][c].planned == false
                                feasible = false              
                            end
                        end
                    end
                    if feasible
                        new_cost_path_loop, delay_cost_path_loop, waiting_cost_path_loop, penalty_cost_path_loop, handling_cost_path_loop, fuel_cost_path_loop = checkSolutionCost(inst, new_sol_path_loop)
                        distance_btw_sols = DistanceSols(inst, new_sol_path_loop, sol_elite)
                        print("||||||||||||||||||||||||||||")
                        print('\n')
                        print("Distance btw sol path and sol elite")
                        print('\n')
                        print(distance_btw_sols)
                        print('\n')
                        print("New cost path loop")
                        print('\n')
                        print(new_cost_path_loop)
                        if distance_btw_sols<delta_distance
                            delta_distance=deepcopy(distance_btw_sols)
                            new_sol_path=deepcopy(new_sol_path_loop)
                        end
                        if new_cost_path_loop<new_cost
                            new_sol=deepcopy(new_sol_path)
                            new_cost=new_cost_path_loop
                            delay_cost=delay_cost_path_loop
                            waiting_cost=waiting_cost_path_loop
                            penalty_cost=penalty_cost_path_loop
                            handling_cost=handling_cost_path_loop
                            fuel_cost=fuel_cost_path_loop
                        end
                    end
                    elapsed_relinking = round((time_ns()-start_time_relinking)/1e9,digits=3)
                end
                if distance_btw_sols<distance_btw_sols_min
                    distance_btw_sols_min=deepcopy(distance_btw_sols)
                end           
            end
            if from_reconstruct==true
                if new_cost<best_reconstruct_cost-paramfixed.RateImproveReconstructOrPath*best_reconstruct_cost
                    best_reconstruct_cost=new_cost
                    reconstruct_no_improve=0
                else
                    reconstruct_no_improve+=1
                end
            end

            #if new_sol.relinking==1
            #    if new_cost<best_relinking_cost-paramfixed.RateImproveReconstructOrPath*best_relinking_cost
            #        best_relinking_cost=deepcopy(new_cost)
            #        relinking_no_improve=0
            #    else
            #        relinking_no_improve+=1
            #    end
            #end
            
            from_reconstruct=false
            
            if nb_iter_reconstruct<paramfixed.Until
                if new_cost<min_cost_heur
                    min_cost_heur=new_cost
                end
            end
            #new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost = new_cost_heur, delay_cost_heur, waiting_cost_heur, penalty_cost_heur, handling_cost_heur, fuel_cost_heur

            elapsed_local = round((time_ns()-start_local)/1e9,digits=3)
            #print('\n')
            #print("Time local search")
            #print('\n')
            #print(elapsed_local)
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
            for n in 1:N
                for (c,p) in enumerate(Pi[n])
                    t = new_sol.visits[n][c].t
                    b = new_sol.visits[n][c].b
                    this_cost, delay_cost, waiting_cost, penalty, handling_cost, fuel_cost, feas = computeCostPosSol(inst, n, c, b, t,new_sol)
                    new_sol.visits[n][c].store.cost =  SplitCosts(ceil(Int, this_cost), ceil(Int,delay_cost), ceil(Int,waiting_cost), ceil(Int,penalty), ceil(Int,handling_cost), ceil(Int,fuel_cost))
                end
            end
            #print('\n')
            #print("The tactics")
            #print('\n')
            #print(usedcplex)
            #print('\n')
            #print(from_reconstruct_iter)
            
            if from_reconstruct_iter==1
                new_sol.reconstruct=1
            end
            print('\n')
            print("(((((((((((((((((())))))))))))))))))")
            print('\n')
            print("Proba tactic one boat")
            print(allparam.Proba.TacticOneBoat)
            print('\n')
            print("Proba tactic all boats")
            print('\n')
            print(allparam.Proba.TacticAllBoats)
            print("Proba tactic reversed")
            print('\n')
            print(allparam.Proba.Reversed)
            print("(((((((((((((((((())))))))))))))))))")
            if new_sol.reconstruct==1
                allparam = UpdateAfterReconstructParameters(paramchosen, allparam, cost, new_cost, paramfixed.lookforconstrained)   
            else
                allparam = UpdateAfterHeurParameters(paramchosen, allparam, cost, new_cost, paramfixed.lookforconstrained)   
            end

            if new_cost<cost
                #print('\n')
                #print("################# lol")
                sol.better=1
                cost=deepcopy(new_cost)
                sol=deepcopy(new_sol)
            end
            if new_sol.reconstruct!=1
                if length(elite_pool)<paramfixed.LengthElite
                    distance_sols_min = 10000000
                    for sol_elite in elite_pool
                        distance_sols = DistanceSols(inst, sol_elite[1], new_sol)
                        if distance_sols<distance_sols_min
                            distance_sols_min=distance_sols
                        end
                    end
                    if distance_sols_min<min_distance_elite
                        min_distance_elite=distance_sols_min
                    end
                    push!(elite_pool, (deepcopy(new_sol),new_cost,distance_sols_min))                   
                else
                    elite_pool = sort(elite_pool, by=x->x[2])
                    if new_cost<elite_pool[end][2]
                        distance_sols_min = 10000000
                        worst_distance = elite_pool[end][3]
                        for sol_elite in elite_pool
                            distance_sols = DistanceSols(inst, sol_elite[1], new_sol)
                            if distance_sols<distance_sols_min
                                distance_sols_min=distance_sols
                            end
                        end
                        if distance_sols_min>min_distance_elite
                            print("__________________")
                            print("new el in pool")
                            elite_pool[end]=(deepcopy(new_sol),new_cost, distance_sols_min)
                            elite_pool = sort(elite_pool, by=x->x[3])
                            min_distance_elite=deepcopy(elite_pool[1][3])
                        end  
                    end
                end
                #if length(elite_pool)<paramfixed.LengthElite
                #    push!(elite_pool, (deepcopy(new_sol),new_cost))
                #else
                #    elite_pool = sort(elite_pool, by=x->x[2])
                #    if new_cost<elite_pool[end][2]
                #        elite_pool[end]=(deepcopy(new_sol),new_cost)
                #        elite_pool = sort(elite_pool, by=x->x[2])
                #    end
                #end
            end
            list_costs_elite= Vector{}()
            list_dist_elite= Vector{}()
            for el in elite_pool
                push!(list_costs_elite, el[2])
                push!(list_dist_elite, el[3])
            end
            sol.average_cost_elite=mean(list_costs_elite)
            sol.average_dist_elite=mean(list_dist_elite)
            d_after = prepareSolIterSoft(seed,N,Nout,qli,nb_iter,inst, new_sol, new_cost, allparam, paramchosen, expname)
            
            print('\n')
            print("###################################")
            print(greedy_no_improve)
            #print('\n')
            #print("#######################")
            #print(reconstruct_no_improve)
            print('\n')
            print("#######################")
            print(relinking_no_improve)
            
            if reconstruct_no_improve>paramfixed.maxNoImprove
                if cost<best_cost
                    best_cost=deepcopy(cost)
                    best_sol=deepcopy(sol)
                end
                sol = greedyrandomizedconstruction(inst, paramchosen, allparam, paramfixed, max_time_heur)
                min_cost_heur=1000000000
                best_reconstruct_cost=1000000000
                cost=1000000000
                nb_iter_reconstruct=0
                reconstruct_no_improve=0
                proba_temperature = 1
                greedy_no_improve=0
                focus_on_remove=false
            end
            from_reconstruct_iter=0
            usedcplex=0
            d_alliter_after[nb_iter]=d_after
            nb_iter+=1
            nb_iter_reconstruct+=1

            if nb_iter_restart_params>40
                allparam= RestartParamNb(allparam)
            end

                        
                         
        else
            print('\n')
            print("We had unfeasability")
        end
        elapsed = round((time_ns()-start)/1e9,digits=3)
    #print('\n')
    #print("End")
    #print('\n')
    #print("###########################")
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

    if best_cost<cost
        cost=deepcopy(best_cost)
        sol=deepcopy(best_sol)
    end

    if feasible && checkSolutionFeasability(inst, sol)
        #print('\n')
        #print("The solution is feasible")
        for iter in keys(d_alliter_after)
            #CSV.write(location*"results_jobs/benchmarks_HEUR/final/$expname"*"/iterations_before_local/sol_$seed"*"_$N"*"_$Nout"*"_$qli"*"/iter_$iter"*"_beforelocal"*".csv", d_alliter_before[iter])
            CSV.write(location*"results_jobs/benchmarks_HEUR/final/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli"*"/iter_$iter"*".csv", d_alliter_after[iter])
        end
        #print('\n')
        #print("The number of iterations is :")
        #print('\n')
        #print(nb_iter)
        nb_iter+=1
        
        d = prepareSolIter(seed,N,Nout,qli,nb_iter,inst, sol, cost, allparam, paramfixed, expname)
        CSV.write(location*"results_jobs/benchmarks_HEUR/final/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli"*"/iter_$nb_iter"*".csv", d)
        return sol, cost, allparam
    else
        cost=1000000000
        allparam = initializeParam(paramfixed)
        sol=initializeSol(inst, allparam)
        return sol, cost, allparam
    end
end




######## Test function :
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
    prop_oneboattime = 0.2
    prop_allboatcost = 0.2
    prop_allboatdist = 0.2
    prop_allboattime = 0.5

    # The window size for the visits to look at :
    window=0.2

    # The prop of ships to remove for the reconstruction :
    proptoremove=0.25

    # Should we push at constraints position :
    pushatconstraint=true

    # Look for constrained at the beginning :
    lookforconstraint=true

    # The proportion of boats to remove for the local search :
    alphaboatmin=0.1
    alphaboatmax=0.4
    alpharandommin=0.1
    alpharandommax=0.4

    # All the parameters :
    paramfixed=FixedParameters(alpharandommin,alpharandommax,alphaboatmin,alphaboatmax,proptoremove,window, pushatconstraint, lookforconstraint)

    # Maximum time for the local search :
    time_local=5

    # Maximum time for the heuristic :
    max_time_heur=30

    # maximum time for the experiment :
    max_time=300

    # the temperature parameter :
    temperature=0.93


    #print("Start")
    #print('\n')
    for instance_name in all_instances
        split_instance = split(instance_name,"_")
        seed=parse(Int64,split_instance[3])
        N=parse(Int64,split_instance[4])
        Nout=parse(Int64,split_instance[5])
        qli=parse(Int64,split(split_instance[6],".")[1])
        print("The instance : $seed"*"_$N"*"_$Nout"*"_$qli")
        inst = readInstFromFile(location*"MCBAP-multi-port-berth-allocation-problem/Large/CP2_Inst_$seed"*"_$N"*"_$Nout"*"_$qli"*".txt")
        sol, cost, allparam = GRASP_reactive(seed,N,Nout,qli, type1, type2, type3,  paramfixed, temperature, time_local, max_time_heur, max_time, expname, location)
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
            #print('\n')
            #print("The probas :")
            #print('\n')
            #print(allparam.Proba)
            #print('\n')
            #print("The solution :")
            #print('\n')
            #print(sol.visits)
            #print('\n')
            print("And the cost is ")
            print('\n')
            print(cost)
        end
    end
end
