using UnPack 
using Statistics
import XLSX
using CSV, Tables
using DataFrames
using Random
using StatsBase
include("../../MBAP_INST.jl")
include("../get_iterations.jl")
include("../MBAP_SOL.jl")
include("../toolsMatrixTimes.jl")
include("../check_solution.jl")
include("../utilInit.jl")
include("localsearch.jl")
include("constrainedPos.jl")

################################ GRAPS algorithm

#################### The functions for the two main steps

## The first step 
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
    random_value = rand()
    if paramfixed.lookforconstrained==false || random_value>=Alpha.RateConstrained[paramchosen.IndexRateConstrained]
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
    end  
    if paramfixed.lookforconstrained==true && random_value<Alpha.RateConstrained[paramchosen.IndexRateConstrained]
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


## The second step
function SelectNewVisitAllShips(inst::Instance, sol::Sol, paramchosen::ChosenParameters,allparam::AllParameters, paramfixed::FixedParameters, list_visits::Vector{Tuple}) #, placed::Vector{Tuple{Int64,Int64}})
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
            return reduce(vcat, (pos_chosen,pos_chosen_constrained))[rand(1:length(pos_chosen)+length(pos_chosen_constrained))], true
        else
            return pos_chosen[rand(1:length(pos_chosen))], true
        end
    else
        return NewVisit(0,0,0,0,0.0,0.0,0.0, false, ToStoreVisit(SplitCosts(0,0,0,0,0,0),0,"","")), false
    end
end


#################### The functions which iterates these two steps
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
    new_visit, feasible = SelectNewVisitAllShips(inst, sol, paramchosen, allparam, paramfixed, list_to_be_visited_first)
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
            new_visit, feasible = SelectNewVisitAllShips(inst, sol, paramchosen, allparam, paramfixed, list_to_be_visited_first)
            elapsed = round((time_ns()-start)/1e9,digits=3)
            if feasible==false
                sol.failed=1
                nb_to_remove_visits=0
                nb_to_remove_ships=0
                return initializeSol(inst, allparam)
            end
        else
            elapsed = round((time_ns()-start)/1e9,digits=3)
        end
    end
    return sol
end