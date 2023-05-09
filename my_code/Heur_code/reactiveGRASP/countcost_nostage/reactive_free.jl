using UnPack 
using Statistics
import XLSX
using CSV, Tables
using DataFrames
using Random
using StatsBase
include("../../../MBAP_INST.jl")
include("../../../get_iterations.jl")
include("../../MBAP_SOL.jl")
include("../../toolsMatrixTimes.jl")
include("../../check_solution.jl")
include("utilInit.jl")
include("../../localSearch/localCPLEX.jl")

mutable struct NewVisit
    n::Int64    # ship number
    c::Int64    # visit number
    b::Int64  #postion
    t::Int64    # port
    cost::Float64    
    distance::Float64   
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
function SelectNewVisitPerShip(inst::Instance, sol::Sol, paramchosen::ChosenParameters, allparam::AllParameters, n, c, p, printdetails) 
    @unpack N, h, qli, shipsIn, Pi, Bp, delta, dist, T, maxC, maxT = inst
    @unpack M = sol
    @unpack Alpha, Proba, AverageCost, Nbexp, Q = allparam
    l = ceil(Int, shipsIn[n].l/qli)
    new_visits = Vector{NewVisit}()
    new_visits_constrained = Vector{NewVisit}()
    max_dist=0
    max_cost=0
    min_dist=1000000000
    min_cost=1000000000
    max_dist_const=0
    max_cost_const=0
    min_dist_const=1000000000
    min_cost_const=1000000000
    find_feasible = false
    t1=sol.visits[n][c].minT
    t2=sol.visits[n][c].maxT
    sch = sol.visits[n]
    bestPos=shipsIn[n].Bi[c]/qli
    ## Check all the possible solutions
    for b in 1:Bp[p]-l
        hand = ceil(Int, h[n][c][b])
        for t in t1:t2
            if t + hand <= t2
                if M[n][c][b,t+1]
                    this_cost, delay_cost, waiting_cost, penalty, handling_cost, fuel_cost, feas = computeCostPosSol(inst, n, c, b, t,sol)
                    if feas
                        distance = abs(bestPos-b)/l
                        if distance>=max_dist
                            max_dist=distance
                        end
                        if this_cost>=max_cost
                            max_cost=this_cost
                        end
                        if distance<=min_dist
                            min_dist=distance
                        end
                        if this_cost<=min_cost
                            min_cost=this_cost
                        end
                        push!(new_visits, NewVisit(n,c,b,t,this_cost,distance, false, ToStoreVisit(SplitCosts(ceil(Int, this_cost), ceil(Int,delay_cost), ceil(Int,waiting_cost), ceil(Int,penalty), ceil(Int,handling_cost), ceil(Int,fuel_cost)),0,"","")))
                    end
                end
            end
        end
    end  
    
    ## Check the constrained ones
    pos = findConstrainedPos(inst, sol, p, t1, t2, n, c, l)
    for (b,t) in pos
        hand = ceil(Int, h[n][c][b])
        if t + hand <= t2 
            if sol.M[n][c][b,t+1]
                this_cost, delay_cost, waiting_cost, penalty, handling_cost, fuel_cost, feas = computeCostPosSol(inst, n, c, b, t ,sol)
                if feas
                    distance = abs(bestPos-b)/l
                    if distance>=max_dist_const
                        max_dist_const=distance
                    end
                    if this_cost>=max_cost_const
                        max_cost_const=this_cost
                    end
                    if distance<=min_dist_const
                        min_dist_const=distance
                    end
                    if this_cost<=min_cost_const
                        min_cost_const=this_cost
                    end
                    push!(new_visits_constrained, NewVisit(n,c,b,t,this_cost,distance, true, ToStoreVisit(SplitCosts(ceil(Int, this_cost), ceil(Int,delay_cost), ceil(Int,waiting_cost), ceil(Int,penalty), ceil(Int,handling_cost), ceil(Int,fuel_cost)),0,"","")))
                end
            end
        end
    end
    if printdetails
        print("lol")
        print('\n')
        print(new_visits)
        print('\n')
        print(new_visits_constrained)
        print('\n')
    end
    pos_chosen=NewVisit(0,0,0,0,0.0,0.0, false, ToStoreVisit(SplitCosts(0,0,0,0,0,0),0,"",""))
    pos_chosen_constrained=NewVisit(0,0,0,0,0.0,0.0, false, ToStoreVisit(SplitCosts(0,0,0,0,0,0),0,"",""))
    tactic=paramchosen.TacticTypeShip
    if length(new_visits)>0
        if tactic=="cost"
            alpha_value= Alpha.CostTypeShip[paramchosen.IndexTypeShip]
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
            alpha_value= Alpha.DistanceTypeShip[paramchosen.IndexTypeShip]
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
    end
    if length(new_visits_constrained)>0
        if tactic=="cost"
            alpha_value= Alpha.CostTypeShip[paramchosen.IndexTypeShip]
            list_available_shippos = Vector{NewVisit}()
            for visit in new_visits_constrained
                if visit.constrained
                    if visit.cost<=min_cost_const+alpha_value*(max_cost_const-min_cost_const)
                        visit.store.tacticBoat="cost"
                        push!(list_available_shippos, visit)
                    end                    
                end
            end
            pos_chosen_constrained = list_available_shippos[rand(1:length(list_available_shippos))]
        end
        if tactic=="dist"
            alpha_value= Alpha.DistanceTypeShip[paramchosen.IndexTypeShip]
            list_available_shippos = Vector{NewVisit}()
            for visit in new_visits_constrained
                if visit.constrained
                    if visit.distance<=min_dist_const+alpha_value*(max_dist_const-min_dist_const)
                        visit.store.tacticBoat="dist"
                        push!(list_available_shippos, visit)
                    end                    
                end
            end
            pos_chosen_constrained = list_available_shippos[rand(1:length(list_available_shippos))]
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
function SelectNewVisitAllShips(inst::Instance, sol::Sol, paramchosen::ChosenParameters, allparam::AllParameters) #, placed::Vector{Tuple{Int64,Int64}})
    #TODO: t1 cannot be greater than t2!
    # adapt to given planned previous visits for each ship
    @unpack N, h, qli, shipsIn, Pi, Bp, delta, dist, T, maxC, maxT = inst
    @unpack M = sol
    @unpack Alpha, Proba, AverageCost, Nbexp, Q = allparam
    # C = maximum(length.(Pi))
    new_visits = Vector{NewVisit}()
    new_visits_constrained = Vector{NewVisit}()
    max_count=0
    max_cost=0
    min_count=1000000000
    min_cost=1000000000
    max_count_const=0
    max_cost_const=0
    min_count_const=1000000000
    min_cost_const=1000000000
    count_available=Dict()
    for n in 1:N
        count_available[n]=Dict()
        for (c,p) in enumerate(Pi[n])
            sch = sol.visits[n]
            if sch[c].planned == false
                visit, visit_constrained, feasible, count_visit=SelectNewVisitPerShip(inst, sol, paramchosen, allparam, n, c, p, false)
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
                        if this_cost>max_cost
                            max_cost=this_cost
                        end
                        if this_cost<min_cost
                            min_cost=this_cost
                        end
                        push!(new_visits,visit)
                    end
                    if visit_constrained.n!=0
                        this_cost = visit_constrained.cost
                        if this_cost>=max_cost_const
                            max_cost_const=this_cost
                        end
                        if this_cost<=min_cost_const
                            min_cost_const=this_cost
                        end
                        push!(new_visits_constrained,visit_constrained)
                    end
                end
            end
        end
    end
    pos_chosen=[]
    if length(new_visits)>0
        tactic=paramchosen.TacticAllShip
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

        if tactic=="count"
            alpha_value = Alpha.CountAllShip[paramchosen.IndexAllShip]
            list_available_shippos = Vector{NewVisit}()
            for visit in new_visits
                if count_available[visit.n][visit.c]<=min_count+alpha_value*(max_count-min_count)
                    visit.store.tacticAll="count"
                    push!(list_available_shippos, visit)                  
                end
            end
            pos_chosen = list_available_shippos
        end
    end
    pos_chosen_constrained=[]
    if length(new_visits_constrained)>0
        tactic=paramchosen.TacticAllShip
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

        if tactic=="count"
            alpha_value = Alpha.CountAllShip[paramchosen.IndexAllShip]
            list_available_shippos = Vector{NewVisit}()
            for visit in new_visits_constrained
                if count_available[visit.n][visit.c]<=min_count_const+alpha_value*(max_count_const-min_count_const)
                    visit.store.tacticAll="count"
                    push!(list_available_shippos, visit)                  
                end
            end
            pos_chosen_constrained = list_available_shippos
        end
    end

    if length(pos_chosen)+length(pos_chosen_constrained)>0
        alpha_value = Alpha.RateConstrained[paramchosen.IndexRateConstrained]
        epsilon=0.0000001
        if length(pos_chosen)>1
            pos_chosen = sample(pos_chosen, ceil(Int,alpha_value*length(pos_chosen)+epsilon))
        end
        if length(pos_chosen_constrained)>1
            pos_chosen_constrained = sample(pos_chosen_constrained, ceil(Int,(1-alpha_value)*length(pos_chosen_constrained)+epsilon))
        end
        all_pos=vcat(pos_chosen,pos_chosen_constrained)
        return all_pos[rand(1:length(all_pos))], true
    else
        return NewVisit(0,0,0,0,0.0,0.0, false, ToStoreVisit(SplitCosts(0,0,0,0,0,0),0,"","")), false
    end
end


function greedyrandomizedconstruction(inst::Instance, paramchosen::ChosenParameters, allparam::AllParameters, max_time)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT = inst
    sol = initializeSol(inst)
    new_visit, feasible = SelectNewVisitAllShips(inst, sol, paramchosen, allparam)
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
            new_visit, feasible = SelectNewVisitAllShips(inst, sol, paramchosen, allparam)
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
                                    new_visit, new_visit_constrained, feasible=SelectNewVisitPerShip(inst, sol, paramchosen, allparam, n, i, Pi[n][i], false)
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


function ChooseParamExp(allparam::AllParameters, N, type1, type2)
    @unpack Alpha, Proba, AverageCost, Nbexp, Q, adjustproba = allparam
    @unpack n_proba_typeship_distance, n_proba_typeship_cost, n_proba_all_count, n_proba_all_cost, n_proba_constrained, n_proba_local = adjustproba
    if type1=="both"
        tactictypeship = sample(["cost","dist"], Weights(Proba.ChooseTacticTypeShip))
    end
    if type1=="cost"
        tactictypeship="cost"
    end
    if type1=="dist"
        tactictypeship="dist"
    end
    if tactictypeship=="cost"
        IndexTypeShip = sample(1:n_proba_typeship_cost+1, Weights(Proba.CostTypeShip))
    else
        IndexTypeShip = sample(1:n_proba_typeship_distance+1, Weights(Proba.DistanceTypeShip))
    end


    if type2=="both"
        tacticallship = sample(["cost","count"], Weights(Proba.ChooseTacticAllShip))
    end
    if type2=="cost"
        tacticallship="cost"
    end
    if type2=="count"
        tacticallship="count"
    end
    if tacticallship=="cost"
        indexallship = sample(1:n_proba_all_cost+1, Weights(Proba.CostAllShip))
    else
        indexallship = sample(1:n_proba_all_count+1, Weights(Proba.CountAllShip))
    end
    indexconstrained = sample(1:n_proba_constrained+1, Weights(Proba.RateConstrained))

    indexlocalsearch = sample(1:n_proba_local+1, Weights(Proba.ChooseTacticLocalSearch))

    return ChosenParameters(tactictypeship,IndexTypeShip,tacticallship,indexallship, indexconstrained, indexlocalsearch)
end


function UpdateParametersExp(paramchosen::ChosenParameters, allparam::AllParameters, N, cost, new_cost)
    @unpack adjustproba = allparam
    @unpack n_proba_typeship_distance, n_proba_typeship_cost, n_proba_all_count, n_proba_all_cost, n_proba_constrained, n_proba_local = adjustproba
    newparam=deepcopy(allparam)
    if paramchosen.TacticTypeShip=="cost"
        nbexp = allparam.Nbexp.ChooseTacticTypeShip[1]
        newparam.AverageCost.ChooseTacticTypeShip[1] = (allparam.AverageCost.ChooseTacticTypeShip[1]*nbexp+new_cost)/(nbexp+1)
        newparam.Nbexp.ChooseTacticTypeShip[1]=nbexp+1
        for i in 1:2
            newav = deepcopy(newparam.AverageCost.ChooseTacticTypeShip[i])
            newparam.Q.ChooseTacticTypeShip[i]=cost^2/newav^2
        end
        for i in 1:2
            newq=deepcopy(newparam.Q.ChooseTacticTypeShip[i])
            newparam.Proba.ChooseTacticTypeShip[i] = newq/sum(newparam.Q.ChooseTacticTypeShip)
        end
        index = paramchosen.IndexTypeShip
        nbexp = allparam.Nbexp.CostTypeShip[index]
        newparam.AverageCost.CostTypeShip[index] = (allparam.AverageCost.CostTypeShip[index]*nbexp+new_cost)/(nbexp+1)
        newparam.Nbexp.CostTypeShip[index]=nbexp+1
        for i in 1:n_proba_typeship_cost+1
            newav = deepcopy(newparam.AverageCost.CostTypeShip[i])
            newparam.Q.CostTypeShip[i]=cost^2/newav
        end
        for i in 1:n_proba_typeship_cost+1
            newq=deepcopy(newparam.Q.CostTypeShip[i])
            newparam.Proba.CostTypeShip[i] = newq/sum(newparam.Q.CostTypeShip)
        end
    end
    if paramchosen.TacticTypeShip=="dist"
        nbexp = newparam.Nbexp.ChooseTacticTypeShip[2]
        newparam.AverageCost.ChooseTacticTypeShip[2] = (allparam.AverageCost.ChooseTacticTypeShip[2]*nbexp+new_cost)/(nbexp+1)
        newparam.Nbexp.ChooseTacticTypeShip[2]=nbexp+1
        for i in 1:2
            newav = deepcopy(newparam.AverageCost.ChooseTacticTypeShip[i])
            newparam.Q.ChooseTacticTypeShip[i]=cost^2/newav^2
        end
        for i in 1:2
            newq=deepcopy(newparam.Q.ChooseTacticTypeShip[i])
            newparam.Proba.ChooseTacticTypeShip[i] = newq/sum(newparam.Q.ChooseTacticTypeShip)
        end
        index = paramchosen.IndexTypeShip
        nbexp = allparam.Nbexp.DistanceTypeShip[index]
        newparam.AverageCost.DistanceTypeShip[index] = (allparam.AverageCost.DistanceTypeShip[index]*nbexp+new_cost)/(nbexp+1)
        newparam.Nbexp.DistanceTypeShip[index]=nbexp+1
        for i in 1:n_proba_typeship_distance+1
            newparam.Q.DistanceTypeShip[i]=cost^2/newparam.AverageCost.DistanceTypeShip[i]^2
        end
        for i in 1:n_proba_typeship_distance+1
            newparam.Proba.DistanceTypeShip[i] = newparam.Q.DistanceTypeShip[i]/sum(newparam.Q.DistanceTypeShip[n])
        end
    end
    
    if paramchosen.TacticAllShip=="cost"
        nbexp = allparam.Nbexp.ChooseTacticAllShip[1]
        newparam.AverageCost.ChooseTacticAllShip[1] = (allparam.AverageCost.ChooseTacticAllShip[1]*nbexp+new_cost)/(nbexp+1)
        newparam.Nbexp.ChooseTacticAllShip[1]=nbexp+1
        for i in 1:2
            newparam.Q.ChooseTacticAllShip[i]=cost/newparam.AverageCost.ChooseTacticAllShip[i]^2
        end
        for i in 1:2
            newparam.Proba.ChooseTacticAllShip[i] = newparam.Q.ChooseTacticAllShip[i]/sum(newparam.Q.ChooseTacticAllShip)
        end
        index = paramchosen.IndexAllShip
        nbexp = allparam.Nbexp.CostAllShip[index]
        newparam.AverageCost.CostAllShip[index] = (allparam.AverageCost.CostAllShip[index]*nbexp+new_cost)/(nbexp+1)
        newparam.Nbexp.CostAllShip[index]=nbexp+1
        for i in 1:n_proba_all_cost+1
            newparam.Q.CostAllShip[i]=cost^2/newparam.AverageCost.CostAllShip[i]^2
        end
        for i in 1:n_proba_all_cost+1
            newparam.Proba.CostAllShip[i] = newparam.Q.CostAllShip[i]/sum(newparam.Q.CostAllShip)
        end
    end
    if paramchosen.TacticAllShip=="count"
        nbexp = allparam.Nbexp.ChooseTacticAllShip[2]
        newparam.AverageCost.ChooseTacticAllShip[2] = (allparam.AverageCost.ChooseTacticAllShip[2]*nbexp+new_cost)/(nbexp+1)
        newparam.Nbexp.ChooseTacticAllShip[2]=nbexp+1
        for i in 1:2
            newparam.Q.ChooseTacticAllShip[i]=cost^2/newparam.AverageCost.ChooseTacticAllShip[i]^2
        end
        for i in 1:2
            newparam.Proba.ChooseTacticAllShip[i] = newparam.Q.ChooseTacticAllShip[i]/sum(newparam.Q.ChooseTacticAllShip)
        end
        index = paramchosen.IndexAllShip
        nbexp = allparam.Nbexp.CountAllShip[index]
        newparam.AverageCost.CountAllShip[index] = (allparam.AverageCost.CountAllShip[index]*nbexp+new_cost)/(nbexp+1)
        newparam.Nbexp.CountAllShip[index]=nbexp+1
        for i in 1:n_proba_all_count+1
            newparam.Q.CountAllShip[i]=cost^2/newparam.AverageCost.CountAllShip[i]^2
        end
        for i in 1:n_proba_all_count+1
            newparam.Proba.CountAllShip[i] = newparam.Q.CountAllShip[i]/sum(newparam.Q.CountAllShip)
        end
    end

    index = paramchosen.IndexRateConstrained
    nbexp = allparam.Nbexp.RateConstrained[index]
    newparam.AverageCost.RateConstrained[index] = (allparam.AverageCost.RateConstrained[index]*nbexp+new_cost)/(nbexp+1)
    newparam.Nbexp.RateConstrained[index]=nbexp+1
    for i in 1:n_proba_constrained+1
        newparam.Q.RateConstrained[i]=cost^2/newparam.AverageCost.RateConstrained[i]^2
    end
    for i in 1:n_proba_constrained+1
        newparam.Proba.RateConstrained[i] = newparam.Q.RateConstrained[i]/sum(newparam.Q.RateConstrained)
    end

    index = paramchosen.ChooseTacticLocalSearch
    nbexp = allparam.Nbexp.ChooseTacticLocalSearch[index]
    newparam.AverageCost.ChooseTacticLocalSearch[index] = (allparam.AverageCost.ChooseTacticLocalSearch[index]*nbexp+new_cost)/(nbexp+1)
    newparam.Nbexp.ChooseTacticLocalSearch[index]=nbexp+1
    for i in 1:n_proba_local+1
        newparam.Q.ChooseTacticLocalSearch[i]=cost^2/newparam.AverageCost.ChooseTacticLocalSearch[i]^2
    end
    for i in 1:n_proba_local+1
        newparam.Proba.ChooseTacticLocalSearch[i] = newparam.Q.ChooseTacticLocalSearch[i]/sum(newparam.Q.ChooseTacticLocalSearch)
    end



    return newparam
end


function UpdateParameters(paramchosen::ChosenParameters, allparam::AllParameters, N, cost, new_cost)
    @unpack adjustproba = allparam
    @unpack n_proba_typeship_distance, n_proba_typeship_cost, n_proba_all_count, n_proba_all_cost, n_proba_constrained, n_proba_local = adjustproba
    newparam=deepcopy(allparam)
    if paramchosen.TacticTypeShip=="cost"
        nbexp = allparam.Nbexp.ChooseTacticTypeShip[1]
        newparam.AverageCost.ChooseTacticTypeShip[1] = (allparam.AverageCost.ChooseTacticTypeShip[1]*nbexp+new_cost)/(nbexp+1)
        newparam.Nbexp.ChooseTacticTypeShip[1]=nbexp+1
        for i in 1:2
            newav = deepcopy(newparam.AverageCost.ChooseTacticTypeShip[i])
            newparam.Q.ChooseTacticTypeShip[i]=cost/newav
        end
        for i in 1:2
            newq=deepcopy(newparam.Q.ChooseTacticTypeShip[i])
            newparam.Proba.ChooseTacticTypeShip[i] = newq/sum(newparam.Q.ChooseTacticTypeShip)
        end
        index = paramchosen.IndexTypeShip
        nbexp = allparam.Nbexp.CostTypeShip[index]
        newparam.AverageCost.CostTypeShip[index] = (allparam.AverageCost.CostTypeShip[index]*nbexp+new_cost)/(nbexp+1)
        newparam.Nbexp.CostTypeShip[index]=nbexp+1
        for i in 1:n_proba_typeship_cost+1
            newav = deepcopy(newparam.AverageCost.CostTypeShip[i])
            newparam.Q.CostTypeShip[i]=cost/newav
        end
        for i in 1:n_proba_typeship_cost+1
            newq=deepcopy(newparam.Q.CostTypeShip[i])
            newparam.Proba.CostTypeShip[i] = newq/sum(newparam.Q.CostTypeShip)
        end
    end
    if paramchosen.TacticTypeShip=="dist"
        nbexp = newparam.Nbexp.ChooseTacticTypeShip[2]
        newparam.AverageCost.ChooseTacticTypeShip[2] = (allparam.AverageCost.ChooseTacticTypeShip[2]*nbexp+new_cost)/(nbexp+1)
        newparam.Nbexp.ChooseTacticTypeShip[2]=nbexp+1
        for i in 1:2
            newav = deepcopy(newparam.AverageCost.ChooseTacticTypeShip[i])
            newparam.Q.ChooseTacticTypeShip[i]=cost/newav
        end
        for i in 1:2
            newq=deepcopy(newparam.Q.ChooseTacticTypeShip[i])
            newparam.Proba.ChooseTacticTypeShip[i] = newq/sum(newparam.Q.ChooseTacticTypeShip)
        end
        index = paramchosen.IndexTypeShip
        nbexp = allparam.Nbexp.DistanceTypeShip[index]
        newparam.AverageCost.DistanceTypeShip[index] = (allparam.AverageCost.DistanceTypeShip[index]*nbexp+new_cost)/(nbexp+1)
        newparam.Nbexp.DistanceTypeShip[index]=nbexp+1
        for i in 1:n_proba_typeship_distance+1
            newparam.Q.DistanceTypeShip[i]=cost/newparam.AverageCost.DistanceTypeShip[i]
        end
        for i in 1:n_proba_typeship_distance+1
            newparam.Proba.DistanceTypeShip[i] = newparam.Q.DistanceTypeShip[i]/sum(newparam.Q.DistanceTypeShip)
        end
    end
    if paramchosen.TacticAllShip=="cost"
        nbexp = allparam.Nbexp.ChooseTacticAllShip[1]
        newparam.AverageCost.ChooseTacticAllShip[1] = (allparam.AverageCost.ChooseTacticAllShip[1]*nbexp+new_cost)/(nbexp+1)
        newparam.Nbexp.ChooseTacticAllShip[1]=nbexp+1
        for i in 1:2
            newparam.Q.ChooseTacticAllShip[i]=cost/newparam.AverageCost.ChooseTacticAllShip[i]
        end
        for i in 1:2
            newparam.Proba.ChooseTacticAllShip[i] = newparam.Q.ChooseTacticAllShip[i]/sum(newparam.Q.ChooseTacticAllShip)
        end
        index = paramchosen.IndexAllShip
        nbexp = allparam.Nbexp.CostAllShip[index]
        newparam.AverageCost.CostAllShip[index] = (allparam.AverageCost.CostAllShip[index]*nbexp+new_cost)/(nbexp+1)
        newparam.Nbexp.CostAllShip[index]=nbexp+1
        for i in 1:n_proba_all_cost+1
            newparam.Q.CostAllShip[i]=cost/newparam.AverageCost.CostAllShip[i]
        end
        for i in 1:n_proba_all_cost+1
            newparam.Proba.CostAllShip[i] = newparam.Q.CostAllShip[i]/sum(newparam.Q.CostAllShip)
        end
    end
    if paramchosen.TacticAllShip=="count"
        nbexp = allparam.Nbexp.ChooseTacticAllShip[2]
        newparam.AverageCost.ChooseTacticAllShip[2] = (allparam.AverageCost.ChooseTacticAllShip[2]*nbexp+new_cost)/(nbexp+1)
        newparam.Nbexp.ChooseTacticAllShip[2]=nbexp+1
        for i in 1:2
            newparam.Q.ChooseTacticAllShip[i]=cost/newparam.AverageCost.ChooseTacticAllShip[i]
        end
        for i in 1:2
            newparam.Proba.ChooseTacticAllShip[i] = newparam.Q.ChooseTacticAllShip[i]/sum(newparam.Q.ChooseTacticAllShip)
        end
        index = paramchosen.IndexAllShip
        nbexp = allparam.Nbexp.CountAllShip[index]
        newparam.AverageCost.CountAllShip[index] = (allparam.AverageCost.CountAllShip[index]*nbexp+new_cost)/(nbexp+1)
        newparam.Nbexp.CountAllShip[index]=nbexp+1
        for i in 1:n_proba_all_count+1
            newparam.Q.CountAllShip[i]=cost/newparam.AverageCost.CountAllShip[i]
        end
        for i in 1:n_proba_all_count+1
            newparam.Proba.CountAllShip[i] = newparam.Q.CountAllShip[i]/sum(newparam.Q.CountAllShip)
        end
    end
    index = paramchosen.IndexRateConstrained
    nbexp = allparam.Nbexp.RateConstrained[index]
    newparam.AverageCost.RateConstrained[index] = (allparam.AverageCost.RateConstrained[index]*nbexp+new_cost)/(nbexp+1)
    newparam.Nbexp.RateConstrained[index]=nbexp+1
    for i in 1:n_proba_constrained+1
        newparam.Q.RateConstrained[i]=cost/newparam.AverageCost.RateConstrained[i]
    end
    for i in 1:n_proba_constrained+1
        newparam.Proba.RateConstrained[i] = newparam.Q.RateConstrained[i]/sum(newparam.Q.RateConstrained)
    end

    index = paramchosen.ChooseTacticLocalSearch
    nbexp = allparam.Nbexp.ChooseTacticLocalSearch[index]
    newparam.AverageCost.ChooseTacticLocalSearch[index] = (allparam.AverageCost.ChooseTacticLocalSearch[index]*nbexp+new_cost)/(nbexp+1)
    newparam.Nbexp.ChooseTacticLocalSearch[index]=nbexp+1
    for i in 1:n_proba_local+1
        newparam.Q.ChooseTacticLocalSearch[i]=cost/newparam.AverageCost.ChooseTacticLocalSearch[i]
    end
    for i in 1:n_proba_local+1
        newparam.Proba.ChooseTacticLocalSearch[i] = newparam.Q.ChooseTacticLocalSearch[i]/sum(newparam.Q.ChooseTacticLocalSearch)
    end


    return newparam
end



function GRASP_reactive(seed,N,Nout,qli, type1, type2, adjustproba, alphaboat, alpharandom, time_local, max_time_heur, max_time, expname)
    inst = readInstFromFile("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_$seed"*"_$N"*"_$Nout"*"_$qli"*".txt")
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp = inst
    cost=1000000000
    worst_cost=1000000000
    sol=initializeSol(inst)
    allparam = initializeParam(inst, adjustproba)
    start = time_ns()
    elapsed = round((time_ns()-start)/1e9,digits=3)
    start_iter=time_ns()
    elapsed_iter = round((time_ns()-start_iter)/1e9,digits=3)
    nb_iter=0
    while elapsed<max_time
        paramchosen = ChooseParamExp(allparam, N, type1, type2)
        start_heur = time_ns()
        new_sol = greedyrandomizedconstruction(inst, paramchosen, allparam, max_time_heur)

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

            d=prepareSolIter(seed,N,Nout,qli,nb_iter,inst, new_sol, new_cost_heur, expname)
            CSV.write("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/reactiveGRASP/countcost_nostage/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli"*"/iter_$nb_iter"*"_beforelocal"*".csv", d)
            
            start_local=time_ns()
            new_sol, new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost = local_search(inst, deepcopy(new_sol), ceil(Int,cost), allparam, paramchosen, alphaboat, alpharandom, time_local)
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

            d = prepareSolIter(seed,N,Nout,qli,nb_iter,inst, new_sol, cost, expname)
            CSV.write("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/reactiveGRASP/countcost_nostage/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli"*"/iter_$nb_iter"*".csv", d)
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
            allparam = UpdateParameters(paramchosen, allparam, N, cost, new_cost)    
            #allparam = UpdateParametersExp(paramchosen, allparam, N, cost, new_cost)                
        else
            allparam = UpdateParameters(paramchosen, allparam, N, cost, worst_cost) 
            #allparam = UpdateParametersExp(paramchosen, allparam, N, cost, worst_cost) 
        end
        elapsed = round((time_ns()-start)/1e9,digits=3)
    end
    #print('\n')
    #print("Cost at the end of heur")
    #print('\n')
    #print(cost)
    d = prepareSolIter(seed,N,Nout,qli,nb_iter,inst, sol, cost, expname)
    nb_iter+=1
    CSV.write("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/reactiveGRASP/countcost_nostage/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli"*"/iter_$nb_iter"*".csv", d)
    return sol, cost, allparam
end

#inst = readInstFromFile("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_1_10_5_80.txt")
#GRASP_reactive(inst::Instance, type1, type2, adjustproba, alphaboat, alpharandom, time_local, max_time_heur, max_time, expname)
#CSV.write("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/reactiveGRASP/sols/HEUR_LOCAL_exp1_sol_$seed"*"_$N"*"_$Nout"*"_$qli"*".csv", d)
#expname="exp2"
#seed=1
#N=10
#Nout=5
#qli=80
#if isdir("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/reactiveGRASP/countcost_nostage/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli")==false
#    mkdir("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/reactiveGRASP/countcost_nostage/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli")
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