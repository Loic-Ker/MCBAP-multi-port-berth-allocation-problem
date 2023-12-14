using UnPack 
using Statistics
import XLSX
using CSV, Tables
using DataFrames
using Random
using StatsBase
using Distributions
using JuMP, CPLEX
include("../../MBAP_INST.jl")
include("../toolsMatrixTimes.jl")
include("../check_solution.jl")
include("../utilInit.jl")



################################ The Local search method

#################### Here these two functions are doing the same two construction steps that the GRASP ; select a position and time for each visit and then select a visit

## The first step   
function SelectNewVisitPerShipLocalSearch(inst::Instance, sol::Sol, n, c, p, paramfixed)
    @unpack N, h, qli, shipsIn, Pi, Bp, delta, dist, T, maxC, maxT, Pc = inst
    @unpack M = sol
    l = ceil(Int, shipsIn[n].l/qli)
    new_visits_constrained = Vector{Tuple}()
    t1=sol.visits[n][c].minT
    t2=sol.visits[n][c].maxT    
    bestPos=shipsIn[n].Bi[c]/qli
    if paramfixed.lookforconstrained 
        pos = findConstrainedPos(inst, sol, p, t1, t2, n, c, l)
        for (b,t) in pos
            hand = ceil(Int, h[n][c][b])
            if t + hand <= t2
                if M[n][c][b,t+1]
                    this_cost, delay_cost, waiting_cost, penalty, handling_cost, fuel_cost, feas = computeCostPosSol(inst, n, c, b, t,sol)
                    time=penalty
                    if feas
                        distance = abs(bestPos-b)/l
                        push!(new_visits_constrained, (NewVisit(n,c,b,t,this_cost,distance,time, true, ToStoreVisit(SplitCosts(ceil(Int, this_cost), ceil(Int,delay_cost), ceil(Int,waiting_cost), ceil(Int,penalty), ceil(Int,handling_cost), ceil(Int,fuel_cost)),0,"","")), t+hand, distance,this_cost))
                    end
                end
            end
        end
        new_visits_constrained = sort(new_visits_constrained, by = x -> x[2])
    else
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
                        time=deepcopy(penalty)
                        if feas
                            distance = abs(bestPos-b)/l
                            find_it=true
                            push!(new_visits_constrained, (NewVisit(n,c,b,t,this_cost,distance,time, true, ToStoreVisit(SplitCosts(ceil(Int, this_cost), ceil(Int,delay_cost), ceil(Int,waiting_cost), ceil(Int,penalty), ceil(Int,handling_cost), ceil(Int,fuel_cost)),0,"","")), t+hand, distance, this_cost))
                        end
                    end
                end
            end
        end
    end
    if paramfixed.LocalSearchOne=="time"
        new_visits_constrained = sort(new_visits_constrained, by = x -> x[2])
    end
    if paramfixed.LocalSearchOne=="dist"
        new_visits_constrained = sort(new_visits_constrained, by = x -> x[3])
    end
    if paramfixed.LocalSearchOne=="cost"
        new_visits_constrained = sort(new_visits_constrained, by = x -> x[4])
    end
    if length(new_visits_constrained)>0
        return new_visits_constrained[1], true
    else
        return NewVisit(0,0,0,0,0.0,0.0,0.0, false, ToStoreVisit(SplitCosts(0,0,0,0,0,0),0,"","")), false
    end
end

## And then the second
function SelectNewVisitLocalSearch(inst::Instance, sol::Sol, list_visits::Vector{Tuple}, paramchosen, paramfixed)
    @unpack N, h, qli, shipsIn, Pi, Bp, delta, dist, T, maxC, maxT = inst
    @unpack M = sol
    @unpack Alpha, Proba, AverageCost, Nbexp, Q = allparam
    new_visits = Vector{NewVisit}()
    new_visits_constrained = Vector{NewVisit}()
    new_visits_constrained = Vector{Tuple}()
    for (n,c) in list_visits
        sch = sol.visits[n]
        if sch[c].planned == false
            visit_constrained, feasible=SelectNewVisitPerShipLocalSearch(inst, sol, n, c, Pi[n][c], paramfixed)
            if feasible
                push!(new_visits_constrained,visit_constrained)
            end
        end
    end

    if paramfixed.LocalSearchAll=="cost"
        new_visits_constrained = sort(new_visits_constrained,  by = x -> x[4])
    end
    if paramfixed.LocalSearchAll=="dist"
        new_visits_constrained = sort(new_visits_constrained,  by = x -> x[3])
    end
    if paramfixed.LocalSearchAll=="time"
        new_visits_constrained = sort(new_visits_constrained,  by = x -> x[2])
    end

    if length(new_visits_constrained)>0
        return new_visits_constrained[1][1] , true
    else
        return NewVisit(0,0,0,0,0.0,0.0,0.0, false, ToStoreVisit(SplitCosts(0,0,0,0,0,0),0,"","")), false
    end
end


## And then we use these two steps to construct completly the solution 
function replaceFromList(inst::Instance, this_sol::Sol, allparam::AllParameters, paramchosen, boats_to_be_visited)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT = inst
    sol=deepcopy(this_sol)
    continue_=true
    
    list_to_be_visited = Vector{Tuple}()
    list_to_be_visited_first = Vector{Tuple}()
    for n in boats_to_be_visited
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
    new_visit, feasible = SelectNewVisitLocalSearch(inst, sol, list_to_be_visited_first, paramchosen, paramfixed)
    while length(list_to_be_visited)>0
        if feasible
            @unpack n,c,b,t,cost,distance,constrained,store= new_visit
            deleteat!(list_to_be_visited, findall(x->x==(n,c),list_to_be_visited))
            deleteat!(list_to_be_visited_first, findall(x->x==(n,c),list_to_be_visited_first))
            if c<length(Pi[n])
                push!(list_to_be_visited_first, (n,c+1))
            end
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
            for n in boats_to_be_visited
                for (c,p) in enumerate(Pi[n])
                    if sol.visits[n][c].planned == false
                        continue_ = true
                    end
                end
            end


            if continue_
                if length(list_to_be_visited_first)>0
                    new_visit, feasible = SelectNewVisitLocalSearch(inst, sol, list_to_be_visited_first, paramchosen, paramfixed)
                else
                    new_visit, feasible = SelectNewVisitLocalSearch(inst, sol, list_to_be_visited, paramchosen, paramfixed)
                end
                if feasible==false
                    return initializeSol(inst, allparam), false
                end
            end
        else
            return initializeSol(inst, allparam), false
        end
    end
    return sol, true
end




#################### Removing the visits we want to remove a cluster of visits closed to each other given a random boat
function get_closestlocal(inst::Instance, sol::Sol, randomBoat, alpha)
    @unpack N, Ntot, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT, Nl, gamma, Hc, Dc, Fc, Ic, Pc, beta, ports = inst
    all_visits_distance=[]
    closest_visits=[]
    for c_boat in 1:length(inst.Pi[randomBoat])
        for n in 1:N
            for c in 1:length(inst.Pi[n])
                if inst.Pi[n][c] == inst.Pi[randomBoat][c_boat] && n!=randomBoat
                    distance_squares=abs(sol.visits[n][c].b-sol.visits[randomBoat][findfirst(item -> item == inst.Pi[n][c], inst.Pi[randomBoat])].b)+abs(sol.visits[n][c].t-sol.visits[randomBoat][findfirst(item -> item == inst.Pi[n][c], inst.Pi[randomBoat])].t)
                    append!(all_visits_distance,[(distance_squares,[n,c])])
                end
            end
        end
        closest=sort(all_visits_distance, by = first)[1:min(alpha,length(all_visits_distance))]
        for el in closest
            append!(closest_visits,[el[2]])
        end
    end
    for c in 1:length(inst.Pi[randomBoat])
        append!(closest_visits,[[randomBoat,c]])
    end
    return closest_visits
end




#################### Here wee define a function to remove some visits and then replace them using the functions above
function localSearchRemovalReplace(inst::Instance, this_sol::Sol, allparam, paramfixed, paramchosen)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT = inst
    sol=deepcopy(this_sol)
    nb_to_remove = ceil(Int, paramfixed.LocalSearchBoat*N)
    start = time_ns()
    boats_to_remove = Vector{Tuple}()
    for n in 1:N
        all_penalty_this_boat = Vector{}()
        for (c,p) in enumerate(Pi[n])
            push!(all_penalty_this_boat, sol.visits[n][c].store.cost.all - sol.visits[n][c].store.cost.fuel)
        end
        push!(boats_to_remove, (n,mean(all_penalty_this_boat)))
    end

    boats_to_remove = sort(boats_to_remove, by= x-> x[2], rev=true)
    boats_to_remove = boats_to_remove[1:min(length(boats_to_remove),nb_to_remove)]
    boats_to_be_visited = Vector{}()
    for (n,penalty) in boats_to_remove
        for (c,p) in enumerate(Pi[n])
            sol.visits[n][c].p = -1
            sol.visits[n][c].b = -1
            sol.visits[n][c].t = -1
            sol.visits[n][c].planned = false
            sol.visits[n][c].minT = max(shipsIn[n].sT[c], T[n,c,1])
            sol.visits[n][c].maxT = maxT
            push!(boats_to_be_visited, n)
        end
    end

    random_removals = Vector{}()
    for n in 1:N
        possible=false
        for (c,p) in enumerate(Pi[n])
            if sol.visits[n][c].planned
                possible=true
            end
        end
        if possible
            push!(random_removals, n)
        end
    end
    random_removals = sample(random_removals, ceil(Int, paramfixed.LocalSearchRandom*length(random_removals)); replace=false)#ceil(Int, rand(Uniform(paramfixed.LocalSearchRandomMin,paramfixed.LocalSearchRandomMax))*length(random_removals)); replace=false)

    for n in random_removals
        for (c,p) in enumerate(Pi[n])
            sol.visits[n][c].p = -1
            sol.visits[n][c].b = -1
            sol.visits[n][c].t = -1
            sol.visits[n][c].planned = false
            sol.visits[n][c].minT = max(shipsIn[n].sT[c], T[n,c,1])
            sol.visits[n][c].maxT = maxT
        end
    end
    start = time_ns()
    sol, feasible = replaceFromList(inst, sol, allparam, paramchosen, boats_to_be_visited)
    if feasible
        sol, feasible = replaceFromList(inst, sol, allparam, paramchosen, random_removals)
        elapsed = round((time_ns()-start)/1e9,digits=3)
        #print('\n')
        #println("Local Search second step: $(elapsed) seconds")
        if feasible
            return sol, true
        else
            return initializeSol(inst, allparam), false
        end
    else
        return initializeSol(inst, allparam), false
    end
end


#################### Finally we do this step in a loop as long as we don't reach the maximum local search time
function manualLocalSearch(inst::Instance, this_sol::Sol, cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost, allparam, paramfixed, paramchosen, timelocal, bestcost)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp = inst
    start = time_ns()
    elapsed = round((time_ns()-start)/1e9,digits=3)
    sol=deepcopy(this_sol)
    old_cost=deepcopy(cost)
    old_waiting_cost=deepcopy(waiting_cost)
    old_delay_cost=deepcopy(delay_cost)
    old_penalty_cost=deepcopy(penalty_cost)
    old_handling_cost=deepcopy(handling_cost)
    old_fuel_cost=deepcopy(fuel_cost)
    count_total=0
    while elapsed<timelocal
        start_step = time_ns()
        new_sol, feasible = localSearchRemovalReplace(inst, sol, allparam, paramfixed, paramchosen)
        feasible1=true
        for n in 1:N
            # The times :
            for (c,p) in enumerate(inst.Pi[n])
                if new_sol.visits[n][c].planned == false
                    feasible1 = false              
                end
            end
        end
        if feasible && feasible1 && checkSolutionFeasability(inst, new_sol)
            count_total=count_total+1
            new_cost1, delay_cost1, waiting_cost1, penalty_cost1, handling_cost1, fuel_cost1=checkSolutionCost(inst, new_sol)
            new_sol1, new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost = pushTime(inst, new_sol, paramfixed, time_local)
            if new_cost1>new_cost
                new_sol.pushimprove=1
                new_sol=deepcopy(new_sol1)
                new_cost=deepcopy(new_cost1)
                delay_cost=deepcopy(delay_cost1)
                waiting_cost=deepcopy(waiting_cost1)
                penalty_cost=deepcopy(penalty_cost1)
                handling_cost=deepcopy(handling_cost1)
                fuel_cost=deepcopy(fuel_cost1)
            end
            if new_cost<old_cost
                sol=deepcopy(new_sol)
                old_cost=deepcopy(new_cost)
                old_delay_cost=deepcopy(delay_cost)
                old_waiting_cost=deepcopy(waiting_cost)
                old_penalty_cost=deepcopy(penalty_cost)
                old_handling_cost=deepcopy(handling_cost)
                old_fuel_cost=deepcopy(fuel_cost)
            end
            #allparam = UpdateAfterLocalParameters(paramchosen, allparam, cost, bestcost)
            #paramchosen = ChooseAfterLocal(allparam, paramchosen, paramfixed)
        end
        elapsed = round((time_ns()-start)/1e9,digits=3)
        sol.countlocal=count_total
    end
    return sol, old_cost, old_delay_cost, old_waiting_cost, old_penalty_cost, old_handling_cost, old_fuel_cost
end



################################ The pushing method

#################### Two functions to get the minimum time available for a visit

## Here without any constraint (whatever position and time available)
function getSpecificAvailableTimes(inst::Instance, sol::Sol, listnc::Vector{Tuple})
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT, qli = inst
    pos_time = Vector{Tuple}()
    for (n,c) in listnc
        this_sol_time = sol.visits[n][c].t+ceil(Int, h[n][c][sol.visits[n][c].b])
        t1=sol.visits[n][c].minT
        l = ceil(Int, shipsIn[n].l/qli)
        ## Check all the possible solutions
        notfound=false
        for b in 1:Bp[Pi[n][c]]-l
            hand = ceil(Int, h[n][c][b])
            while t1<=this_sol_time && notfound==false
                if sol.M[n][c][b,t1+1]
                    notfound = true
                    push!(pos_time,(n,c,b,t1,t1+hand))
                end
                t1+=1
            end
        end
        if notfound==false
            push!(pos_time,(n,c,sol.visits[n][c].b,sol.visits[n][c].t,this_sol_time))
        end
    end
    pos_time = sort(pos_time, by = x -> x[5])
    return pos_time
end

## Here it's the same function but we add the constraint that the ship must be close to another ship in the 2D graph
function getSpecificConstrainedAvailableTimes(inst::Instance, sol::Sol, listnc::Vector{Tuple}, paramfixed)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT, qli = inst
    pos_time = Vector{Tuple}()
    for (n,c) in listnc
        this_sol_time = sol.visits[n][c].t+ceil(Int, h[n][c][sol.visits[n][c].b])
        t1=sol.visits[n][c].minT
        l = ceil(Int, shipsIn[n].l/qli)
        ## Check all the possible solutions
        notfound=false
        for b in 1:Bp[Pi[n][c]]-l
            hand = ceil(Int, h[n][c][b])
            while t1<=this_sol_time && notfound==false
                if sol.M[n][c][b,t1+1]
                    if paramfixed.PushOnlyConstrained
                        ## We chek if we place it at the borders of the berth
                        if b==1 || b>=Bp[Pi[n][c]]-l-1 
                            notfoud = true
                            push!(pos_time,(n,c,b,t1,t1+hand))
                        else
                            if false in sol.M[n][c][b+l:b+l+1,t1+1:t1+hand]  || false in sol.M[n][c][max(b-1,1):b,t1+1:t1+hand]
                                notfoud = true
                                push!(pos_time,(n,c,b,t1,t1+hand))
                            end
                        end
                    else
                        notfoud = true
                        push!(pos_time,(n,c,b,t1,t1+hand))
                    end
                end
                t1+=1
            end
        end
        if notfound==false
            push!(pos_time,(n,c,sol.visits[n][c].b,sol.visits[n][c].t,this_sol_time))
        end
    end
    pos_time = sort(pos_time, by = x -> x[5])
    return pos_time
end



#################### Function to squeeze all the visits together
## Here we squeeze all the visits until we reach the maximum local search time

function pushTime(inst::Instance, sol::Sol, paramfixed, timelocal)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT = inst
    listnc = Vector{Tuple}()
    cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost=checkSolutionCost(inst, sol)
    new_sol=deepcopy(sol)
    for n in 1:N
        for c in 1:length(Pi[n])
            push!(listnc,(n,c))
        end
    end
    start = time_ns()
    elapsed = round((time_ns()-start)/1e9,digits=3)
    first=true
    while length(listnc)>0 && elapsed<timelocal
        if first
            pos_time = getSpecificAvailableTimes(inst, new_sol, listnc)
            first=false
        else
            if paramfixed.PushOnlyConstrained
                pos_time = getSpecificConstrainedAvailableTimes(inst, new_sol, listnc, paramfixed)
            else
                pos_time = getSpecificAvailableTimes(inst, new_sol, listnc)
            end
        end
        chosennc = pos_time[1]
        n = chosennc[1]
        c = chosennc[2]
        b = chosennc[3]
        t = chosennc[4]
        l = ceil(Int, shipsIn[n].l/qli)
        hand = ceil(Int, h[n][c][chosennc[3]])
        @unpack M, visits = new_sol
        new_sol.visits = updateTimesAfterVisit(inst, visits, n, c, b, t)
        new_sol.M = updateMpositions(inst, n, c, b, t, l, hand, M)
        new_sol.visits[n][c].p = Pi[n][c]
        new_sol.visits[n][c].b = b
        new_sol.visits[n][c].t = t
        deleteat!(listnc,findall(x->x[1]==n && x[2]==c, listnc))
        elapsed = round((time_ns()-start)/1e9,digits=3)
    end
    if checkSolutionFeasability(inst, new_sol)
        new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost=checkSolutionCost(inst, new_sol)
        if new_cost<cost
            return new_sol, new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost
        else
            return sol, cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost
        end
    else
        return sol, cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost
    end
end
