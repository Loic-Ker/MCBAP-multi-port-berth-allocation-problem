using UnPack 
using Statistics
import XLSX
using CSV, Tables
using DataFrames
using Random
using StatsBase
using Distributions
## Packages needed for solving MIPs
using JuMP, CPLEX
#m = Model(GLPK.Optimizer)
#m = Model(HiGHS.Optimizer)
include("../../MBAP_INST.jl")
include("MBAP_SOL.jl")
include("toolsMatrixTimes.jl")
include("check_solution.jl")


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

function getOnevisitPushingTime(inst::Instance, sol::Sol, n,c,tactic)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT, qli = inst
    pos_time = Vector{Tuple}()
    if tactic
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
    else
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
    end
    if notfound==false
        push!(pos_time,(n,c,sol.visits[n][c].b,sol.visits[n][c].t,this_sol_time))
    end
    pos_time = sort(pos_time, by = x -> x[5])
    return pos_time
end

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
    #print('\n')
    #print("We start the local search")
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
        #if listnc==[]
        #    if checkSolutionFeasability(inst, new_sol)
        #        new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost=checkSolutionCost(inst, new_sol)
        #        if new_cost<cost
        #            cost=new_cost
        #            sol=new_sol
        #            new_sol = deepcopy(sol)
        #        end
        #    end
        #    for n in 1:N
        #        for c in 1:length(Pi[n])
        #            push!(listnc,(n,c))
        #        end
        #    end
        #end
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

function pushTime2(inst::Instance, sol::Sol, paramfixed, timelocal)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT = inst
    listnc = Vector{Tuple}()
    cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost=checkSolutionCost(inst, sol)
    new_sol=deepcopy(sol)
    for n in 1:N
        for c in 1:length(Pi[n])
            push!(listnc,(n,c,sol.visits[n][c].t))
        end
    end
    listnc=sort(listnc, by = x -> x[3])
    #print('\n')
    #print("We start the local search")
    start = time_ns()
    elapsed = round((time_ns()-start)/1e9,digits=3)
    first=true
    for (n,c,time) in listnc
        if first
            pos_time = getOnevisitPushingTime(inst, new_sol, n,c,true)
            first=false
        else
            pos_time = getOnevisitPushingTime(inst, new_sol, n,c,paramfixed.PushOnlyConstrained)
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
        #if listnc==[]
        #    if checkSolutionFeasability(inst, new_sol)
        #        new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost=checkSolutionCost(inst, new_sol)
        #        if new_cost<cost
        #            cost=new_cost
        #            sol=new_sol
        #            new_sol = deepcopy(sol)
        #        end
        #    end
        #    for n in 1:N
        #        for c in 1:length(Pi[n])
        #            push!(listnc,(n,c))
        #        end
        #    end
        #end
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





##################################################################################################################################
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

function SelectNewVisitPerShipLocalSearch(inst::Instance, sol::Sol, n, c, p)
    @unpack N, h, qli, shipsIn, Pi, Bp, delta, dist, T, maxC, maxT, Pc = inst
    @unpack M = sol
    l = ceil(Int, shipsIn[n].l/qli)
    new_visits_constrained = Vector{Tuple}()
    t1=sol.visits[n][c].minT
    t2=sol.visits[n][c].maxT    
    bestPos=shipsIn[n].Bi[c]/qli
    only_constrained=true
    if only_constrained==true
        pos = findConstrainedPos(inst, sol, p, t1, t2, n, c, l)
        for (b,t) in pos
            hand = ceil(Int, h[n][c][b])
            if t + hand <= t2
                if M[n][c][b,t+1]
                    this_cost, delay_cost, waiting_cost, penalty, handling_cost, fuel_cost, feas = computeCostPosSol(inst, n, c, b, t,sol)
                    time=penalty
                    if feas
                        distance = abs(bestPos-b)/l
                        push!(new_visits_constrained, (NewVisit(n,c,b,t,this_cost,distance,time, true, ToStoreVisit(SplitCosts(ceil(Int, this_cost), ceil(Int,delay_cost), ceil(Int,waiting_cost), ceil(Int,penalty), ceil(Int,handling_cost), ceil(Int,fuel_cost)),0,"","")), t+hand, this_cost))
                    end
                end
            end
        end
    end
    new_visits_constrained = sort(new_visits_constrained, by = x -> x[2])
    if length(new_visits_constrained)<=0
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
                        time=penalty
                        if feas
                            distance = abs(bestPos-b)/l
                            find_it=true
                            push!(new_visits_constrained, (NewVisit(n,c,b,t,this_cost,distance,time, true, ToStoreVisit(SplitCosts(ceil(Int, this_cost), ceil(Int,delay_cost), ceil(Int,waiting_cost), ceil(Int,penalty), ceil(Int,handling_cost), ceil(Int,fuel_cost)),0,"","")), t+hand, this_cost))
                        end
                    end
                end
            end
        end
    end
    new_visits_constrained = sort(new_visits_constrained, by = x -> x[2])
    if length(new_visits_constrained)>0
        return new_visits_constrained[1], true
    else
        return NewVisit(0,0,0,0,0.0,0.0,0.0, false, ToStoreVisit(SplitCosts(0,0,0,0,0,0),0,"","")), false
    end
end

function SelectNewVisitLocalSearch(inst::Instance, sol::Sol, list_visits::Vector{Tuple}) #, placed::Vector{Tuple{Int64,Int64}})
    #TODO: t1 cannot be greater than t2!
    # adapt to given planned previous visits for each ship
    @unpack N, h, qli, shipsIn, Pi, Bp, delta, dist, T, maxC, maxT = inst
    @unpack M = sol
    @unpack Alpha, Proba, AverageCost, Nbexp, Q = allparam
    # C = maximum(length.(Pi))
    new_visits = Vector{NewVisit}()
    new_visits_constrained = Vector{NewVisit}()
    count_available=Dict()
    for n in 1:N
        count_available[n]=Dict()
    end
    new_visits_constrained = Vector{Tuple}()
    for (n,c) in list_visits
        sch = sol.visits[n]
        if sch[c].planned == false
            visit_constrained, feasible=SelectNewVisitPerShipLocalSearch(inst, sol, n, c, Pi[n][c])
            if feasible
                push!(new_visits_constrained,visit_constrained)
            end
        end
    end
    new_visits_constrained = sort(new_visits_constrained,  by = x -> x[3])
    if length(new_visits_constrained)>0
        return new_visits_constrained[1][1] , true
    else
        return NewVisit(0,0,0,0,0.0,0.0,0.0, false, ToStoreVisit(SplitCosts(0,0,0,0,0,0),0,"","")), false
    end
end

function replaceFromList(inst::Instance, this_sol::Sol, allparam::AllParameters, boats_to_be_visited)
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
    new_visit, feasible = SelectNewVisitLocalSearch(inst, sol, list_to_be_visited_first)
    while length(list_to_be_visited)>0
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
                new_visit, feasible = SelectNewVisitLocalSearch(inst, sol, list_to_be_visited_first)
            else
                new_visit, feasible = SelectNewVisitLocalSearch(inst, sol, list_to_be_visited)
            end
            if feasible==false
                return initializeSol(inst, allparam), false
            end
        end
    end
    return sol, true
end

function pertubateCost(n,c,sol,inst, listindex)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT = inst
    @unpack M = sol
    @unpack Alpha, Proba, AverageCost, Nbexp, Q = allparam
    listcost = [sol.visits[n][c].store.cost.fuel, sol.visits[n][c].store.cost.delay, sol.visits[n][c].store.cost.waiting, sol.visits[n][c].store.cost.penalty, sol.visits[n][c].store.cost.handling]
    pertubatedcost = 0
    for i in listindex
        pertubatedcost += listcost[i]
    end
    return pertubatedcost
end

function localSearchRemovalReplace(inst::Instance, this_sol::Sol, allparam, paramfixed)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT = inst
    sol=deepcopy(this_sol)
    nb_to_remove = ceil(Int, rand(Uniform(paramfixed.LocalSearchBoatMin,paramfixed.LocalSearchBoatMax))*N)
    start = time_ns()
    boats_to_remove = Vector{Tuple}()
    #listindexcost=sample([1,2,3,4,5,5], 3; replace=false)
    for n in 1:N
        all_penalty_this_boat = Vector{}()
        for (c,p) in enumerate(Pi[n])
            push!(all_penalty_this_boat, sol.visits[n][c].store.cost.all - sol.visits[n][c].store.cost.fuel) #pertubateCost(n,c,sol, inst,listindexcost)) 
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
    random_removals = sample(random_removals, ceil(Int, rand(Uniform(paramfixed.LocalSearchRandomMin,paramfixed.LocalSearchRandomMax))*length(random_removals)); replace=false)

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
    
    elapsed = round((time_ns()-start)/1e9,digits=3)
    print('\n')
    println("Local Search first step: $(elapsed) seconds")
    start = time_ns()
    sol, feasible = replaceFromList(inst, sol, allparam, boats_to_be_visited)
    if feasible
        sol, feasible = replaceFromList(inst, sol, allparam, random_removals)
        elapsed = round((time_ns()-start)/1e9,digits=3)
        print('\n')
        println("Local Search second step: $(elapsed) seconds")
        if feasible
            return sol, true
        else
            return initializeSol(inst, allparam), false
        end
    else
        return initializeSol(inst, allparam), false
    end
end


function manualLocalSearch(inst::Instance, this_sol::Sol, cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost, allparam, paramfixed, timelocal)
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
    print('\n')
    print("Initial cost: ")
    print('\n')
    print(old_cost)

    while elapsed<timelocal
        start_step = time_ns()
        new_sol, feasible = localSearchRemovalReplace(inst, sol, allparam, paramfixed)
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
            new_cost1, delay_cost1, waiting_cost1, penalty_cost1, handling_cost1, fuel_cost1=checkSolutionCost(inst, new_sol)
            new_sol1, new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost = pushTime(inst, new_sol, paramfixed, time_local)
            print('\n')
            print("Local search old cost: ")
            print('\n')
            print(old_cost)
            print('\n')
            print("Local search new cost: ")
            print('\n')
            print(new_cost)
            if new_cost1>new_cost
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
        end
        elapsed = round((time_ns()-start)/1e9,digits=3)
    end
    return sol, old_cost, old_delay_cost, old_waiting_cost, old_penalty_cost, old_handling_cost, old_fuel_cost
end

#################################################################################################################################

function CPLEXoptimizeLocalSearch(inst::Instance, sol::Sol, not_taken)
    m = Model(CPLEX.Optimizer)
    set_optimizer_attribute(m, "CPX_PARAM_EPINT", 1e-8)
    set_optimizer_attribute(m, "CPXPARAM_Threads", 1)
    set_optimizer_attribute(m, "CPX_PARAM_TILIM", 40)
    #set_silent(m)

    @unpack N, Ntot, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT, Nl, gamma, Hc, Dc, Fc, Ic, Pc, beta, ports = inst

    shipsOutFlatt=collect(Iterators.flatten(shipsOut))
    txt_error=""
    feasible=true
    ## Now let's compare with the other boats
    for n in N+1:Ntot-1
        l = ceil(Int, Nl[n]/qli)
        for n_ in n+1:Ntot
            l_ = ceil(Int, Nl[n_]/qli)
            for (c,p) in enumerate(inst.Pi[n])
                t1 = shipsOutFlatt[n-N].time
                t2 = t1 +  shipsOutFlatt[n-N].hand
                x1 = round(Int, shipsOutFlatt[n-N].berth/qli)
                x2 = x1 + l
                for (c_,p_) in enumerate(inst.Pi[n_])
                    if p==p_
                        t1_ = shipsOutFlatt[n_-N].time
                        t2_ = t1_ +  shipsOutFlatt[n_-N].hand
                        x1_ = round(Int, shipsOutFlatt[n_-N].berth/qli)
                        x2_ = x1_ + l_
                    end
                end
            end
        end
    end
    if feasible
        Pi_extend =  deepcopy(inst.Pi)
        for n in 1:N
            append!(Pi_extend[n],[P+1])
        end
        dist= hcat(dist,[0 for p in 1:P])
        dist = [dist;transpose([0 for p in 1:(P+1)])]

        @variable(m,x[n=1:Ntot, 1:length(inst.Pi[n])] >= 1, Int)
        @variable(m,y[n=1:Ntot, 1:length(inst.Pi[n])] >= 0, Int)

        ## The ones to check the positions
        @variable(m,sig[n1=1:Ntot,n2=1:Ntot, 1:length(inst.Pi[n1]), 1:length(inst.Pi[n2])], Bin)
        @variable(m,del[n1=1:Ntot,n2=1:Ntot, 1:length(inst.Pi[n1]), 1:length(inst.Pi[n2])], Bin)
        @variable(m,hand[n=1:Ntot, 1:length(inst.Pi[n])] >= 0, Int)

        
        for n in N+1:Ntot
            for c in 1:length(inst.Pi[n])
                @constraint(m,x[n,c]==round(Int, shipsOutFlatt[n-N].berth/qli)+1)
                @constraint(m,y[n,c]==shipsOutFlatt[n-N].time)
                @constraint(m,hand[n,c]==shipsOutFlatt[n-N].hand)
            end
        end

        
        ## Specific optimization
        @variable(m,v[n=1:N, 1:length(inst.Pi[n]), 1:length(delta)], Bin)
        @variable(m,d[n=1:N, 1:length(inst.Pi[n])] >= 0)
        @variable(m,u[n=1:N, 1:length(inst.Pi[n])] >= 0)

        ## Arriving times
        @variable(m,a[n=1:N, 1:length(inst.Pi[n])] >= 0)
        @variable(m,r[n=1:N, 1:length(inst.Pi[n])] >= 0)
        

        #### Fixing variables to compare to heur
        for n in 1:N
            for c in 1:length(inst.Pi[n])
                if [n,c] ∉ not_taken
                    fix(x[n,c],sol.visits[n][c].b; force=true)
                    fix(y[n,c],sol.visits[n][c].t; force=true)
                end
            end
        end

        ## Berth constraints
        @constraint(m,c1[n=1:N, c=1:length(inst.Pi[n])], x[n,c]+ceil(Int, Nl[n]/qli) <= inst.Bp[inst.Pi[n][c]])

        ##inst.Bp[inst.Pi[n1][c1]]
        ##maxT
        ## Boat constraints
        for n1 in 1:Ntot
            for c1 in 1:length(inst.Pi[n1])
                for n2 in 1:Ntot
                    for c2 in 1:length(inst.Pi[n2])
                        if (n1!=n2 && inst.Pi[n1][c1]==inst.Pi[n2][c2] && ((n1<=N) || (n2<=N)))
                            @constraint(m, x[n1,c1]+ceil(Int, Nl[n1]/qli) <= x[n2,c2]+10000*(1-sig[n1,n2,c1,c2]))
                            @constraint(m, y[n1,c1]+hand[n1,c1] <= y[n2,c2]+10000*(1-del[n1,n2,c1,c2]))
                        end
                    end
                end
            end
        end  
        
        ## Boat constraints dependent variables
        for n1 in 1:Ntot
            for c1 in 1:length(inst.Pi[n1])
                for n2 in n1:Ntot
                    for c2 in 1:length(inst.Pi[n2])
                        if (n1!=n2 && inst.Pi[n1][c1]==inst.Pi[n2][c2] && ((n1<=N) || (n2<=N)))
                            @constraint(m, del[n1,n2,c1,c2]+del[n2,n1,c2,c1]+sig[n1,n2,c1,c2]+sig[n2,n1,c2,c1] >=1)
                        end
                    end
                end
            end
        end  
        

        ## Time constraints
        speeds = length(delta)
        for n in 1:N
            if length(inst.Pi[n])>1
                l=length(inst.Pi[n])-1
                for c in 1:l
                    @constraint(m, y[n,c]+hand[n,c]+sum(v[n,c,s]*delta[s] for s in 1:speeds)*dist[inst.Pi[n][c],inst.Pi[n][c+1]] == a[n,c+1])
                    @constraint(m,sum(v[n,c,s] for s in 1:speeds)==1)
                end
            end
        end

        @constraint(m,c6[n=1:N, c=1:length(inst.Pi[n])], a[n,c] <= y[n,c])
        @constraint(m,c7[n=1:N, c=1:length(inst.Pi[n])],  T[n,c,1]<= y[n,c])
        @constraint(m,c8[n=1:N, c=1:length(inst.Pi[n])],  y[n,c]+hand[n,c]-shipsIn[n].eT[c]<= d[n,c])
        @constraint(m,c9[n=1:N, c=1:length(inst.Pi[n])],  y[n,c]+hand[n,c]-T[n,c,2]<= u[n,c])
        @constraint(m,c10[n=1:N, c=1:length(inst.Pi[n])],   hand[n,c] == (1 + beta*r[n,c]*qli)*ports[inst.Pi[n][c]].minH[shipsIn[n].type])
        @constraint(m,c11[n=1:N, c=1:length(inst.Pi[n])],   x[n,c]-(inst.shipsIn[n].Bi[c]/qli)-1 <= r[n,c])
        @constraint(m,c12[n=1:N, c=1:length(inst.Pi[n])],   (inst.shipsIn[n].Bi[c]/qli)+1-x[n,c] <= r[n,c])

        speeds = length(delta)
        @objective(m, Min, sum(sum(u[n,c]*Pc+d[n,c]*Dc+hand[n,c]*Hc+(y[n,c]-a[n,c])*Ic for c in 1:length(inst.Pi[n]))
        +sum(sum(v[n,c,s]*gamma[shipsIn[n].type,s] for s in 1:speeds)*dist[Pi_extend[n][c],Pi_extend[n][c+1]]*Fc for c in 1:length(inst.Pi[n])) 
        for n in 1:N)) 
        optimize!(m)                  #Solving the model
        
        d = Dict{Any,Any}(
                k => value.(v) for 
                (k, v) in object_dictionary(m) if v isa AbstractArray{VariableRef})
        d["inst"]= inst.Pi
        Nlen=Vector{Int64}()

        for n in Nl
            push!(Nlen,ceil(Int, n/qli))
        end
        d["length_boats"]= Nlen

        T_vec=Vector{Vector{Vector{Int64}}}()
        for n in 1:N
            push!(T_vec,Vector{Vector{Int64}}(undef, 0))
            for c in 1:length(inst.Pi[n])
                push!(T_vec[n],T[n,c,:])
            end
        end

        d["calls"]=T_vec

        d["objectif"]=JuMP.objective_value(m)

        return txt_error, d, JuMP.objective_value(m)
    else
        return txt_error, Dict{Any,Any}(), txt_error
    end
end

function local_search(inst::Instance, sol::Sol, cost::Int64, paramchosen::ChosenParameters, paramfixed::FixedParameters, max_time)
    @unpack N, Ntot, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT, Nl, gamma, Hc, Dc, Fc, Ic, Pc, beta, ports = inst
    new_sol=deepcopy(sol)
    new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost=checkSolutionCost(inst, new_sol)
    list_initialize_boats = []
    list_initialize_all = []
    for n in 1:N
        cost_penalty=0
        for c in 1:length(inst.Pi[n])
            cost_penalty += sol.visits[n][c].store.cost.penalty
            append!(list_initialize_all,[(sol.visits[n][c].store.cost.penalty,[n,c])])
        end
        append!(list_initialize_boats,[(cost_penalty,n)])
    end
    list_initialize_boats = sort(list_initialize_boats, by= x-> x[1], rev=true)
    list_initialize_all = sort(list_initialize_all, by= x-> x[1], rev=true)
    #list_initialize_boats = shuffle(list_initialize_boats)
    #list_initialize_all = shuffle(list_initialize_all)
    start = time_ns()
    elapsed = round((time_ns()-start)/1e9,digits=3)
    tactic = paramchosen.TacticLocalSearch
    first=true
    while elapsed<max_time
        if tactic=="boat"
            new_boat_removed = list_initialize_boats
            new_boat=list_initialize_boats[1][2]
            other_visits=get_random_same_port(inst,new_sol,new_boat, paramfixed.CPLEXBoat)
        else
            if first==true
                other_visits_removed=list_initialize_all[1:min(paramfixed.CPLEXRandom, length(list_initialize_all))]
                other_visits=list_initialize_all[1:min(paramfixed.CPLEXRandom, length(list_initialize_all))]
                for i in 1:length(other_visits)
                    other_visits[i]=other_visits[i][2]
                end
            else
                other_visits_removed=sample(list_initialize_all, min(paramfixed.CPLEXRandom, length(list_initialize_all)))
                other_visits=sample(list_initialize_all, min(paramfixed.CPLEXRandom, length(list_initialize_all)))
                for i in 1:length(other_visits)
                    other_visits[i]=other_visits[i][2]
                end
            end
        end
        new_solCPLEX=deepcopy(new_sol)
        txt_error, d, cost2 = CPLEXoptimizeLocalSearch(inst, new_sol, other_visits)
        for visit in other_visits
            #print('\n')
            #print(visit[1])
            #print('\n')
            #print(visit[2])
            #print('\n')
            #print(d[:x])
            #print('\n')
            #print(d[:y])
            #print('\n')
            #print("############")
            #print('\n')
            #print(d[:x][visit[1],visit[2]])
            #print('\n')
            #print(d[:y][visit[1],visit[2]])
            #print('\n')
            new_solCPLEX.visits[visit[1]][visit[2]].t=round(d[:y][visit[1],visit[2]])           
            new_solCPLEX.visits[visit[1]][visit[2]].b=round(d[:x][visit[1],visit[2]])
        end
        if checkSolutionFeasability(inst, new_solCPLEX)
            new_sol=deepcopy(new_solCPLEX)
            new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost=checkSolutionCost(inst, new_sol)
            #if new_cost<cost
            #    return new_sol, new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost
            #end
        end
        elapsed = round((time_ns()-start)/1e9,digits=3)
    end
    
    return new_sol, new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost
end

function get_random_same_port(inst::Instance, sol::Sol, randomBoat, alpha)
    @unpack N, Ntot, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT, Nl, gamma, Hc, Dc, Fc, Ic, Pc, beta, ports = inst
    visits_to_remove=[]
    for c_boat in 1:length(inst.Pi[randomBoat])
        all_visits_distance=[]
        for n in 1:N
            for c in 1:length(inst.Pi[n])
                if inst.Pi[n][c] == inst.Pi[randomBoat][c_boat] && n!=randomBoat
                    append!(all_visits_distance,[(n,c)])
                end
            end
        end
        all_visits_distance=sample(all_visits_distance, min(alpha,length(all_visits_distance)); replace=false)
        append!(visits_to_remove, all_visits_distance)
    end
    return visits_to_remove
end


function get_closest(inst::Instance, sol::Sol, randomBoat, alpha)
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

#seed=1
#N=4
#Nout=3
#qli=10
#time_heur_max=20
#time_local_search=6
#inst = readInstFromFile("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_$seed"*"_$N"*"_$Nout"*"_$qli"*".txt")
#sol, cost, allparam = GRASP_reactive(instance, "cost","cost", AdjustProba(2,3,3,3,3), 3,10)
#new_sol, cost2 = local_search(inst, sol,4,7)















