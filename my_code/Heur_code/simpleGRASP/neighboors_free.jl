using UnPack 
using Statistics
import XLSX
using CSV, Tables
using DataFrames
using Random
include("../../MBAP_INST.jl")
include("../MBAP_SOL.jl")
include("../toolsMatrixTimes.jl")
include("../check_solution.jl")



function sortSchedule(best_per_visit_by_schedule)
    best_per_visit_by_schedule1 = filter(t -> t[2]==1, best_per_visit_by_schedule)
    if length(best_per_visit_by_schedule1)>0
        return best_per_visit_by_schedule1
    end
    
    best_per_visit_by_schedule2 = filter(t -> t[2]==2, best_per_visit_by_schedule)
    if length(best_per_visit_by_schedule2)>0
        return best_per_visit_by_schedule2
    end
    
    best_per_visit_by_schedule3 = filter(t -> t[2]==3, best_per_visit_by_schedule)
    if length(best_per_visit_by_schedule3)>0
        return best_per_visit_by_schedule3
    end

    
    return best_per_visit_by_schedule 
end


function initialize(inst::Instance)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT = inst
    sol = Sol(inst)
    # solution = [Vector{Tuple{Int64,Int64,Int64}}() for n in 1:N] # (p,b,t)
    for n in 1:N, (c,p) in enumerate(Pi[n])
        sol.visits[n][c].minT = max(shipsIn[n].sT[c], T[n,c,1])
        sol.visits[n][c].maxT = min(5*T[n,c,2],maxT)
    end
    sol.M = generateOccupiedMx(inst, sol.visits)
    sol = updateTimesInitialization(inst,sol)
    return sol
end


function findConstrainedPos(inst::Instance, sol::Sol, port::Int, t1::Int, t2::Int, n::Int, c::Int, l::Int)
    # find ship visits during t1-t2
    pos = Tuple{Int, Int, Float64}[]
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
                            if 0.5 < x && x + l <= inst.Bp[p]
                                h = ceil(Int, inst.h[n][c][x])
                                for y in [t - h, t + hand]
                                    if t1 <= y && y + h <= t2
                                        distance = abs(x-b)/(b+1)
                                        push!(pos, (x, y, distance))
                                    end
                                end
                            end
                        end
                        for x in [b - l, b - l-1, b + len, b+len+1]
                            if 0.5 < x && x + l <= inst.Bp[p]
                                h = ceil(Int, inst.h[n][c][x])
                                for y in t-h:t+hand
                                    if t1 <= y && y + h <= t2
                                        distance = abs(x-b)/(b+1)
                                        push!(pos, (x, y, distance))
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
                if 0.5 < x && x + l <= inst.Bp[port]
                    h = ceil(Int, inst.h[n][c][x])
                    if ix in [2,3]
                        if t1 <= time - h && time <= t2
                            distance = (abs(x-berth)/(berth+1))
                            push!(pos, (x, time - h, distance))
                        end
                        if t1 <= time + hand && time + hand + h <= t2
                            distance = (abs(x-berth)/(berth+1))
                            push!(pos, (x, time + hand, distance))
                        end
                    end
                    if ix in [1,4]
                        if t1 <= time && time + h <= t2
                            distance = (abs(x-berth)/(berth+1))
                            push!(pos, (x, time, distance))
                        end
                        if t1 <= time + hand - h && time + hand <= t2
                            distance = (abs(x-berth)/(berth+1))
                            push!(pos, (x, time + hand - h, distance))
                        end
                    end
                end
            end
        end
    end
    for t in t1:t2
        h0 = ceil(Int, inst.h[n][c][1])
        hB = ceil(Int, inst.h[n][c][inst.Bp[port]-l-1])
        if t + h0 <= t2
            push!(pos, (1, t, 0))
            push!(pos, (inst.Bp[port]-l-1, t, 0))
        end
    end
    return pos
end

# compute feasible positions (b,t) in different ways :  
## Order by best costs
## Get random ones
## Get the ones with less feasible solutions
## Get the ones closer to other ports visit
## Ordered by port visit number et the end
function findOneFeasiblePositions(inst::Instance, sol::Sol, n, c, p) 
    @unpack N, h, qli, shipsIn, Pi, Bp, delta, dist, T, maxC, maxT = inst
    @unpack M = sol
    l = ceil(Int, shipsIn[n].l/qli)
    sch = sol.visits[n]
    this_visit_bt = Tuple{Int64,Int64,Int64,Int64,Float64,Int64}[]
    best_this_visit_bt = Tuple{Int64,Int64,Int64,Int64,Float64,Int64}[]
    this_visit_constrained_bt = Tuple{Int64,Int64,Int64,Int64,Float64,Float64}[]
    best_this_visit_constrained_bt = Tuple{Int64,Int64,Int64,Int64,Float64,Float64}[]
    find_feasible = false
    t1=sol.visits[n][c].minT
    t2=sol.visits[n][c].maxT
    if sch[c].planned == false
        ## Check all the possible solutions
        for b in 1:Bp[p]-l
            hand = ceil(Int, h[n][c][b])
            for t in t1:t2
                if t + hand <= t2
                    if M[n][c][b,t+1]
                        this_cost, feas = computeCostPosOptim(inst, n, c, b, t,sol)
                        if feas
                            push!(this_visit_bt, (n,c,b,t,this_cost,0))
                        end
                    end
                end
            end
        end  
        
        ## Check the constrained ones
        pos = findConstrainedPos(inst, sol, p, t1, t2, n, c, l)
        for (b,t,distance) in pos
            hand = ceil(Int, h[n][c][b])
            if t + hand <= t2
                if sol.M[n][c][b,t+1]
                    this_cost, feas = computeCostPosOptim(inst, n, c, b, t ,sol)
                    if feas
                        push!(this_visit_constrained_bt, (n,c,b,t,this_cost,distance))
                    end
                end
            end
        end
    end

    #if length(best_this_visit_constrained)>0
    #    return best_this_visit_constrained
    #end
    if length(this_visit_constrained_bt)>0.1
        ### Not free loop
        ## Get best of this visit by minimum distance
        best_this_visit_constrained_bt = sort(this_visit_constrained_bt, by= x-> x[5])[1:min(3,length(this_visit_constrained_bt))]
    end

    if length(this_visit_bt)>0.1
        best_this_visit_bt = sort(this_visit_bt, by= x-> x[5])[1:min(3,length(this_visit_bt))]
    end
    
    return [best_this_visit_bt;best_this_visit_constrained_bt]
end

# compute feasible positions (b,t) in different ways :  
## Order by best costs
## Get random ones
## Get the ones with less feasible solutions
## Get the ones closer to other ports visit
## Ordered by port visit number et the end
function findFeasiblePositions(inst::Instance, sol::Sol) #, placed::Vector{Tuple{Int64,Int64}})
    #TODO: t1 cannot be greater than t2!
    # adapt to given planned previous visits for each ship
    @unpack N, h, qli, shipsIn, Pi, Bp, delta, dist, T, maxC, maxT = inst
    @unpack M = sol
    # C = maximum(length.(Pi))
    random_per_visit = Tuple{Int64,Int64,Int64,Int64,Float64,Int64}[]
    best_per_visit = Tuple{Int64,Int64,Int64,Int64,Float64,Int64}[]
    best_per_visit_constrained = Tuple{Int64,Int64,Int64,Int64,Float64,Float64}[]
    for n in 1:N
        l = ceil(Int, shipsIn[n].l/qli)
        sch = sol.visits[n]
        for (c,p) in enumerate(Pi[n])
            this_visit_bt = Tuple{Int64,Int64,Int64,Int64,Float64,Int64}[]
            this_visit_constrained_bt = Tuple{Int64,Int64,Int64,Int64,Float64,Float64}[]
            find_feasible = false
            t1=sol.visits[n][c].minT
            t2=sol.visits[n][c].maxT
            if sch[c].planned == false
                count=0
                ## Check all the possible solutions
                for b in 1:Bp[p]-l
                    hand = ceil(Int, h[n][c][b])
                    for t in t1:t2
                        if t + hand <= t2
                            if M[n][c][b,t+1]
                                this_cost, feas = computeCostPosOptim(inst, n, c, b, t,sol)
                                if feas
                                    count+=1
                                    push!(this_visit_bt, (n,c,b,t,this_cost,0))
                                end
                            end
                        end
                    end
                end  
                
                ## Check the constrained ones
                pos = findConstrainedPos(inst, sol, p, t1, t2, n, c, l)
                for (b,t,distance) in pos
                    hand = ceil(Int, h[n][c][b])
                    if t + hand <= t2
                        if sol.M[n][c][b,t+1]
                            this_cost, feas = computeCostPosOptim(inst, n, c, b, t ,sol)
                            if feas
                                push!(this_visit_constrained_bt, (n,c,b,t,this_cost,distance))
                            end
                        end
                    end
                end


                if length(this_visit_bt)>0.1
                    this_visit_bt =  [sort(this_visit_bt, by= x-> x[5])[1]]
                    if length(this_visit_bt)>0
                        for visit in this_visit_bt
                            n = visit[1]
                            c = visit[2]
                            b = visit[3]
                            t = visit[4]
                            this_cost = visit[5]
                            push!(best_per_visit, (n,c,b,t,this_cost,count))
                        end
                    end
                end
                
                if length(this_visit_constrained_bt)>0.1
                    ### Not free loop
                    ## Get best of this visit by minimum distance
                    best_this_visit_constrained_bt = [sort(this_visit_constrained_bt, by= x-> x[5])[1]]

                    if length(best_this_visit_constrained_bt)>0
                        for visit in best_this_visit_constrained_bt
                            push!(best_per_visit_constrained, visit)
                        end
                    end
                end
            end
        end
    end

    if length(best_per_visit_constrained)>0
    #    best_per_visit_by_schedule = sortSchedule(best_per_visit_constrained)
    #    best_per_visit_by_constrained = sort(best_per_visit_by_schedule, by= x-> x[end]) #best_per_visit_by_schedule #sort(best_per_visit_by_schedule, by= x-> x[end])
    #    best_per_visit_by_constrained = best_per_visit_constrained
    #    return best_per_visit_by_constrained
        #best_per_visit_constrained = sort(best_per_visit_constrained, by= x-> (x[end],x[5]))[1:min(3,length(best_per_visit_constrained))]
        best_per_visit_constrained = sort(best_per_visit_constrained, by= x-> x[5])[1:min(3,length(best_per_visit_constrained))]
    end

    if length(best_per_visit)>0
        #best_per_visit = sort(best_per_visit, by= x-> x[end])
        #best_per_visit_by_schedule = best_per_visit
        #best_per_visit_by_schedule = sortSchedule(best_per_visit_by_schedule)
        #best_per_visit_by_schedule = sort(best_per_visit_by_schedule, by= x-> (x[end],x[5]))
        #return best_per_visit_by_schedule

        best_per_visit = sort(best_per_visit, by= x-> x[5])[1:min(3,length(best_per_visit))]
    end
    return [best_per_visit_constrained;best_per_visit]
    #return best_per_visit
end

# compute the cost of ship n at port visit c berthing at pos b and time t given solution sol for the optimization process
function computeCostPosOptim(inst::Instance, n::Int64, c::Int64, b::Int64, t::Int64 ,sol::Sol)
    @unpack P, h, Pc, shipsIn, T, delta, gamma, dist, Hc, Dc, Fc, Ic = inst
    @unpack type, eT = shipsIn[n]
    p = inst.Pi[n][c]
    cost = Hc*h[n][c][b]
    delay = t + h[n][c][b] > eT[c] ? t + h[n][c][b] - eT[c] : 0.
    cost += Dc*delay
    penalty = t + h[n][c][b] > T[n,c,2] ? Pc*(t + h[n][c][b] - T[n,c,2]) : 0.
    cost += penalty
    #########
    ## Remove this part for later : do not take into account the fuel
    if c < length(sol.visits[n]) && sol.visits[n][c+1].planned
        pN = sol.visits[n][c+1].p
        bN = sol.visits[n][c+1].b
        tN = sol.visits[n][c+1].t
        # (pN,bN,tN) = sch[c+1]
        deltaT = tN - (t + h[n][c][b])
        s = findLowestSpeed(deltaT, delta, dist[p,pN])
        if s == -1
            # this should only happen if we quit the port too late
            s = length(delta)
            return 1000000, false
        end
        ## Change here, we do not want to take into account in the optimization the 
        ## differences in distances between the different ports
        cost += Fc*gamma[type,s]*dist[p,pN]
    end
    if c > 1 && sol.visits[n][c-1].planned
        pP = sol.visits[n][c-1].p
        bP = sol.visits[n][c-1].b
        tP = sol.visits[n][c-1].t
        hP = h[n][c-1][bP]
        deltaT = t - (tP + hP)
        s = findLowestSpeed(deltaT, delta, dist[pP,p])
        if s == -1
            # this should only happen if we quit the port too late
            s = length(delta)
            return 1000000, false
        end
        ## Here we can remove this part too if we do not want to take the fuel into account
        cost += Fc*gamma[type,s]*dist[p,pP]
        wait = t - (tP + hP + delta[s]*dist[pP,p])
        cost += Ic*wait
    end
    return cost, true
end


function greedyrandomizedconstruction(inst::Instance, max_time)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT = inst
    sol = initialize(inst)
    ordered_visits = findFeasiblePositions(inst, sol)
    continue_=true
    start = time_ns()
    elapsed = round((time_ns()-start)/1e9,digits=3)

    while continue_ && elapsed<max_time
        print('\n')
        print("Set of visits found")
        print('\n')
        print(ordered_visits[1:min(length(ordered_visits),6)])
        print('\n')
        n,c,b,t,this_cost,count = ordered_visits[rand(1:min(length(ordered_visits),5))]
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
        
        #print('\n')
        #print(sol.visits)
        #print('\n')
        #print("Solution at this point")
        #print('\n')
        
        continue_=false
        for boat_vis in sol.visits
            for vis in boat_vis
                if vis.planned == false
                    continue_ = true
                end
            end
        end

        ordered_visits = findFeasiblePositions(inst, sol)
        if continue_
            if length(ordered_visits)<1 
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
                                    ordered_visits = findOneFeasiblePositions(inst, sol, n, i, Pi[n][i])
                                    if length(ordered_visits)>0
                                        n_,c_,b_,t_,this_cost_,count_ = ordered_visits[rand(1:min(length(ordered_visits),6))]
                                        l_ = ceil(Int, shipsIn[n_].l/qli)
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
                                    else
                                        return initialize(inst)
                                    end
                                end
                            else
                                return initialize(inst)
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



function GRASP_V1(inst::Instance, max_time_heur, max_time)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp = inst
    cost=10000000
    sol = initialize(inst)
    start = time_ns()
    elapsed = round((time_ns()-start)/1e9,digits=3)
    while elapsed<max_time
        new_sol = greedyrandomizedconstruction(inst, max_time_heur)
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
            new_cost=checkSolutionCost(inst, new_sol)
            if new_cost<cost
                sol=new_sol
                cost=new_cost
            end
        end
        elapsed = round((time_ns()-start)/1e9,digits=3)
    end
    return sol, cost
end




#instance = readInstFromFile("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_2_4_3_10.txt")
#sol = greedyrandomizedconstruction(instance, false, 5)

function printSol(instance, t1, t2)
    sol, cost = GRASP_V1(instance, t1, t2)

    print('\n')
    print("The solution :")
    print('\n')
    print(sol.visits)
    print('\n')
    print("And the cost is ")
    print('\n')
    print(cost)
    d=prepareSol(instance, sol, cost)
    return cost,d
end

function prepareSol(inst, sol, cost)
    @unpack N, Ntot, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT, Nl = inst
    d = Dict{Any,Any}()
    
    ## Get the instances 
    d["inst"]= inst.Pi

    ## Get the port calls
    T_vec=Vector{Vector{Vector{Int64}}}()
    for n in 1:N
        push!(T_vec,Vector{Vector{Int64}}(undef, 0))
        for c in 1:length(inst.Pi[n])
            push!(T_vec[n],T[n,c,:])
        end
    end
    d["calls"]=T_vec

    ## Get the boats length :
    Nlen=Vector{Int64}()
        for n in Nl
            push!(Nlen,ceil(Int, n/qli))
        end
    d["length_boats"]= Nlen
    
    x = Vector{Vector{Any}}(undef, 0)
    y = Vector{Vector{Any}}(undef, 0)
    hand = Vector{Vector{Any}}(undef, 0)
    v = Vector{Vector{Any}}(undef, 0)
    a = Vector{Vector{Any}}(undef, 0)
    for n in 1:N
        x_thisboat= Vector{Any}(undef, 0)
        y_thisboat= Vector{Any}(undef, 0)
        hand_thisboat= Vector{Any}(undef, 0)
        v_thisboat= Vector{Any}(undef, 0)
        a_thisboat= Vector{Any}(undef, 0)
        for (c,p) in enumerate(inst.Pi[n])
            push!(x_thisboat, sol.visits[n][c].b)
            push!(y_thisboat, sol.visits[n][c].t)
            push!(hand_thisboat, h[n][c][sol.visits[n][c].b])
            if c < length(inst.Pi[n])
                pN = sol.visits[n][c+1].p
                bN = sol.visits[n][c+1].b
                tN = sol.visits[n][c+1].t
                # (pN,bN,tN) = sch[c+1]
                deltaT = tN - (sol.visits[n][c].t +  h[n][c][sol.visits[n][c].b])
                s = findLowestSpeed(deltaT, delta, dist[sol.visits[n][c].p,pN])
                push!(v_thisboat, ceil(Int, s))
            end
            if c>1
                pN = sol.visits[n][c-1].p
                bN = sol.visits[n][c-1].b
                tN = sol.visits[n][c-1].t
                # (pN,bN,tN) = sch[c+1]
                deltaT = sol.visits[n][c].t - (tN +  h[n][c-1][bN])
                s = findLowestSpeed(deltaT, delta, dist[pN,sol.visits[n][c].p])
                push!(a_thisboat, tN + h[n][c-1][bN] + delta[s]*dist[pN,sol.visits[n][c].p])
            end
            if c==1
                push!(a_thisboat,sol.visits[n][c].t)
            end
        end
        push!(x,x_thisboat)
        push!(y,y_thisboat)
        push!(hand,hand_thisboat)
        push!(v,v_thisboat)
        push!(a,a_thisboat)
    end

    shipsOutFlatt=collect(Iterators.flatten(shipsOut))
    for n in N+1:Ntot
        x_thisboat= Vector{Any}(undef, 0)
        y_thisboat= Vector{Any}(undef, 0)
        hand_thisboat= Vector{Any}(undef, 0)
        for (c,p) in enumerate(inst.Pi[n])
            push!(x_thisboat, round(shipsOutFlatt[n-N].berth/qli)+1)
            push!(y_thisboat, shipsOutFlatt[n-N].time)
            push!(hand_thisboat, shipsOutFlatt[n-N].hand)
        end
        push!(x,x_thisboat)
        push!(y,y_thisboat)
        push!(hand,hand_thisboat)
    end          
    d["x"]=x
    d["y"]=y
    d["a"]=a
    d["v"]=v
    d["hand"]=hand
    d["objectif"]=cost
    return d
end

newbenchmark = [Any["Seed","N","Nout","qli","OldLB","OldUB","OldTime","NewCost","MyTime"]]

import XLSX



function makeSolHeur()
    xf = CSV.read("D:/DTU-Courses/DTU-Thesis/berth_allocation/bernardo_bench/Small_Inst_Res.csv", DataFrame)
    newbenchmark = DataFrame(Seed= [0],N= [0],Nout= [0],qli= [0],OldLB= [0],OldUB= [0],OldTime= [0],HeurCost= [0])
    for N in 10:10
        for qli in [80]
            for Nout in 5:5
                for seed in 1:2
                    print("The instance : $seed"*"_$N"*"_$Nout"*"_$qli")
                    inst = readInstFromFile("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_$seed"*"_$N"*"_$Nout"*"_$qli"*".txt")
                    cost, d = printSol(inst, 30, 90)
                    CSV.write("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/simpleGRASP/sols/HEUR_sol_$seed"*"_$N"*"_$Nout"*"_$qli"*".csv", d)
                    filtering=xf[(xf.Seed.==seed) .& (xf.N.==N) .& (xf.qli.==qli) .& (xf.Nout.==Nout),:]
                    LB=filtering.LB[1]
                    UB=filtering.UB[1]
                    OldTime = filtering.Time[1]
                    this_benchmark=DataFrame(Seed= [seed],N= [N],Nout= [Nout],qli= [qli],OldLB= [LB],OldUB= [UB],OldTime= [OldTime],HeurCost= [ ceil(Int, cost)])
                    newbenchmark=append!(newbenchmark,this_benchmark)
                end
            end
        end
    end
    return newbenchmark
end

#newbenchmark = makeSolHeur()

#CSV.write("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/greedy_only/HEUR_N4_N12.csv", newbenchmark)

inst = readInstFromFile("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_1_4_3_10.txt")
sol, cost= GRASP_V1(inst, 10,60)
print('\n')
print(cost)