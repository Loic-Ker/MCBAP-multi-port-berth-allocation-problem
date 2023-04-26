using UnPack 
using Statistics
import XLSX
include("MBAP_INST.jl")
include("../MBAP_SOL.jl")
include("../check_solution.jl")


# check if two periods overlap
function doOverlap(t, tF, t_, tF_)
    if t_ < t < tF_
        return true
    elseif t_ < tF < tF_
        return true
    elseif t < t_ < tF
        return true
    elseif t < tF_ < tF
        return true
    elseif t == t_ || tF == tF_
        return true
    end
    return false
end

# check if two rectangles overlap : for the final check function
function doOverlapRectangle(x1l, x1r, y1d, y1u, x2l, x2r, y2d, y2u)
    # check if one element is a line
    if x1l == x1r || y1u == y1d || x2l == x2r || y2u == y2d
        return false
    end
    # If one rectangle is on left side of other
    if x1l >= x2r || x2l >= x1r
        return false
    end
    # If one rectangle is above other
    if y1d >= y2u || y2d >= y1u
        return false
    end
    return true
end

# compute the lowest feasible speed to satisfy travel time t
function findLowestSpeed(t::Float64, delta::Vector{Float64}, dist::Int64)
    S = length(delta)
    for s in 1:S
        tt = delta[s]*dist
        if tt < t
            return s
        end
    end
    return -1
end



# update matrix M of occupied positions
function updateM(inst::Instance, n::Int64, c::Int64, b::Int64, t::Int64, nq::Int64, hand::Int64, M)
    @unpack N, Pi, h, shipsIn, Bp, qli = inst
    ## We look at all the other boats and all their visits
    for n_ in 1:N, (c_,p_) in enumerate(Pi[n_])
        ## If it's the same port, and we have a different boat or a different visit for the same boat
        if ((n_ == n && c != c_) || (n_ != n)) && p_ == Pi[n][c]
            l = ceil(Int, shipsIn[n_].l/qli)
            ## We can not not position our new boat between b_boat-length_newboat
            bMin = b - l + 1 < 0.5 ? 1 : b - l + 1
            ## And b_boat+length_boat if the have the same time line
            bMax = b + nq - 1 > Bp[p_]-l+1 ?  Bp[p_]-l+1 : b + nq - 1
            for b_ in bMin:bMax
                ## Then here we check the time line
                ## For each of these position we need the time to not be between :
                handling = ceil(Int, h[n_][c_][b_])
                ## t-handling_newboat
                tMin = t - handling + 1 < -0.1 ? 0 : t - handling + 1
                ## and t+handling_boat
                tMax = t + hand - 1
                ## Then we set everything to false
                M[n_][c_][b_,tMin+1:tMax+1] .= false
            end
        end
    end
    return M
end



# generate matrix M of occupied positions
function generateOccupiedMx(inst::Instance, sol) #, unassigned::Vector{Tuple{Int64, Int64}})
    @unpack N, P, Pi, T, Bp, h, dist, qli, shipsIn, shipsOut, delta, maxT, maxC = inst
    ## The matrix we create contains all the time/distance discretization
    M = [[zeros(Bool, Bp[Pi[n][c]]-ceil(Int, shipsIn[n].l/qli)+1, maxT+1) for c in 1:length(Pi[n])] for n in 1:N]
    ## The between our bounded time intervalls we set the positions as possible (true)
    for i in 1:N
        nq = ceil(Int,shipsIn[i].l/qli)
        for (c,p) in enumerate(Pi[i])
            t1 = max(shipsIn[i].sT[c], T[i,c,1])
            t2 = maxT 
            for b in 1:Bp[p]-nq+1
                hand = ceil(Int, h[i][c][b])
                M[i][c][b,t1+1:t2-hand+2] .= true
            end
        end
    end

    ## Here we remove the positions if the places are already taken by other boats
    ## We use the same tactic that the updateMatrix function 
    for port in 1:P
        for el in shipsOut[port]
            @unpack time, berth, hand, length = el
            b = round(Int, berth/qli) + 1
            nq = ceil(Int, length/qli)
            for n in 1:N, (c,p) in enumerate(Pi[n])
                ## If it's in the same port
                if port == p
                    l = ceil(Int, shipsIn[n].l/qli)
                    ## We can not not position our new boat between b_boat-length_newboat
                    bMin = b - l + 1 < 0.5 ? 1 : b - l + 1
                    ## And b_boat+length_boat if the have the same time line
                    bMax = b + nq - 1 > Bp[p]-l+1 ?  Bp[p]-l+1 : b + nq - 1
                    for b_ in bMin:bMax
                        ## Then here we check the time line
                        ## For each of these position we need the time to not be between :
                        handling = ceil(Int, h[n][c][b_])
                        ## t-handling_newboat
                        tMin = time - handling + 1 < -0.1 ? 0 : time - handling + 1
                        ## and t+handling_boat
                        tMax = time + hand - 1
                        ## Then we set everything to false
                        M[n][c][b_,tMin+1:tMax+1] .= false
                    end
                end
            end
        end
    end

    ## Here we can also create the matrix from a previous constructed solution
    ## This could be great to use for memory algorithm ?
    for n in 1:N
        for (c,val) in enumerate(sol[n])
            if val.planned # !((n,c) in unassigned)
                @unpack p,b,t = val
                # (p,b,t) = val
                nq = ceil(Int,shipsIn[n].l/qli)
                hand = ceil(Int, h[n][c][b])
                M = updateM(inst, n, c, b, t, nq, hand, M)
            end
        end
    end
    return M
end



## This functions permits to find constrained possible solutions for a given port visit of boat n
function findNodes(inst::Instance, sol::Sol, port::Int, t1::Int, t2::Int, n::Int, c::Int, l::Int)
    ## We are gonna store all the positions/time coordinates
    pos = Tuple{Int, Int}[]

    ## Let's get all the visits planned and see if we can plan our new one close to an exisiting one
    for n_ in 1:inst.N
        for (c_,p) in enumerate(inst.Pi[n_])
            if p == port
                if sol.visits[n_][c_].planned
                    @unpack b,t = sol.visits[n_][c_]
                    hand = ceil(Int, inst.h[n_][c_][b])
                    if doOverlap(t1, t2, t, t + hand)
                        len = ceil(Int, inst.shipsIn[n_].l/inst.qli)

                        # Let's fix the times and move in position
                        for x in b-l:b+len
                            if 0.5 < x && x + l <= inst.Bp[p]
                                h = ceil(Int, inst.h[n][c][x])
                                for y in [t - h, t + hand]
                                    if t1 <= y && y + h <= t2
                                        push!(pos, (x, y))
                                    end
                                end
                            end
                        end

                        # Let's fix the positions and move in time
                        for x in [b - l, b + len]
                            if 0.5 < x && x + l <= inst.Bp[p]
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

    # Let's do the same but this we considered external port visits fixed since the beginning
    # We get also all the positions feasible here

    for extS in inst.shipsOut[port]
        @unpack time, berth, hand, length = extS
        len = length
        t = time
        b = berth
        if doOverlap(t1, t2, time, time + hand)
            # Let's fix the times and move in position
            for x in b-l:b+len
                if 0.5 < x && x + l <= inst.Bp[port]
                    h = ceil(Int, inst.h[n][c][x])
                    for y in [t - h, t + hand]
                        if t1 <= y && y + h <= t2
                            push!(pos, (x, y))
                        end
                    end
                end
            end

            # Let's fix the positions and move in time
            for x in [b - l, b + len]
                if 0.5 < x && x + l <= inst.Bp[port]
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

    # We add here the berth limits
    for t in t1:t2
        h0 = ceil(Int, inst.h[n][c][1])
        hB = ceil(Int, inst.h[n][c][inst.Bp[port]-l])
        if t + h0 <= t2
            push!(pos, (1, t))
            push!(pos, (inst.Bp[port]-l, t))
        end
    end

    # And then here the solution when we start at t1
    for x in 1:inst.Bp[port]-l
        h = ceil(Int, inst.h[n][c][x])
        if t1 + h <= t2
            push!(pos, (x,t1))
        end
    end
    return pos
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
        end
        ## Change here, we do not want to take into account in the optimization the 
        ## differences in distances between the different ports
        cost += Fc*(sum(gamma[:,s])/3)*(sum(dist)/2/P)/2
    end
    if c > 1 && sol.visits[n][c-1].planned
        pP = sol.visits[n][c-1].p
        bP = sol.visits[n][c-1].b
        tP = sol.visits[n][c-1].t
        hP = h[n][c-1][bP]
        deltaT = t - (tP + hP)
        s = findLowestSpeed(deltaT, delta, dist[pP,p])
        ## Here we can remove this part too if we do not want to take the fuel into account
        cost += Fc*(sum(gamma[:,s])/3)*(sum(dist)/2/P)/2
        wait = t - (tP + hP + delta[s]*dist[pP,p])
        cost += Ic*wait
    end
    return cost
end





function checkPreviousNextVisit(inst::Instance, sch::Vector{Visit}, n::Int, c::Int, port::Int)
    @unpack N, h, qli, shipsIn, Pi, Bp, delta, dist, T, maxC, maxT = inst
    t1 = max(T[n,c,1], shipsIn[n].sT[c])
    t2 = T[n,c,2]
    # if previous visit is planned, update earliest tA
    if c > 1 && sch[c-1].planned #(n,c-1) in placed
        pP = sch[c-1].p
        bP = sch[c-1].b
        tP = sch[c-1].t
        # (pP,bP,tP) = sol.visits[n][c-1]
        hP = h[n][c-1][bP]
        tA = ceil(Int, tP + hP + delta[end]*dist[pP,port])
        if tA>t1 
            t1=tA
        end
    end
    # if next visit is planned, update latest tD
    if c < length(Pi[n]) && sch[c+1].planned #(n,c+1) in placed
        pN = sch[c+1].p
        bN = sch[c+1].b
        tN = sch[c+1].t
        # (pN,bN,tN) = sol.solution[n][c+1]
        tD = ceil(Int, tN - delta[end]*dist[port,pN])
        if tD<0
            tD=0
        end
        if t2 > tD 
            t2 = tD
        end
    end
    return t1, t2
end

# compute feasible positions (b,t) in two different ways :  
## Only looking at the cost and returning the best
## Get only the most constrained and return only the best most constrained
function findBestFeasiblePositions(inst::Instance, sol::Sol, correction::Bool) #, placed::Vector{Tuple{Int64,Int64}})
    #TODO: t1 cannot be greater than t2!
    # adapt to given planned previous visits for each ship
    @unpack N, h, qli, shipsIn, Pi, Bp, delta, dist, T, maxC, maxT = inst
    @unpack M = sol
    # C = maximum(length.(Pi))
    best_per_visit = Tuple{Int64,Int64,Int64,Int64,Float64}[]
    best_constrained_per_visit = Tuple{Int64,Int64,Int64,Int64,Float64}[]
    for n in 1:N
        l = ceil(Int, shipsIn[n].l/qli)
        sch = sol.visits[n]
        for (c,p) in enumerate(Pi[n])
            this_visit_bt = Tuple{Int64,Int64,Int64,Int64,Float64}[]
            this_visit_constrained_bt = Tuple{Int64,Int64,Int64,Int64,Float64}[]
            find_feasible = false
            if sch[c].planned == false 
                ## We change t1 if you have a visit before 
                ## And t2 if we planned a visit after
                t1, t2 = checkPreviousNextVisit(inst, sch, n, c, p)
                ## Check all the possible solutions
                for b in 1:Bp[p]-l
                    hand = ceil(Int, h[n][c][b])
                    for t in t1:t2
                        if t + hand <= t2
                            if M[n][c][b,t+1]
                                find_feasible=true
                                this_cost = computeCostPosOptim(inst, n, c, b, t,sol)
                                push!(this_visit_bt, (n,c,b,t,this_cost))
                            end
                        end
                    end
                end

                ## Check the constrained ones
                pos = findNodes(inst, sol, p, t1, t2, n, c, l)
                for (b,t) in pos
                    hand = ceil(Int, h[n][c][b])
                    if t + hand <= t2
                        if sol.M[n][c][b,t+1]
                            find_feasible=true
                            this_cost = computeCostPosOptim(inst, n, c, b, t ,sol)
                            push!(this_visit_constrained_bt, (n,c,b,t,this_cost))
                        end
                    end
                end


                if find_feasible==false && correction
                    for b in 1:inst.Bp[p]-l
                        hand = ceil(Int, inst.h[n][c][b])
                        t1 = max(T[n,c,2]-hand+1, t1)
                        t2 = inst.maxT
                        c_ = length(inst.Pi[n])
                        while c_ != c
                            t2 -= ceil(Int, minimum(inst.h[n][c_])) + ceil(Int, inst.dist[inst.Pi[n][c_-1], inst.Pi[n][c_]]*inst.delta[end])
                            c_ -= 1
                        end
                    end
                    for b in 1:Bp[p]-l #b in 1:Bp[p]-l+1
                        hand = ceil(Int, h[n][c][b])
                        for t in t1:t2-1 # this can be tightened
                            if t + hand <= t2
                                #TODO: replace format of M, instead M[p][n][b,t] = true, ship n can ship at (p,b,t)
                                if M[n][c][b,t+1]
                                # if !(false in M[p][b:b+l-1, t+1:t+hand])
                                    find_feasible=true
                                    this_cost = computeCostPosOptim(inst, n, c, b, t+1 ,sol)
                                    push!(this_visit_bt, (n,c,b,t+1,this_cost))
                                end
                            end
                        end
                    end
                end
                if length(this_visit_bt)>0.1
                    best_this_visit = sort(this_visit_bt, by= x-> x[end])[1]
                    push!(best_per_visit, best_this_visit)
                end
                if length(this_visit_constrained_bt)>0.1
                    best_this_visit_constrained = sort(this_visit_constrained_bt, by= x-> x[end])[1]
                    push!(best_constrained_per_visit, best_this_visit_constrained)
                end
            end
        end
    end

    best_per_visit = sort(best_per_visit, by= x-> x[end])
    #print('\n')
    #print("Set of best visits found")
    #print('\n')
    #print(best_per_visit)
    #print('\n')
    best_constrained_per_visit = sort(best_constrained_per_visit, by= x-> x[end])
    #print(best_constrained_per_visit)
    #print('\n')
    #print('\n')
    return best_per_visit, best_constrained_per_visit
end

function initialize(inst::Instance)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp = inst
    sol = Sol(inst)
    # solution = [Vector{Tuple{Int64,Int64,Int64}}() for n in 1:N] # (p,b,t)
    for n in 1:N, (c,p) in enumerate(Pi[n])
        sol.visits[n][c].minT = max(shipsIn[n].sT[c], T[n,c,1])
        sol.visits[n][c].maxT = T[n,c,2]
    end
    sol.M = generateOccupiedMx(inst, sol.visits)
    return sol
end

function greedyrandomizedconstruction(inst::Instance, correction, max_time)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp = inst
    sol = initialize(inst)
    not_find_total = 0
    ordered_visits, constrained_ordered_visits = findBestFeasiblePositions(inst, sol, correction)
    continue_=true
    start = time_ns()
    elapsed = round((time_ns()-start)/1e9,digits=3)
    while continue_ && elapsed<max_time
        if length(ordered_visits)<1
            ordered_visits=[]
        end
        if length(constrained_ordered_visits)<1
            constrained_ordered_visits=[]
        end
        visits_to_select = vcat(ordered_visits[1:min(length(ordered_visits),3)],constrained_ordered_visits[1:min(length(constrained_ordered_visits),3)])
        n,c,b,t,this_cost = visits_to_select[rand(1:min(length(visits_to_select),6))]
        l = ceil(Int, shipsIn[n].l/qli)
        hand = ceil(Int, h[n][c][b])
        M = sol.M
        sol.M = updateM(inst, n, c, b, t, l, hand, M)
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

        ordered_visits, constrained_ordered_visits = findBestFeasiblePositions(inst, sol, correction)
        if continue_
            if length(ordered_visits)<1 && length(constrained_ordered_visits)<1
                not_find_total = not_find_total+1
                print("###################")
                print('\n')
                print("We do not find a solution")
                print('\n')
                print(sol.visits)
                print('\n')

                sol = initialize(inst)
                ordered_visits, constrained_ordered_visits = findBestFeasiblePositions(inst, sol, correction)
            end
        end
        elapsed = round((time_ns()-start)/1e9,digits=3)
    end

    #print("###################")
    #print('\n')
    #print("Number of time we had to stop the heuristic")
    #print('\n')
    #print(not_find_total)
    #print('\n')
    return sol #, solution
end



#instance = readInstFromFile("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_1_4_3_10.txt")
#sol = greedyrandomizedconstruction(instance, false, 5)

function GRASP_V1(inst::Instance, correction, max_time_heur, max_time)
    cost=1000000
    sol = initialize(inst)
    start = time_ns()
    elapsed = round((time_ns()-start)/1e9,digits=3)
    while elapsed<max_time
        new_sol = greedyrandomizedconstruction(inst, correction, max_time_heur)
        feasible = true
        for boat_vis in new_sol.visits
            for vis in boat_vis
                if vis.planned == false
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

function printSol(instance)
    sol, cost = GRASP_V1(instance, true, 15, 50)

    print('\n')
    print("The solution :")
    print('\n')
    print(sol.visits)
    print('\n')
    print("And the cost is ")
    print('\n')
    print(cost)
    return(cost)
    return cost
end

newbenchmark = [Any["Seed","N","Nout","qli","Time","NewCost"]]
print("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_1_14_3_40.txt")
inst = readInstFromFile("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_1_14_3_40.txt")
cost = printSol(inst)
this_benchmark=[2,4,3,40,60, ceil(Int,cost)]
push!(newbenchmark,this_benchmark)
newbenchmark

#for N in 14:14
#    for Nout in 3:5
#        for qli in [10,20,40,80]
#            for seed in 1:5
#                print("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_$seed"*"_$N"*"_$Nout"*"_$qli"*".txt")
#                inst = readInstFromFile("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_$seed"*"_$N"*"_$Nout"*"_$qli"*".txt")
#                cost = printSol(inst)
#                this_benchmark=[seed,N,Nout,qli,60, ceil(Int,cost)]
#                push!(newbenchmark,this_benchmark)
#            end
#        end
#    end
#end

#print(mapreduce(permutedims, vcat, newbenchmark))
#XLSX.writetable("N14.xlsx",
#    sheet1=(collect(eachcol(mapreduce(permutedims, vcat, newbenchmark))), string.("col", 1:6)) )


