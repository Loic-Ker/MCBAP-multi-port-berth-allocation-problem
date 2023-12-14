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
include("GRASP.jl")



#################### Here we define the pathrelinking step 
## It's very similar to the greedy construction, but start from an already constructed solution
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