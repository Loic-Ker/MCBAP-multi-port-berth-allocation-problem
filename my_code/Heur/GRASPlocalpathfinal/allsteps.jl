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
include("PathrelinkingStep.jl")

#################### The complete algorithm

## This function puts all the steps together
## Plus I implemented the paht relinking iteration and elite pool construction here
## So first a method to compute the distance between the solutions based on their positions in the 2D space
function DistanceSols(inst::Instance, new_sol, sol)
    @unpack N, P, Pi, shipsIn, shipsOut, h, dist, delta, qli, T, Bp = inst
    tot_distance = 0
    for n in 1:N
        # The times :
        for (c,p) in enumerate(inst.Pi[n])
            tot_distance = tot_distance + ((new_sol.visits[n][c].t-sol.visits[n][c].t)/(sol.visits[n][c].t+1))+((new_sol.visits[n][c].b-sol.visits[n][c].b)/(sol.visits[n][c].b+1))
        end
    end
    return tot_distance
end

## Then the main function
function GRASP_reactive(seed,N,Nout,qli, type1, type2, type3, paramfixed, time_local, max_time_heur, max_time, expname, location)
    inst = readInstFromFile(location*"MCBAP-multi-port-berth-allocation-problem/Large/CP2_Inst_$seed"*"_$N"*"_$Nout"*"_$qli"*".txt")
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp = inst
    bestcost=1000000000
    cost=1000000000
    allparam = initializeParam(paramfixed)
    sol=initializeSol(inst, allparam)
    start = time_ns()
    elapsed = round((time_ns()-start)/1e9,digits=3)
    start_iter=time_ns()
    elapsed_iter = round((time_ns()-start_iter)/1e9,digits=3)
    nb_iter=0
    d_alliter_after=Dict()
    when_list = Vector{Tuple}()
    when_dict=Dict()
    min_cost_heur=1000000000
    greedy_no_improve=0
    usedcplex=0
    nb_iter_restart_params = 0
    first_relinking = true
    elite_pool = Vector{}()
    distance_sols_min = 10000000000
    old_elite_cost_mean=1000000000
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
    min_distance_elite=10000000000

    while elapsed<max_time

        ## First the greedy randomized construction
        nb_iter_restart_params = nb_iter_restart_params + 1
        paramchosen = ChooseAllParam(allparam, paramfixed)
        start_heur = time_ns()
        random_value = rand()

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
        elapsed_heur = round((time_ns()-start_heur)/1e9,digits=3)
        if feasible && checkSolutionFeasability(inst, new_sol)

            new_cost_heur, delay_cost_heur, waiting_cost_heur, penalty_cost_heur, handling_cost_heur, fuel_cost_heur = checkSolutionCost(inst, new_sol)
            
            start_local=time_ns()
            
            ## We check if the GRASP solution is better than the best GRASP one found so far
            if new_cost_heur<min_cost_heur
                min_cost_heur=deepcopy(new_cost_heur)
                best_sol_heur=deepcopy(new_sol)
            end
            if new_cost_heur>=min_cost_heur 
                greedy_no_improve=greedy_no_improve+1
            else
                greedy_no_improve=0
            end

            ## Then we implement the local search
            ## Here I only keep the 20% best GRASP solutions to apply the local search
            if new_cost_heur<min_cost_heur+paramfixed.windowLocalSearch*min_cost_heur
                new_sol, new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost = manualLocalSearch(inst, new_sol, new_cost_heur, delay_cost_heur, waiting_cost_heur, penalty_cost_heur, handling_cost_heur, fuel_cost_heur, allparam, paramfixed, paramchosen, time_local, cost)
                new_sol.usedLocalSearch=1
		new_sol.store.costLocal=SplitCosts(ceil(Int,new_cost), ceil(Int,delay_cost), ceil(Int,waiting_cost), ceil(Int,penalty_cost), ceil(Int, handling_cost), ceil(Int,fuel_cost))
            else
                new_cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost = new_cost_heur, delay_cost_heur, waiting_cost_heur, penalty_cost_heur, handling_cost_heur, fuel_cost_heur
            end
            
            ## If we actually applied the local search we try a path relinking method
            if paramfixed.pathRelinking=="yes" && new_sol.usedLocalSearch==1 && length(elite_pool)==paramfixed.LengthElite
                ## First we check the distance of the new solution between the minimum distance in the elite pool
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

                ## Then we apply the path relinking method
                start_time_relinking = time_ns()
                elapsed_relinking = round((time_ns()-start_time_relinking)/1e9,digits=3)
                new_sol_path=deepcopy(new_sol)
                delta_distance = deepcopy(distance_btw_sols)
                countpath = 0
                while elapsed_relinking<paramfixed.MaxTimeRelinking && distance_btw_sols>distance_btw_sols_min
                    new_sol.pathcost = new_cost
                    countpath = countpath+1
                    ## We use our pseudo LNS method
                    new_sol_path_loop = greedyremoverandomconstruction(inst, new_sol_path, paramchosen, allparam, paramfixed, true, max_time_heur)
                    feasible=true
                    for n in 1:N
                        # The times :
                        for (c,p) in enumerate(inst.Pi[n])
                            if new_sol.visits[n][c].planned == false || new_sol.visits[n][c].planned == false || sol.visits[n][c].b <= 0 || sol.visits[n][c].t <= 0
                                feasible = false              
                            end
                        end
                    end
                    ## Then if we find a feasible solution we check its cost and distance to the elite pool
                    if feasible                           
                        new_cost_path_loop, delay_cost_path_loop, waiting_cost_path_loop, penalty_cost_path_loop, handling_cost_path_loop, fuel_cost_path_loop = checkSolutionCost(inst, new_sol_path_loop)
                        distance_btw_sols = DistanceSols(inst, new_sol_path_loop, sol_elite)
                        if distance_btw_sols<delta_distance
                            delta_distance=deepcopy(distance_btw_sols)
                            new_sol_path=deepcopy(new_sol_path_loop)
                        end
                        if new_cost_path_loop<new_cost
                            new_sol=deepcopy(new_sol_path)
                            new_cost=deepcopy(new_cost_path_loop)
                            delay_cost=deepcopy(delay_cost_path_loop)
                            waiting_cost=deepcopy(waiting_cost_path_loop)
                            penalty_cost=deepcopy(penalty_cost_path_loop)
                            handling_cost=deepcopy(handling_cost_path_loop)
                            fuel_cost=deepcopy(fuel_cost_path_loop)
                            new_sol.pathcost = deepcopy(new_cost)
                        end
                    end
                    elapsed_relinking = round((time_ns()-start_time_relinking)/1e9,digits=3)
                end
                if distance_btw_sols<distance_btw_sols_min
                    distance_btw_sols_min=deepcopy(distance_btw_sols)
                end
                new_sol.countpath = countpath           
            end
            
            ## Here we update the elite pool if needed (and the minimum distance between the solutions in the elite pool)
            if length(elite_pool)<paramfixed.LengthElite
                distance_sols_min = 10000000000
                for sol_elite in elite_pool
                    distance_sols = DistanceSols(inst, sol_elite[1], new_sol)
                    if distance_sols<distance_sols_min
                        distance_sols_min=deepcopy(distance_sols)
                    end
                end
                if distance_sols_min<min_distance_elite
                    min_distance_elite=deepcopy(distance_sols_min)
                end
                push!(elite_pool, (deepcopy(new_sol),new_cost,distance_sols_min))                   
            else
                elite_pool = sort(elite_pool, by=x->x[2])
                distance_sols_min = 10000000000
                if new_cost<elite_pool[end][2]
                    for sol_elite in elite_pool
                        distance_sols = DistanceSols(inst, sol_elite[1], new_sol)
                        if distance_sols<distance_sols_min
                            distance_sols_min=deepcopy(distance_sols)
                        end
                    end
                    if distance_sols_min>min_distance_elite
                        elite_pool[end]=(deepcopy(new_sol),new_cost, distance_sols_min)
                        elite_pool = sort(elite_pool, by=x->x[3])
                        min_distance_elite=deepcopy(elite_pool[1][3])
                    end  
                end
            end
            list_costs_elite= Vector{}()
            list_dist_elite= Vector{}()
            for el in elite_pool
                push!(list_costs_elite, el[2])
                push!(list_dist_elite, el[3])
            end
            if length(elite_pool)>=paramfixed.LengthElite
                if old_elite_cost_mean==10000000000
                    old_elite_cost_mean=deepcopy(mean(list_costs_elite))
                end
            end
            new_sol.average_cost_elite=mean(list_costs_elite)
            new_sol.average_dist_elite=min_distance_elite

            
            ## At the end we update the parameters 
            allparam = UpdateAfterHeurParameters(paramchosen, allparam, cost, new_cost, paramfixed.lookforconstrained)
            if new_sol.usedLocalSearch==1
                allparam = UpdateAfterLocalParameters(paramchosen, allparam, cost, new_cost)
            end
            
            ## Check if we have a better cost
            new_sol.better=0
            if new_cost<bestcost
                bestcost=deepcopy(new_cost)
                cost=deepcopy(new_cost)
                sol=deepcopy(new_sol)
                new_sol.better=1
            end

            ## and store the solution constructed (for analysis)
            new_sol.store.costHeur=SplitCosts(ceil(Int,new_cost_heur), ceil(Int,delay_cost_heur), ceil(Int,waiting_cost_heur), ceil(Int,penalty_cost_heur), ceil(Int,handling_cost_heur), ceil(Int,fuel_cost_heur))
            new_sol.store.timeHeur=elapsed_heur
            new_sol.store.timeLocalSearch=time_local
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

            d_after = prepareSolIterSoft(seed,N,Nout,qli,nb_iter,inst, new_sol, new_cost, allparam, paramchosen, expname)
            usedcplex=0
            d_alliter_after[nb_iter]=d_after
            nb_iter+=1
            if nb_iter_restart_params>paramfixed.restartParams
                allparam= RestartParamNb(allparam)
            end

                        
                         
        else
            print('\n')
            print("We had unfeasability")
        end
        elapsed = round((time_ns()-start)/1e9,digits=3)
    end


    
    ## Here the last steps for the final solution : check its cost and feasability and then save the solution

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
        for iter in keys(d_alliter_after)
            CSV.write(location*"results_jobs/benchmarks_HEUR/finalGRASP/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli"*"/iter_$iter"*".csv", d_alliter_after[iter])
        end
        nb_iter+=1
        
        d = prepareSolIter(seed,N,Nout,qli,nb_iter,inst, sol, cost, allparam, paramfixed, expname)
        CSV.write(location*"results_jobs/benchmarks_HEUR/finalGRASP/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli"*"/iter_$nb_iter"*".csv", d)
        return sol, cost, allparam
    else
        cost=1000000000
        allparam = initializeParam(paramfixed)
        sol=initializeSol(inst, allparam)
        return sol, cost, allparam
    end
end



