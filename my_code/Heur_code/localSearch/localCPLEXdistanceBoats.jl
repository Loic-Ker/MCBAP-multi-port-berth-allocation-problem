using UnPack 
using Statistics
import XLSX
using CSV, Tables
using DataFrames
using Random
using StatsBase
## Packages needed for solving MIPs
using JuMP, CPLEX
#m = Model(GLPK.Optimizer)
#m = Model(HiGHS.Optimizer)
include("../../MBAP_INST.jl")
include("../MBAP_SOL.jl")
include("../toolsMatrixTimes.jl")
include("../check_solution.jl")

function CPLEXoptimizeLocalSearch(inst::Instance, sol::Sol, not_taken)
    m = Model(CPLEX.Optimizer)
    set_optimizer_attribute(m, "CPX_PARAM_EPINT", 1e-8)

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
                        if doOverlapRectangle(t1, t2, x1, x2, t1_, t2_, x1_, x2_)
                            txt_error = txt_error*"####"*"The port $p, the first box : [$t1,$t2] and [$x1,$x2], the second box : [$t1_,$t2_] and [$x1_,$x2_]"
                            print('\n')
                            print("The port $p")
                            print('\n')
                            print("The first box : [$t1,$t2] and [$x1,$x2]")
                            print('\n')
                            print("The second box : [$t1_,$t2_] and [$x1_,$x2_]")
                            print('\n')
                        end
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
                    @constraint(m,x[n,c]==sol.visits[n][c].b)
                    @constraint(m,y[n,c]==sol.visits[n][c].t)
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
        
        println("Objective: ", JuMP.objective_value(m)) #Printing the objective value
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

function local_search(inst::Instance, sol::Sol, alpha, max_time)
    @unpack N, Ntot, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT, Nl, gamma, Hc, Dc, Fc, Ic, Pc, beta, ports = inst
    new_sol=deepcopy(sol)
    new_cost=checkSolutionCost(inst, new_sol)
    list_initialize = []
    for n in 1:N
        append!(list_initialize,[n])
    end
    start = time_ns()
    elapsed = round((time_ns()-start)/1e9,digits=3)
    while elapsed<max_time && length(list_initialize)>0.2*length(inst.Pi)
        new_boat=rand(list_initialize)
        other_visits=get_closest(inst,new_sol,new_boat,alpha)
        new_solCPLEX=deepcopy(new_sol)
        txt_error, d, cost = CPLEXoptimizeLocalSearch(inst, new_sol, other_visits)
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
            print("The new boat :")
            print('\n')
            print(new_boat)
            print('\n')
            print("During local search")
            print('\n')
            print(new_solCPLEX.visits)
            new_sol=deepcopy(new_solCPLEX)
            new_cost=checkSolutionCost(inst, new_sol)
            list_initialize=setdiff(list_initialize, [new_boat])
            print('\n')
            print(list_initialize)
        end
        elapsed = round((time_ns()-start)/1e9,digits=3)
    end
    
    return new_sol, new_cost
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

#print('\n')
#print("First cost solution")
#print('\n')
#print(cost)
#print('\n')


#print('\n')
#print("Second cost solution")
#print('\n')
#print(cost2)
#print('\n')














