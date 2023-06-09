
include("../MBAP_INST.jl")
import XLSX
## Package to save the dict
using CSV, Tables
using DataFrames

## Packages needed for solving MIPs
using JuMP, CPLEX
#m = Model(GLPK.Optimizer)
#m = Model(HiGHS.Optimizer)

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

function doOverlapRectangleCorrect(x1l, x1r, y1d, y1u, x2l, x2r, y2d, y2u)
    # check if one element is a line
    if x1l == x1r || y1u == y1d || x2l == x2r || y2u == y2d
        return false, ""
    end
    # If one rectangle is on left side of other
    if x1l >= x2r || x2l >= x1r
        return false, ""
    end
    # If one rectangle is above other
    if y1d >= y2u || y2d >= y1u
        return false, ""
    end
    if abs(y1d-y2d)==1
        y1d=y2u
    end
    if abs(y1d-y2u)==1
        to_change="Change x1"
    end
    if abs(y2d-y1u)==1
        to_change="Change x2"
    end
    return true, to_change
end


function CPLEXoptimize(N,Nout,seed,qli, time, location)
    m = Model(CPLEX.Optimizer)
    set_optimizer_attribute(m, "CPX_PARAM_EPINT", 1e-8)
    set_optimizer_attribute(m, "CPX_PARAM_TILIM", time)
    set_optimizer_attribute(m, "CPXPARAM_Threads", 1)


    #inst = readInstFromFile("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_$seed"*"_$N"*"_$Nout"*"_$qli"*".txt")
    inst = readInstFromFile(location*"MCBAP-multi-port-berth-allocation-problem/Large/CP2_Inst_$seed"*"_$N"*"_$Nout"*"_$qli"*".txt")
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

        @variable(m,x[n=1:Ntot, 1:length(inst.Pi[n])] >= 0, Int)
        @variable(m,y[n=1:Ntot, 1:length(inst.Pi[n])] >= 0)

        ## The ones to check the positions
        @variable(m,sig[n1=1:Ntot,n2=1:Ntot, 1:length(inst.Pi[n1]), 1:length(inst.Pi[n2])], Bin)
        @variable(m,del[n1=1:Ntot,n2=1:Ntot, 1:length(inst.Pi[n1]), 1:length(inst.Pi[n2])], Bin)
        @variable(m,hand[n=1:Ntot, 1:length(inst.Pi[n])] >= 0)

        
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
        #heur_results = CSV.File("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/greedy_only/sols/HEUR_sol_$seed"*"_$N"*"_$Nout"*"_$qli"*".csv") |> Dict
        #x_heur = eval(Meta.parse(heur_results["x"]))
        #y_heur = eval(Meta.parse(heur_results["y"]))
        #hand_heur = eval(Meta.parse(heur_results["hand"]))
        #a_heur = eval(Meta.parse(heur_results["a"]))
        #v_heur = eval(Meta.parse(heur_results["v"]))

        #for n in 1:N
        #    for c in 1:length(inst.Pi[n])
        #        @constraint(m,x[n,c]==x_heur[n][c])
        #        @constraint(m,y[n,c]==y_heur[n][c])
        #        @constraint(m,hand[n,c]==hand_heur[n][c])
        #        @constraint(m,a[n,c]==a_heur[n][c])
        #    end
        #end

        #speeds = length(delta)
        #for n in 1:N
        #    if length(inst.Pi[n])>1
        #        l=length(inst.Pi[n])-1
        #        for c in 1:l
        #            for s in 1:speeds
        #                if v_heur[n][c]==s
        #                    @constraint(m, v[n,c,s]==1)
        #                end
        #            end
        #        end
        #    end
        #end

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
        @constraint(m,c11[n=1:N, c=1:length(inst.Pi[n])],   x[n,c] - (inst.shipsIn[n].Bi[c]/qli) <= r[n,c])
        @constraint(m,c12[n=1:N, c=1:length(inst.Pi[n])],   (shipsIn[n].Bi[c]/qli)-x[n,c] <= r[n,c])

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

    



function makeSoltest(time, location)
    newbenchmark = DataFrame(Seed= [""],N= [""],Nout= [""],qli= [""], Time= [0], CPLEX= [0], Box= [""]) #HeurCost= [0],
    all_instances = readdir(location*"MCBAP-multi-port-berth-allocation-problem/Large")
    for instance_name in all_instances[1:2]
        split_instance = split(instance_name,"_")
        seed=split_instance[3]
        N=split_instance[4]
        Nout=split_instance[5]
        qli=split(split_instance[6],".")[1]
        print("The instance : $seed"*"_$N"*"_$Nout"*"_$qli")
        start = time_ns()
        box, d, cost = CPLEXoptimize(N,Nout,seed,qli, time, location) 
        elapsed = ceil(Int, round((time_ns()-start)/1e9,digits=3))
        print('\n')
        print(cost)
        print('\n')
        CSV.write(location*"results_jobs/benchmarks_CPLEX/sols_5min/CPLEX_sol_$seed"*"_$N"*"_$Nout"*"_$qli"*".csv", d)
        #CSV.write("D:/DTU-Courses/DTU-Thesis/berth_allocation/MCBAP-multi-port-berth-allocation-problem/results_jobs/benchmarks_CPLEX/sols/CPLEX_sol_$seed"*"_$N"*"_$Nout"*"_$qli"*".csv", d)
        this_benchmark=DataFrame(Seed= [seed],N= [N],Nout= [Nout],qli= [qli], Time= [elapsed], CPLEX=[ceil(Int, cost)],  Box= [box]) #HeurCost= [costHeur],
        newbenchmark=append!(newbenchmark,this_benchmark)
    end
    return newbenchmark
end

location = "D:/DTU-Courses/DTU-Thesis/berth_allocation/"
#location="/zhome/c3/6/164957/code_git/"
#time = parse(Int64,ARGS[1])
#minN = 8
#maxN = 8
time = 300
newbenchmark = makeSoltest(time, location)
CSV.write(location*"results_jobs/benchmarks_CPLEX/CPLEX_NLarge_results_5min.csv", newbenchmark)
    

## At each iteration :
    # Get the different cost right after greedy heurisitics
    # And after the local searh and compare
    # Compare also solutions visually before and after ?
    # Get successfull and not successfull iteration
    # Get the probability of each coefficient
    # For each visit : when was it added ? What was the tactic ? What was the cost to add it ?

## At the end : 
    # Get the cost of each visit
    # Get the cost for each type of boat 
