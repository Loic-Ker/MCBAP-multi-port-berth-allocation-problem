include("reactive_free.jl")


#instance = readInstFromFile("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_2_4_3_10.txt")
#sol = greedyrandomizedconstruction(instance, false, 5)


function prepareSol(inst, sol, cost)
    @unpack N, Ntot, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT, Nl = inst
    d = Dict{Any,Any}()
    
    ### Here add: 
        #the type of the boats with every visit
        #The cost divided in total and for each visit added
        #The position of each visit added
        # Parameters at the end of the solution make
        # Which tactic for each boat


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
                if s!=-1
                    push!(a_thisboat, tN + h[n][c-1][bN] + delta[s]*dist[pN,sol.visits[n][c].p])
                end
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


function makeExpText(type1, type2, adjustproba, alphaboat, alpharandom, time_local, max_time_heur, max_time, expname)
    filename = "D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/reactiveGRASP/countcost_nostage/$expname"*"/explanations.txt"
    txttype1="The tactic for each ship : "*type1
    txttype2="The tactic between the ships : "*type2
    @unpack n_proba_typeship_distance, n_proba_typeship_cost, n_proba_all_count, n_proba_all_cost, n_proba_constrained, n_proba_local = adjustproba
    txttype1n_proba = "For the tactic for each ship, the number of initial probabilities is $n_proba_typeship_distance for the distance and $n_proba_typeship_cost for the cost"
    txttype2n_proba = "For the tactic for all ships, the number of initial probabilities is $n_proba_all_count for the count and $n_proba_all_count for the cost"
    txtalphaboat = "For the local search with the boats we take $alphaboat boats"
    txtalpharandom = "For the local search with random visits we take $alpharandom visits"
    txttimelocal = "The maximum time for the local search : $time_local "
    txttimeheur = "The maximum time for the heuristic : $max_time_heur "
    txttimemax = "The time of the exp : $max_time"


    open(filename, "w") do file
        write(file, txttype1)
        write(file,'\n')
        write(file, txttype2)
        write(file,'\n')
        write(file, txttype1n_proba)
        write(file,'\n')
        write(file, txttype2n_proba)
        write(file,'\n')
        write(file, txtalphaboat)
        write(file,'\n')
        write(file, txtalpharandom)
        write(file,'\n')
        write(file, txttimeheur)
        write(file,'\n')
        write(file, txttimelocal)
        write(file,'\n')
        write(file, txttimemax)
        write(file,'\n')
    end
end


function makeSolHeur(type1, type2, adjustproba, alphaboat, alpharandom, time_local, max_time_heur, max_time, expname)
    xf = CSV.read("D:/DTU-Courses/DTU-Thesis/berth_allocation/bernardo_bench/Small_Inst_Res.csv", DataFrame)
    newbenchmark = DataFrame(Seed= [0],N= [0],Nout= [0],qli= [0],OldLB= [0],OldUB= [0],OldTime= [0],HeurCost= [0])
    for N in 10:10
        for qli in [10]#,20,40,80]
            for Nout in 5:5
                for seed in 1:5
                    inst = readInstFromFile("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_$seed"*"_$N"*"_$Nout"*"_$qli"*".txt")
                    print("The instance : $seed"*"_$N"*"_$Nout"*"_$qli")
                    if isdir("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/reactiveGRASP/countcost_nostage/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli")==false
                        mkdir("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/reactiveGRASP/countcost_nostage/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli")
                    end
                    sol, cost, allparam = GRASP_reactive(seed,N,Nout,qli,type1, type2, adjustproba, alphaboat, alpharandom, time_local, max_time_heur, max_time, expname)
                    print('\n')
                    print("The solution :")
                    print('\n')
                    print(sol.visits)
                    print('\n')
                    print("And the cost is ")
                    print('\n')
                    print(cost)
                    d=prepareSol(inst, sol, cost)
                    CSV.write("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/reactiveGRASP/countcost_nostage/$expname"*"/final_sols/sol_$seed"*"_$N"*"_$Nout"*"_$qli"*".csv", d)
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


expname="exp1"
type1="both" 
type2="both" 
adjustproba=AdjustProba(3,3,3,3,3,3)
alphaboat=5
alpharandom=18
time_local=20
max_time_heur=120
max_time=20
makeExpText(type1, type2, adjustproba, alphaboat, alpharandom, time_local, max_time_heur, max_time, expname)
newbenchmark = makeSolHeur(type1, type2, adjustproba, alphaboat, alpharandom, time_local, max_time_heur, max_time, expname)
CSV.write("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/reactiveGRASP/countcost_nostage/$expname"*"/N5_qli10_Nout5.csv", newbenchmark)


