include("reactive_free.jl")


#instance = readInstFromFile("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_2_4_3_10.txt")
#sol = greedyrandomizedconstruction(instance, false, 5)

function printSol(instance,  type1, type2, adjustproba, alphaboat, alpharandom, time_local, max_time_heur, max_time)
    sol, cost, allparam = GRASP_reactive(instance, type1,type2, adjustproba, alphaboat, alpharandom, time_local, max_time_heur, max_time)

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



function makeSolHeur(type1, type2, adjustproba, alphaboat, alpharandom, time_local, max_time_heur, max_time)
    xf = CSV.read("D:/DTU-Courses/DTU-Thesis/berth_allocation/bernardo_bench/Small_Inst_Res.csv", DataFrame)
    newbenchmark = DataFrame(Seed= [0],N= [0],Nout= [0],qli= [0],OldLB= [0],OldUB= [0],OldTime= [0],HeurCost= [0])
    for N in 5:10
        for qli in [10,20,40,80]
            for Nout in 3:5
                for seed in 1:5
                    print("The instance : $seed"*"_$N"*"_$Nout"*"_$qli")
                    inst = readInstFromFile("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_$seed"*"_$N"*"_$Nout"*"_$qli"*".txt")
                    cost, d = printSol(inst,  type1, type2, adjustproba, alphaboat, alpharandom, time_local, max_time_heur, max_time)
                    CSV.write("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/reactiveGRASP/countcost_nostage/sols/HEUR_LOCAL_exp2_sol_$seed"*"_$N"*"_$Nout"*"_$qli"*".csv", d)
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

newbenchmark = makeSolHeur("both","both", AdjustProba(3,3,3,3,3,3), 5, 18, 20, 20, 120)
CSV.write("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/reactiveGRASP/countcost_nostage/HEUR_LOCAL_N5_N7_exp2.csv", newbenchmark)


