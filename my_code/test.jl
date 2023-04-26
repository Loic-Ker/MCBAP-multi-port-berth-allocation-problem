include("MBAP_INST.jl")
using CSV, Tables
using DataFrames
inst = readInstFromFile("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_1_4_3_10.txt")
@unpack N, Ntot, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT, Nl, gamma, Hc, Dc, Fc, Ic, Pc, beta, ports = inst
heur_results = CSV.File("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/greedy_only/sols/HEUR_sol_1_4_4_10.csv") |> Dict
x_heur = eval(Meta.parse(heur_results["x"]))
y_heur = eval(Meta.parse(heur_results["y"]))
print('\n')
print(y_heur)
hand_heur = eval(Meta.parse(heur_results["hand"]))
print('\n')
print(hand_heur)
a_heur = eval(Meta.parse(heur_results["a"]))
print('\n')
print(a_heur)
v_heur = eval(Meta.parse(heur_results["v"]))

# compute the lowest feasible speed to satisfy travel time t
function findLowestSpeed(t, delta::Vector{Float64}, dist::Int64)
    S = length(delta)
    for s in 1:S
        tt = delta[s]*dist
        if tt < t
            return s
        end
    end
    return -1
end


for n in 1:N
    for (c,p) in enumerate(inst.Pi[n])
        if c>1
            deltaT = y_heur[n][c] - (y_heur[n][c-1] + hand_heur[n][c-1])
            s = findLowestSpeed(deltaT, delta, dist[inst.Pi[n][c-1],inst.Pi[n][c]])
            print('\n')
            print("#################")
            print('\n')
            print("Arrival times")
            print('\n')
            print(y_heur[n][c-1] + hand_heur[n][c-1] + delta[s]*dist[inst.Pi[n][c-1],inst.Pi[n][c]])
            print('\n')
            print(a_heur[n][c])
            print('\n')
            print(y_heur[n][c])
            print('\n')
            print("Speeds")
            print('\n')
            print(s)
            print('\n')
            print(v_heur[n][c-1])
        end
    end
end

newbenchmark = DataFrame(Seed= [0],N= [0],Nout= [0],qli= [0],OldLB= [0],OldUB= [0],OldTime= [0],HeurCost= [0])
this_benchmark=DataFrame(Seed= [1],N= [1],Nout= [1],qli= [1],OldLB= [1],OldUB= [1],OldTime= [1],HeurCost= [1])
newbenchmark=append!(newbenchmark,this_benchmark)
xf=CSV.read("D:/DTU-Courses/DTU-Thesis/berth_allocation/bernardo_bench/Small_Inst_Res.csv", DataFrame)
filtering = xf[(xf.Seed.==1) .& (xf.N.==4) .& (xf.qli.==10) .& (xf.Nout.==3),:]
CSV.write("D:/DTU-Courses/DTU-Thesis/berth_allocation/bernardo_bench/looooool.csv", filtering)
#CSV.write("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks_HEUR/greedy_only/test.csv",  Tables.table(newbenchmark), writeheader=false)