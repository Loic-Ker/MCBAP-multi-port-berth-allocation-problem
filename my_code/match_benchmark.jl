import XLSX
include("neighbours.jl")

function get_scoresN(xf,lookN,lookNout,lookqli)
    newbenchmark = [Any["Seed","N","Nout","qli","LB","UB","Time","NewCost"]]
    for i in 2:3
        this_data=xf[i,1:7]
        seed = this_data[1]
        N = this_data[2]
        Nout = this_data[3]
        qli= this_data[4]
        LB = this_data[5]
        time=this_data[7]
        if N==lookN && qli==lookqli && Nout==lookNout
            print("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_$seed"*"_$N"*"_$Nout"*"_$qli"*".txt")
            instance = readInstFromFile("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_$seed"*"_$N"*"_$Nout"*"_$qli"*".txt")
            this_benchmark = [0,0,0,0,0,0,0,0]
            this_benchmark[1:7]=xf[i,1:7]
            
            sol, cost = GRASP_V1(instance, false, 15, 25)
            if checkSolutionFeasability(instance, sol)
                this_benchmark[8]=floor(Int,cost)
            end
            if checkSolutionFeasability(instance, sol)==false
                this_benchmark[8]="nope"
                push!(newbenchmark,this_benchmark)
            end
        end
    end
    return newbenchmark
end

xf = XLSX.readdata("D:/DTU-Courses/DTU-Thesis/berth_allocation/benchmarks/Small_Inst_Res.xlsx", "Sheet1!A1:G720")
newbenchmark=get_scoresN(xf,4,3,10)
results = mapreduce(permutedims, vcat, newbenchmark)
print(results)
XLSX.writetable("N4_Nout3_qli10.xlsx",
    sheet1=(collect(eachcol(results)), string.("col", 1:8)) )