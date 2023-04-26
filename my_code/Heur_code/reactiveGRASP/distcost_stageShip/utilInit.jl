using UnPack 
using Statistics
import XLSX
using CSV, Tables
using DataFrames
using Random
include("../../../MBAP_INST.jl")
include("../../MBAP_SOL.jl")
include("../../toolsMatrixTimes.jl")
include("../../check_solution.jl")


mutable struct AlphaParameters
    DistancePerShip::Dict{Int64,Vector{Float64}}
    CostPerShip::Dict{Int64,Vector{Float64}}

    DistanceAllShip::Vector{Float64}
    CostAllShip::Vector{Float64}

    RateConstrained::Vector{Float64}
end


mutable struct Probabilities
    ChooseTacticPerShip::Dict{Int64,Vector{Float64}}
    DistancePerShip::Dict{Int64,Vector{Float64}}
    CostPerShip::Dict{Int64,Vector{Float64}}

    ChooseTacticAllShip::Vector{Float64}
    DistanceAllShip::Vector{Float64}
    CostAllShip::Vector{Float64}

    RateConstrained::Vector{Float64}
end

mutable struct AverageSolCost
    ChooseTacticPerShip::Dict{Int64,Vector{Float64}}
    DistancePerShip::Dict{Int64,Vector{Float64}}
    CostPerShip::Dict{Int64,Vector{Float64}}

    ChooseTacticAllShip::Vector{Float64}
    DistanceAllShip::Vector{Float64}
    CostAllShip::Vector{Float64}

    RateConstrained::Vector{Float64}
end

mutable struct NbExperiments
    ChooseTacticPerShip::Dict{Int64,Vector{Int64}}
    DistancePerShip::Dict{Int64,Vector{Int64}}
    CostPerShip::Dict{Int64,Vector{Int64}}

    ChooseTacticAllShip::Vector{Int64}
    DistanceAllShip::Vector{Int64}
    CostAllShip::Vector{Int64}

    RateConstrained::Vector{Float64}
end


mutable struct QParameters
    ChooseTacticPerShip::Dict{Int64,Vector{Float64}}
    DistancePerShip::Dict{Int64,Vector{Float64}}
    CostPerShip::Dict{Int64,Vector{Float64}}

    ChooseTacticAllShip::Vector{Float64}
    DistanceAllShip::Vector{Float64}
    CostAllShip::Vector{Float64}

    RateConstrained::Vector{Float64}
end

mutable struct AdjustProba
    n_proba_pership::Int64
    n_proba_constrained::Int64
end

mutable struct AllParameters
    Alpha::AlphaParameters
    Proba::Probabilities
    AverageCost::AverageSolCost
    Nbexp::NbExperiments
    Q::QParameters
    adjustproba::AdjustProba
end

mutable struct ChosenParameters
    TacticPerShip::Vector{String}
    IndexPerShip::Vector{Int64}
    TacticAllShip::String
    IndexAllShip::Int64
    IndexRateConstrained::Int64
end

function initializeSol(inst::Instance)
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


function initializeParam(inst::Instance, adjustproba::AdjustProba )
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT = inst
    @unpack n_proba_pership, n_proba_constrained = adjustproba

    AlphaVectorsPerShip = [(i-1)/(2*n_proba_pership) for i in 1:n_proba_pership]
    ProbaVectorsPerShip = [i/n_proba_pership for i in 1:n_proba_pership]
    AverageCostVectorsPerShip = [1.0 for i in 1:n_proba_pership]
    NbVectorsPerShip = [0 for i in 1:n_proba_pership]
    QVectorsPerShip = [1.0 for i in 1:n_proba_pership]
    
    AlphaVectorsRateConstrained = [(i-1)/(n_proba_constrained) for i in 1:n_proba_constrained]
    ProbaVectorsRateConstrained = [i/n_proba_constrained for i in 1:n_proba_constrained]
    AverageCostVectorsRateConstrained = [1.0 for i in 1:n_proba_constrained]
    NbVectorsRateConstrained = [0 for i in 1:n_proba_constrained]
    QVectorsRateConstrained = [1.0 for i in 1:n_proba_constrained]

    ProbaVectorsTactic = [1/2,1/2]
    AverageCostVectorsTactic = [1.0,1.0]
    NbVectorsTactic = [0,0]
    QVectorsTactic = [1.0,1.0]
    
    
    AlphaDistancePerShip=Dict{Any,Vector{Float64}}()
    AlphaCostPerShip=Dict{Any,Vector{Float64}}()
    AlphaDistanceAllShip=deepcopy(AlphaVectorsPerShip)
    AlphaCostAllShip=deepcopy(AlphaVectorsPerShip)
    AlphaRateConstrained=deepcopy(AlphaVectorsRateConstrained)

    ProbaChooseTacticPerShip=Dict{Any,Vector{Float64}}()
    ProbaDistancePerShip=Dict{Any,Vector{Float64}}()
    ProbaCostPerShip=Dict{Any,Vector{Float64}}()
    ProbaChooseTacticAllShip=deepcopy(ProbaVectorsTactic)
    ProbaDistanceAllShip=deepcopy(ProbaVectorsPerShip)
    ProbaCostAllShip=deepcopy(ProbaVectorsPerShip)
    ProbaRateConstrained=deepcopy(ProbaVectorsRateConstrained)

    AverageCostChooseTacticPerShip=Dict{Any,Vector{Float64}}()
    AverageCostDistancePerShip=Dict{Any,Vector{Float64}}()
    AverageCostCostPerShip=Dict{Any,Vector{Float64}}()
    AverageCostChooseTacticAllShip=deepcopy(AverageCostVectorsTactic)
    AverageCostDistanceAllShip=deepcopy(AverageCostVectorsPerShip)
    AverageCostCostAllShip=deepcopy(AverageCostVectorsPerShip)
    AverageRateConstrained=deepcopy(AverageCostVectorsRateConstrained)

    NbChooseTacticPerShip=Dict{Any,Vector{Float64}}()
    NbDistancePerShip=Dict{Any,Vector{Float64}}()
    NbCostPerShip=Dict{Any,Vector{Float64}}()
    NbChooseTacticAllShip=deepcopy(NbVectorsTactic)
    NbDistanceAllShip=deepcopy(NbVectorsPerShip)
    NbCostAllShip=deepcopy(NbVectorsPerShip)
    NbRateConstrained=deepcopy(NbVectorsRateConstrained)

    QChooseTacticPerShip=Dict{Any,Vector{Float64}}()
    QDistancePerShip=Dict{Any,Vector{Float64}}()
    QCostPerShip=Dict{Any,Vector{Float64}}()
    QChooseTacticAllShip=deepcopy(QVectorsTactic)
    QDistanceAllShip=deepcopy(QVectorsPerShip)
    QCostAllShip=deepcopy(QVectorsPerShip)
    QRateConstrained=deepcopy(QVectorsRateConstrained)

    for n in 1:inst.N
        AlphaDistancePerShip[n]=deepcopy(AlphaVectorsPerShip)
        AlphaCostPerShip[n]=deepcopy(AlphaVectorsPerShip)

        ProbaChooseTacticPerShip[n]=deepcopy(ProbaVectorsTactic)
        ProbaDistancePerShip[n]=deepcopy(ProbaVectorsPerShip)
        ProbaCostPerShip[n]=deepcopy(ProbaVectorsPerShip)

        AverageCostChooseTacticPerShip[n]=deepcopy(AverageCostVectorsTactic)
        AverageCostDistancePerShip[n]=deepcopy(AverageCostVectorsPerShip)
        AverageCostCostPerShip[n]=deepcopy(AverageCostVectorsPerShip)

        NbChooseTacticPerShip[n]=deepcopy(NbVectorsTactic)
        NbDistancePerShip[n]=deepcopy(NbVectorsPerShip)
        NbCostPerShip[n]=deepcopy(NbVectorsPerShip)

        QChooseTacticPerShip[n]=deepcopy(QVectorsTactic)
        QDistancePerShip[n]=deepcopy(QVectorsPerShip)
        QCostPerShip[n]=deepcopy(QVectorsPerShip)
    end

    AlphaAll = AlphaParameters(AlphaDistancePerShip,AlphaCostPerShip,AlphaDistanceAllShip,AlphaCostAllShip, AlphaRateConstrained)
    ProbaAll = Probabilities(ProbaChooseTacticPerShip,ProbaDistancePerShip,ProbaCostPerShip,ProbaChooseTacticAllShip,ProbaDistanceAllShip,ProbaCostAllShip, ProbaRateConstrained)
    AverageCostAll = AverageSolCost(AverageCostChooseTacticPerShip,AverageCostDistancePerShip,AverageCostCostPerShip,AverageCostChooseTacticAllShip,AverageCostDistanceAllShip,AverageCostCostAllShip, AverageRateConstrained)
    NbAll = NbExperiments(NbChooseTacticPerShip,NbDistancePerShip,NbCostPerShip,NbChooseTacticAllShip,NbDistanceAllShip,NbCostAllShip, NbRateConstrained)
    QAll = QParameters(QChooseTacticPerShip,QDistancePerShip,QCostPerShip,QChooseTacticAllShip,QDistanceAllShip,QCostAllShip, QRateConstrained)
    
    All = AllParameters(AlphaAll,ProbaAll,AverageCostAll,NbAll,QAll,adjustproba)
    
    return All
end



