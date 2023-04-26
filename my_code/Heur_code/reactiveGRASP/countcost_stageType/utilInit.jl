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
    DistanceTypeShip::Dict{Int64,Vector{Float64}}
    CostTypeShip::Dict{Int64,Vector{Float64}}

    CountAllShip::Vector{Float64}
    CostAllShip::Vector{Float64}
    RateConstrained::Vector{Float64}
    
end


mutable struct Probabilities
    ChooseTacticTypeShip::Dict{Int64,Vector{Float64}}
    DistanceTypeShip::Dict{Int64,Vector{Float64}}
    CostTypeShip::Dict{Int64,Vector{Float64}}

    ChooseTacticAllShip::Vector{Float64}
    CountAllShip::Vector{Float64}
    CostAllShip::Vector{Float64}

    RateConstrained::Vector{Float64}
end

mutable struct AverageSolCost
    ChooseTacticTypeShip::Dict{Int64,Vector{Float64}}
    DistanceTypeShip::Dict{Int64,Vector{Float64}}
    CostTypeShip::Dict{Int64,Vector{Float64}}

    ChooseTacticAllShip::Vector{Float64}
    CountAllShip::Vector{Float64}
    CostAllShip::Vector{Float64}

    RateConstrained::Vector{Float64}
end

mutable struct NbExperiments
    ChooseTacticTypeShip::Dict{Int64,Vector{Int64}}
    DistanceTypeShip::Dict{Int64,Vector{Float64}}
    CostTypeShip::Dict{Int64,Vector{Float64}}

    ChooseTacticAllShip::Vector{Int64}
    CountAllShip::Vector{Int64}
    CostAllShip::Vector{Int64}

    RateConstrained::Vector{Float64}
end


mutable struct QParameters
    ChooseTacticTypeShip::Dict{Int64,Vector{Float64}}
    DistanceTypeShip::Dict{Int64,Vector{Float64}}
    CostTypeShip::Dict{Int64,Vector{Float64}}

    ChooseTacticAllShip::Vector{Float64}
    CountAllShip::Vector{Float64}
    CostAllShip::Vector{Float64}

    RateConstrained::Vector{Float64}
end

mutable struct AdjustProba
    n_proba_typeship_distance::Int64
    n_proba_typeship_cost::Int64
    n_proba_all_count::Int64
    n_proba_all_cost::Int64
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
    TacticTypeShip::Vector{String}
    IndexTypeShip::Vector{Int64}
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
    @unpack n_proba_typeship_distance, n_proba_typeship_cost, n_proba_all_count, n_proba_all_cost, n_proba_constrained = adjustproba

    AlphaVectorsTypeShipDistance = [(i-1)/(n_proba_typeship_distance) for i in 1:n_proba_typeship_distance+1]
    ProbaVectorsTypeShipDistance = [1/(n_proba_typeship_distance+1) for i in 1:n_proba_typeship_distance+1]
    AverageCostVectorsTypeShipDistance = [1.0 for i in 1:n_proba_typeship_distance+1]
    NbVectorsTypeShipDistance = [0 for i in 1:n_proba_typeship_distance+1]
    QVectorsTypeShipDistance = [1.0 for i in 1:n_proba_typeship_distance+1]

    AlphaVectorsTypeShipCost = [(i-1)/(n_proba_typeship_cost) for i in 1:n_proba_typeship_cost+1]
    ProbaVectorsTypeShipCost = [1/(n_proba_typeship_cost+1) for i in 1:n_proba_typeship_cost+1]
    AverageCostVectorsTypeShipCost = [1.0 for i in 1:n_proba_typeship_cost+1]
    NbVectorsTypeShipCost = [0 for i in 1:n_proba_typeship_cost+1]
    QVectorsTypeShipCost = [1.0 for i in 1:n_proba_typeship_cost+1]

    AlphaVectorsAllCost = [(i-1)/(n_proba_all_cost) for i in 1:n_proba_all_cost+1]
    ProbaVectorsAllCost = [1/(n_proba_all_cost+1) for i in 1:n_proba_all_cost+1]
    AverageCostVectorsAllCost = [1.0 for i in 1:n_proba_all_cost+1]
    NbVectorsAllCost = [0 for i in 1:n_proba_all_cost+1]
    QVectorsAllCost = [1.0 for i in 1:n_proba_all_cost+1]


    AlphaVectorsAllCount = [(i-1)/(n_proba_all_count) for i in 1:n_proba_all_count+1]
    ProbaVectorsAllCount = [1/(n_proba_all_count+1) for i in 1:n_proba_all_count+1]
    AverageCostVectorsAllCount = [1.0 for i in 1:n_proba_all_count+1]
    NbVectorsAllCount = [0 for i in 1:n_proba_all_count+1]
    QVectorsAllCount = [1.0 for i in 1:n_proba_all_count+1]
    
    AlphaVectorsRateConstrained = [(i-1)/(n_proba_constrained) for i in 1:n_proba_constrained+1]
    ProbaVectorsRateConstrained = [1/(n_proba_constrained+1) for i in 1:n_proba_constrained+1]
    AverageCostVectorsRateConstrained = [1.0 for i in 1:n_proba_constrained+1]
    NbVectorsRateConstrained = [0 for i in 1:n_proba_constrained+1]
    QVectorsRateConstrained = [1.0 for i in 1:n_proba_constrained+1]

    ProbaVectorsTactic = [1/2,1/2]
    AverageCostVectorsTactic = [1.0,1.0]
    NbVectorsTactic = [0,0]
    QVectorsTactic = [1.0,1.0]
    
    
    AlphaDistanceTypeShip=Dict{Any,Vector{Float64}}()
    AlphaCostTypeShip=Dict{Any,Vector{Float64}}()
    AlphaCountAllShip=deepcopy(AlphaVectorsAllCount)
    AlphaCostAllShip=deepcopy(AlphaVectorsAllCost)
    AlphaRateConstrained=deepcopy(AlphaVectorsRateConstrained)

    ProbaChooseTacticTypeShip=Dict{Any,Vector{Float64}}()
    ProbaDistanceTypeShip=Dict{Any,Vector{Float64}}()
    ProbaCostTypeShip=Dict{Any,Vector{Float64}}()
    ProbaChooseTacticAllShip=deepcopy(ProbaVectorsTactic)
    ProbaCountAllShip=deepcopy(ProbaVectorsAllCount)
    ProbaCostAllShip=deepcopy(ProbaVectorsAllCost)
    ProbaRateConstrained=deepcopy(ProbaVectorsRateConstrained)

    AverageCostChooseTacticTypeShip=Dict{Any,Vector{Float64}}()
    AverageCostDistanceTypeShip=Dict{Any,Vector{Float64}}()
    AverageCostCostTypeShip=Dict{Any,Vector{Float64}}()
    AverageCostChooseTacticAllShip=deepcopy(AverageCostVectorsTactic)
    AverageCostCountAllShip=deepcopy(AverageCostVectorsAllCount)
    AverageCostCostAllShip=deepcopy(AverageCostVectorsAllCost)
    AverageRateConstrained=deepcopy(AverageCostVectorsRateConstrained)

    NbChooseTacticTypeShip=Dict{Any,Vector{Float64}}()
    NbDistanceTypeShip=Dict{Any,Vector{Float64}}()
    NbCostTypeShip=Dict{Any,Vector{Float64}}()
    NbChooseTacticAllShip=deepcopy(NbVectorsTactic)
    NbCountAllShip=deepcopy(NbVectorsAllCount)
    NbCostAllShip=deepcopy(NbVectorsAllCost)
    NbRateConstrained=deepcopy(NbVectorsRateConstrained)

    QChooseTacticTypeShip=Dict{Any,Vector{Float64}}()
    QDistanceTypeShip=Dict{Any,Vector{Float64}}()
    QCostTypeShip=Dict{Any,Vector{Float64}}()
    QChooseTacticAllShip=deepcopy(QVectorsTactic)
    QCountAllShip=deepcopy(QVectorsAllCount)
    QCostAllShip=deepcopy(QVectorsAllCost)
    QRateConstrained=deepcopy(QVectorsRateConstrained)

    for n in 1:3
        AlphaDistanceTypeShip[n]=deepcopy(AlphaVectorsTypeShipDistance)
        AlphaCostTypeShip[n]=deepcopy(AlphaVectorsTypeShipCost)

        ProbaChooseTacticTypeShip[n]=deepcopy(ProbaVectorsTactic)
        ProbaDistanceTypeShip[n]=deepcopy(ProbaVectorsTypeShipDistance)
        ProbaCostTypeShip[n]=deepcopy(ProbaVectorsTypeShipCost)

        AverageCostChooseTacticTypeShip[n]=deepcopy(AverageCostVectorsTactic)
        AverageCostDistanceTypeShip[n]=deepcopy(AverageCostVectorsTypeShipDistance)
        AverageCostCostTypeShip[n]=deepcopy(AverageCostVectorsTypeShipCost)

        NbChooseTacticTypeShip[n]=deepcopy(NbVectorsTactic)
        NbDistanceTypeShip[n]=deepcopy(NbVectorsTypeShipDistance)
        NbCostTypeShip[n]=deepcopy(NbVectorsTypeShipCost)

        QChooseTacticTypeShip[n]=deepcopy(QVectorsTactic)
        QDistanceTypeShip[n]=deepcopy(QVectorsTypeShipDistance)
        QCostTypeShip[n]=deepcopy(QVectorsTypeShipCost)
    end

    AlphaAll = AlphaParameters(AlphaDistanceTypeShip,AlphaCostTypeShip,AlphaCountAllShip,AlphaCostAllShip,AlphaRateConstrained)
    ProbaAll = Probabilities(ProbaChooseTacticTypeShip,ProbaDistanceTypeShip,ProbaCostTypeShip,ProbaChooseTacticAllShip,ProbaCountAllShip,ProbaCostAllShip, ProbaRateConstrained)
    AverageCostAll = AverageSolCost(AverageCostChooseTacticTypeShip,AverageCostDistanceTypeShip,AverageCostCostTypeShip,AverageCostChooseTacticAllShip,AverageCostCountAllShip,AverageCostCostAllShip, AverageRateConstrained)
    NbAll = NbExperiments(NbChooseTacticTypeShip,NbDistanceTypeShip,NbCostTypeShip,NbChooseTacticAllShip,NbCountAllShip,NbCostAllShip, NbRateConstrained)
    QAll = QParameters(QChooseTacticTypeShip,QDistanceTypeShip,QCostTypeShip,QChooseTacticAllShip,QCountAllShip,QCostAllShip, QRateConstrained)
    
    All = AllParameters(AlphaAll,ProbaAll,AverageCostAll,NbAll,QAll,adjustproba)
    
    return All
end



