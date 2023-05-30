using UnPack 
using Statistics
import XLSX
using CSV, Tables
using DataFrames
using Random
include("MBAP_SOL.jl")
include("toolsMatrixTimes.jl")
include("check_solution.jl")


mutable struct AverageSolCost
    TacticOneBoat::Vector{Float64}
    TacticAllBoats::Vector{Float64}
    TacticLocalSearch::Vector{Float64}
end

mutable struct NbExperiments
    TacticOneBoat::Vector{Float64}
    TacticAllBoats::Vector{Float64}
    TacticLocalSearch::Vector{Float64}
end


mutable struct QParameters
    TacticOneBoat::Vector{Float64}
    TacticAllBoats::Vector{Float64}
    TacticLocalSearch::Vector{Float64}
end

mutable struct AllParameters
    Proba::Probabilities
    AverageCost::AverageSolCost
    Nbexp::NbExperiments
    Q::QParameters
end

mutable struct ChosenParameters
    TacticOneBoat::String
    TacticAllBoats::String
    TacticLocalSearch::String
end

mutable struct FixedParameters
    OneBoatCost::Float64
    OneBoatDistance::Float64
    OneBoatTime::Float64
    AllBoatsCost::Float64
    AllBoatsCount::Float64
    AllBoatsTime::Float64
    LocalSearchRandom::Int64
    LocalSearchBoat::Int64
end

function initializeSol(inst::Instance)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT = inst
    sol = Sol(inst)
    # solution = [Vector{Tuple{Int64,Int64,Int64}}() for n in 1:N] # (p,b,t)
    for n in 1:N, (c,p) in enumerate(Pi[n])
        sol.visits[n][c].minT = max(shipsIn[n].sT[c], T[n,c,1])
        sol.visits[n][c].maxT = maxT
    end
    sol.M = generateOccupiedMx(inst, sol.visits)
    #sol = updateTimesInitialization(inst,sol)
    return sol
end


function initializeParam(inst::Instance)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT = inst

    AverageTacticOneBoat = [1.0,1.0,1.0]
    AverageTacticAllBoats = [1.0,1.0,1.0]
    AverageTacticLocalSearch = [1.0,1.0]

    NbTacticOneBoat = [0,0,0]
    NbTacticAllBoats = [0,0,0]
    NbTacticLocalSearch = [0,0]

    QTacticOneBoat = [1.0,1.0,1.0]
    QTacticAllBoats = [1.0,1.0,1.0]
    QTacticLocalSearch = [1.0,1.0]

    ProbaTacticOneBoat = [1/3,1/3,1/3]
    ProbaTacticAllBoats = [1/3,1/3,1/3]
    ProbaTacticLocalSearch = [1/2,1/2]

    AverageCostAll = AverageSolCost(AverageTacticOneBoat,AverageTacticAllBoats,AverageTacticLocalSearch)
    NbAll = NbExperiments(NbTacticOneBoat,NbTacticAllBoats,NbTacticLocalSearch)
    QAll = QParameters(QTacticOneBoat,QTacticAllBoats,QTacticLocalSearch)
    ProbaAll = Probabilities(ProbaTacticOneBoat,ProbaTacticAllBoats,ProbaTacticLocalSearch)

    All = AllParameters(ProbaAll,AverageCostAll,NbAll,QAll)
    
    return All
end



