using UnPack 
using Statistics
import XLSX
using CSV, Tables
using DataFrames
using Random
include("MBAP_SOL.jl")
include("toolsMatrixTimes.jl")
include("check_solution.jl")


mutable struct AlphaParameters
    DistOneShip::Vector{Float64}
    CostOneShip::Vector{Float64}
    TimeOneShip::Vector{Float64}

    TimeAllShip::Vector{Float64}
    CountAllShip::Vector{Float64}
    CostAllShip::Vector{Float64}

    RateConstrained::Vector{Float64}
    
end

mutable struct AverageSolCost
    TacticOneBoat::Vector{Float64}
    TacticAllBoats::Vector{Float64}
    TacticLocalSearch::Vector{Float64}

    TimeOneShip::Vector{Float64}
    DistOneShip::Vector{Float64}
    CostOneShip::Vector{Float64}

    TimeAllShip::Vector{Float64}
    CountAllShip::Vector{Float64}
    CostAllShip::Vector{Float64}

    RateConstrained::Vector{Float64}

end

mutable struct NbExperiments
    TacticOneBoat::Vector{Float64}
    TacticAllBoats::Vector{Float64}
    TacticLocalSearch::Vector{Float64}

    TimeOneShip::Vector{Float64}
    DistOneShip::Vector{Float64}
    CostOneShip::Vector{Float64}

    TimeAllShip::Vector{Float64}
    CountAllShip::Vector{Float64}
    CostAllShip::Vector{Float64}

    RateConstrained::Vector{Float64}
end


mutable struct QParameters
    TacticOneBoat::Vector{Float64}
    TacticAllBoats::Vector{Float64}
    TacticLocalSearch::Vector{Float64}

    TimeOneShip::Vector{Float64}
    DistOneShip::Vector{Float64}
    CostOneShip::Vector{Float64}

    TimeAllShip::Vector{Float64}
    CountAllShip::Vector{Float64}
    CostAllShip::Vector{Float64}

    RateConstrained::Vector{Float64}
end

mutable struct AllParameters
    Proba::Probabilities
    AverageCost::AverageSolCost
    Nbexp::NbExperiments
    Q::QParameters
    Alpha::AlphaParameters
end

mutable struct ChosenParameters
    TacticOneBoat::String
    TacticAllBoats::String
    TacticLocalSearch::String
    IndexOneShip::Int64
    IndexAllShip::Int64
    IndexRateConstrained::Int64
end


mutable struct FixedParameters
    LocalSearchRandom::Int64
    LocalSearchBoat::Int64
    PropToRemove::Float64
    WindowSize::Float64
    PushOnlyConstrained::Bool
end


function initializeSol(inst::Instance, allparam)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT = inst
    sol = Sol(inst)
    # solution = [Vector{Tuple{Int64,Int64,Int64}}() for n in 1:N] # (p,b,t) 
    for n in 1:N, (c,p) in enumerate(Pi[n])
        sol.visits[n][c].minT = max(shipsIn[n].sT[c], T[n,c,1])
        sol.visits[n][c].maxT = maxT
    end
    sol.M = generateOccupiedMx(inst, sol.visits)
    #sol = updateTimesInitialization(inst,sol)
    sol.store.parameters = deepcopy(allparam.Proba)
    return sol
end


function initializeParam()

    alphanormal = [0.15,0.2,0.25]
    #alphanormal = [0.2,0.2,0.2,0.2]
    alphaconstrained = [0.2,0.4,0.6]
    alphaconstrained = [0,0,0]
    averagenormal = [1.0,1.0,1.0]
    nbnormal = [0,0,0]
    probanormal = [1/3,1/3,1/3]

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

    AverageCostAll = AverageSolCost(AverageTacticOneBoat,AverageTacticAllBoats,AverageTacticLocalSearch,deepcopy(averagenormal),deepcopy(averagenormal),deepcopy(averagenormal),deepcopy(averagenormal),deepcopy(averagenormal),deepcopy(averagenormal),deepcopy(averagenormal))
    NbAll = NbExperiments(NbTacticOneBoat,NbTacticAllBoats,NbTacticLocalSearch, deepcopy(nbnormal),deepcopy(nbnormal),deepcopy(nbnormal),deepcopy(nbnormal),deepcopy(nbnormal),deepcopy(nbnormal),deepcopy(nbnormal))
    QAll = QParameters(QTacticOneBoat,QTacticAllBoats,QTacticLocalSearch, deepcopy(averagenormal),deepcopy(averagenormal),deepcopy(averagenormal),deepcopy(averagenormal),deepcopy(averagenormal),deepcopy(averagenormal),deepcopy(averagenormal))
    ProbaAll = Probabilities(ProbaTacticOneBoat,ProbaTacticAllBoats,ProbaTacticLocalSearch, deepcopy(probanormal), deepcopy(probanormal), deepcopy(probanormal), deepcopy(probanormal), deepcopy(probanormal), deepcopy(probanormal), deepcopy(probanormal))
    AlphaAll = AlphaParameters(deepcopy(alphanormal),deepcopy(alphanormal),deepcopy(alphanormal),deepcopy(alphanormal),deepcopy(alphanormal),deepcopy(alphanormal),deepcopy(alphaconstrained))
    
    All = AllParameters(ProbaAll,AverageCostAll,NbAll,QAll,AlphaAll)
    
    return All
end


function ChooseParam(allparam::AllParameters, type1, type2, type3)
    @unpack Proba, AverageCost, Nbexp, Q = allparam
    
    if type1=="all"
        tacticoneship = sample(["cost","dist","time"], Weights(Proba.TacticOneBoat))
    end
    if type1=="cost"
        tacticoneship="cost"
    end
    if type1=="dist"
        tacticoneship="dist"
    end
    if type1=="time"
        tacticoneship="time"
    end

    if type2=="all"
        tacticallship = sample(["cost","count","time"], Weights(Proba.TacticAllBoats))
    end
    if type2=="cost"
        tacticallship="cost"
    end
    if type2=="count"
        tacticallship="count"
    end
    if type2=="time"
        tacticallship="time"
    end

    if type3=="all"
        tacticlocal = sample(["random","boat"], Weights(Proba.TacticLocalSearch))
    end
    if type3=="random"
        tacticlocal="random"
    end
    if type3=="boat"
        tacticlocal="boat"
    end
    
    indexone=0
   
    if tacticoneship=="cost"
        indexone = sample(1:length(allparam.Nbexp.CostOneShip), Weights(Proba.CostOneShip))
    end
    if tacticoneship=="dist"
        indexone = sample(1:length(allparam.Nbexp.CostOneShip), Weights(Proba.DistOneShip))
    end
    if tacticoneship=="time"
        indexone = sample(1:length(allparam.Nbexp.CostOneShip), Weights(Proba.TimeOneShip))
    end

    indexall=0
    if tacticallship=="cost"
        indexall = sample(1:length(allparam.Nbexp.CostOneShip), Weights(Proba.CostAllShip))
    end
    if tacticallship=="count"
        indexall = sample(1:length(allparam.Nbexp.CostOneShip), Weights(Proba.CountAllShip))
    end
    if tacticallship=="time"
        indexall = sample(1:length(allparam.Nbexp.CostOneShip), Weights(Proba.TimeAllShip))
    end

    indexconstraint = sample(1:length(allparam.Nbexp.CostOneShip), Weights(Proba.RateConstrained))
    return ChosenParameters(tacticoneship,tacticallship,tacticlocal, indexone, indexall, indexconstraint)
end


function UpdateParameters(paramchosen::ChosenParameters, allparam::AllParameters, cost, new_cost)
    newparam=deepcopy(allparam)
    index=0
    if paramchosen.TacticOneBoat=="cost"
        index=1
    end
    if paramchosen.TacticOneBoat=="dist"
        index=2
    end
    if paramchosen.TacticOneBoat=="time"
        index=3
    end
    nbexp = allparam.Nbexp.TacticOneBoat[index]
    newparam.AverageCost.TacticOneBoat[index] = deepcopy((allparam.AverageCost.TacticOneBoat[index]*nbexp+new_cost)/(nbexp+1))
    newparam.Nbexp.TacticOneBoat[index]=deepcopy(nbexp+1)
    for i in 1:3
        newav = deepcopy(newparam.AverageCost.TacticOneBoat[i])
        newparam.Q.TacticOneBoat[i]=deepcopy(cost/newav)
    end
    for i in 1:3
        newq=deepcopy(newparam.Q.TacticOneBoat[i])
        newparam.Proba.TacticOneBoat[i] = deepcopy(newq/sum(newparam.Q.TacticOneBoat))
    end

    index=0
    if paramchosen.TacticAllBoats=="cost"
        index=1
    end
    if paramchosen.TacticAllBoats=="count"
        index=2
    end
    if paramchosen.TacticAllBoats=="time"
        index=3
    end
    nbexp = deepcopy(allparam.Nbexp.TacticAllBoats[index])
    newparam.AverageCost.TacticAllBoats[index] = deepcopy((allparam.AverageCost.TacticAllBoats[index]*nbexp+new_cost)/(nbexp+1))
    newparam.Nbexp.TacticAllBoats[index]=deepcopy(nbexp+1)
    for i in 1:3
        newav = deepcopy(newparam.AverageCost.TacticAllBoats[i])
        newparam.Q.TacticAllBoats[i]=deepcopy(cost/newav)
    end
    for i in 1:3
        newq=deepcopy(newparam.Q.TacticAllBoats[i])
        newparam.Proba.TacticAllBoats[i] = deepcopy(newq/sum(newparam.Q.TacticAllBoats))
    end

    if paramchosen.TacticOneBoat=="cost"
        nbexp=deepcopy(allparam.Nbexp.CostOneShip[paramchosen.IndexOneShip])
        newparam.AverageCost.CostOneShip[paramchosen.IndexOneShip] = deepcopy((allparam.AverageCost.CostOneShip[paramchosen.IndexOneShip]*nbexp+new_cost)/(nbexp+1))
        newparam.Nbexp.CostOneShip[paramchosen.IndexOneShip]=deepcopy(nbexp+1)
        for i in 1:length(allparam.Nbexp.CostOneShip)
            newav = deepcopy(newparam.AverageCost.CostOneShip[i])
            newparam.Q.CostOneShip[i]=deepcopy(cost/newav)
        end
        for i in 1:length(allparam.Nbexp.CostOneShip)
            newq=deepcopy(newparam.Q.CostOneShip[i])
            newparam.Proba.CostOneShip[i] = deepcopy(newq/sum(newparam.Q.CostOneShip))
        end
    end

    if paramchosen.TacticOneBoat=="dist"
        nbexp=deepcopy(allparam.Nbexp.DistOneShip[paramchosen.IndexOneShip])
        newparam.AverageCost.DistOneShip[paramchosen.IndexOneShip] = deepcopy((allparam.AverageCost.DistOneShip[paramchosen.IndexOneShip]*nbexp+new_cost)/(nbexp+1))
        newparam.Nbexp.DistOneShip[paramchosen.IndexOneShip]=deepcopy(nbexp+1)
        for i in 1:length(allparam.Nbexp.CostOneShip)
            newav = deepcopy(newparam.AverageCost.DistOneShip[i])
            newparam.Q.DistOneShip[i]=deepcopy(cost/newav)
        end
        for i in 1:length(allparam.Nbexp.CostOneShip)
            newq=deepcopy(newparam.Q.DistOneShip[i])
            newparam.Proba.DistOneShip[i] = deepcopy(newq/sum(newparam.Q.DistOneShip))
        end
    end

    if paramchosen.TacticOneBoat=="time"
        nbexp=deepcopy(allparam.Nbexp.TimeOneShip[paramchosen.IndexOneShip])
        newparam.AverageCost.TimeOneShip[paramchosen.IndexOneShip] = deepcopy((allparam.AverageCost.TimeOneShip[paramchosen.IndexOneShip]*nbexp+new_cost)/(nbexp+1))
        newparam.Nbexp.TimeOneShip[paramchosen.IndexOneShip]=deepcopy(nbexp+1)
        for i in 1:length(allparam.Nbexp.CostOneShip)
            newav = deepcopy(newparam.AverageCost.TimeOneShip[i])
            newparam.Q.TimeOneShip[i]=deepcopy(cost/newav)
        end
        for i in 1:length(allparam.Nbexp.CostOneShip)
            newq=deepcopy(newparam.Q.TimeOneShip[i])
            newparam.Proba.TimeOneShip[i] = deepcopy(newq/sum(newparam.Q.TimeOneShip))
        end
    end


    if paramchosen.TacticAllBoats=="cost"
        nbexp=allparam.Nbexp.CostAllShip[paramchosen.IndexAllShip]
        newparam.AverageCost.CostAllShip[paramchosen.IndexAllShip] = (allparam.AverageCost.CostAllShip[paramchosen.IndexAllShip]*nbexp+new_cost)/(nbexp+1)
        newparam.Nbexp.CostAllShip[paramchosen.IndexAllShip]=nbexp+1
        for i in 1:length(allparam.Nbexp.CostOneShip)
            newav = deepcopy(newparam.AverageCost.CostAllShip[i])
            newparam.Q.CostAllShip[i]=cost/newav
        end
        for i in 1:length(allparam.Nbexp.CostOneShip)
            newq=deepcopy(newparam.Q.CostAllShip[i])
            newparam.Proba.CostAllShip[i] = newq/sum(newparam.Q.CostAllShip)
        end
    end

    if paramchosen.TacticAllBoats=="count"
        nbexp=deepcopy(allparam.Nbexp.CountAllShip[paramchosen.IndexAllShip])
        newparam.AverageCost.CountAllShip[paramchosen.IndexAllShip] = deepcopy((allparam.AverageCost.CountAllShip[paramchosen.IndexAllShip]*nbexp+new_cost)/(nbexp+1))
        newparam.Nbexp.CountAllShip[paramchosen.IndexAllShip]=deepcopy(nbexp+1)
        for i in 1:length(allparam.Nbexp.CostOneShip)
            newav = deepcopy(newparam.AverageCost.CountAllShip[i])
            newparam.Q.CountAllShip[i]=deepcopy(cost/newav)
        end
        for i in 1:length(allparam.Nbexp.CostOneShip)
            newq=deepcopy(newparam.Q.CountAllShip[i])
            newparam.Proba.CountAllShip[i] = deepcopy(newq/sum(newparam.Q.CountAllShip))
        end
    end

    if paramchosen.TacticAllBoats=="time"
        nbexp=deepcopy(allparam.Nbexp.TimeAllShip[paramchosen.IndexAllShip])
        newparam.AverageCost.TimeAllShip[paramchosen.IndexAllShip] = deepcopy((allparam.AverageCost.TimeAllShip[paramchosen.IndexAllShip]*nbexp+new_cost)/(nbexp+1))
        newparam.Nbexp.TimeAllShip[paramchosen.IndexAllShip]=deepcopy(nbexp+1)
        for i in 1:length(allparam.Nbexp.CostOneShip)
            newav = deepcopy(newparam.AverageCost.TimeAllShip[i])
            newparam.Q.TimeAllShip[i]=deepcopy(cost/newav)
        end
        for i in 1:length(allparam.Nbexp.CostOneShip)
            newq=deepcopy(newparam.Q.TimeAllShip[i])
            newparam.Proba.TimeAllShip[i] = newq/sum(newparam.Q.TimeAllShip)
        end
    end

    nbexp=deepcopy(allparam.Nbexp.RateConstrained[paramchosen.IndexRateConstrained])
    newparam.AverageCost.RateConstrained[paramchosen.IndexRateConstrained] = deepcopy((allparam.AverageCost.RateConstrained[paramchosen.IndexRateConstrained]*nbexp+new_cost)/(nbexp+1))
    newparam.Nbexp.RateConstrained[paramchosen.IndexRateConstrained]=deepcopy(nbexp+1)
    for i in 1:length(allparam.Nbexp.CostOneShip)
        newav = deepcopy(newparam.AverageCost.RateConstrained[i])
        newparam.Q.RateConstrained[i]=deepcopy(cost/newav)
    end
    for i in 1:length(allparam.Nbexp.CostOneShip)
        newq=deepcopy(newparam.Q.RateConstrained[i])
        newparam.Proba.RateConstrained[i] = deepcopy(newq/sum(newparam.Q.RateConstrained))
    end

    index=0
    if paramchosen.TacticLocalSearch=="random"
        index=1
    end
    if paramchosen.TacticLocalSearch=="boat"
        index=2
    end
    nbexp = deepcopy(allparam.Nbexp.TacticLocalSearch[index])
    newparam.AverageCost.TacticLocalSearch[index] = deepcopy((allparam.AverageCost.TacticLocalSearch[index]*nbexp+new_cost)/(nbexp+1))
    newparam.Nbexp.TacticLocalSearch[index]=deepcopy(nbexp+1)
    for i in 1:2
        newav = deepcopy(newparam.AverageCost.TacticLocalSearch[i])
        newparam.Q.TacticLocalSearch[i]=deepcopy(cost/newav)
    end
    for i in 1:2
        newq=deepcopy(newparam.Q.TacticLocalSearch[i])
        newparam.Proba.TacticLocalSearch[i] = deepcopy(newq/sum(newparam.Q.TacticLocalSearch))
    end
    return newparam
end




