using UnPack 
using Statistics
import XLSX
using CSV, Tables
using DataFrames
using Random
include("toolsMatrixTimes.jl")
include("check_solution.jl")


mutable struct AlphaParameters
    DistOneShip::Vector{Float64}
    CostOneShip::Vector{Float64}
    TimeOneShip::Vector{Float64}

    TimeAllShip::Vector{Float64}
    DistAllShip::Vector{Float64}
    CostAllShip::Vector{Float64}
    
    ReversedTimeAllShip::Vector{Float64}
    ReversedDistAllShip::Vector{Float64}
    ReversedCostAllShip::Vector{Float64}

    RateConstrained::Vector{Float64}
    
    PropToRemove::Vector{Float64}
end

mutable struct AverageSolCost
    TacticOneBoat::Vector{Float64}
    Reversed::Vector{Float64}
    TacticAllBoats::Vector{Float64}
    ReversedTacticAllBoats::Vector{Float64}
    TacticLocalSearch::Vector{Float64}

    TimeOneShip::Vector{Float64}
    DistOneShip::Vector{Float64}
    CostOneShip::Vector{Float64}

    TimeAllShip::Vector{Float64}
    DistAllShip::Vector{Float64}
    CostAllShip::Vector{Float64}

    ReversedTimeAllShip::Vector{Float64}
    ReversedDistAllShip::Vector{Float64}
    ReversedCostAllShip::Vector{Float64}

    RateConstrained::Vector{Float64}

    PropToRemove::Vector{Float64}

end

mutable struct NbExperiments

    TacticOneBoat::Vector{Float64}
    Reversed::Vector{Float64}
    TacticAllBoats::Vector{Float64}
    ReversedTacticAllBoats::Vector{Float64}
    TacticLocalSearch::Vector{Float64}

    TimeOneShip::Vector{Float64}
    DistOneShip::Vector{Float64}
    CostOneShip::Vector{Float64}

    TimeAllShip::Vector{Float64}
    DistAllShip::Vector{Float64}
    CostAllShip::Vector{Float64}

    ReversedTimeAllShip::Vector{Float64}
    ReversedDistAllShip::Vector{Float64}
    ReversedCostAllShip::Vector{Float64}

    RateConstrained::Vector{Float64}

    PropToRemove::Vector{Float64}
end


mutable struct QParameters
    TacticOneBoat::Vector{Float64}
    Reversed::Vector{Float64}
    TacticAllBoats::Vector{Float64}
    ReversedTacticAllBoats::Vector{Float64}
    TacticLocalSearch::Vector{Float64}

    TimeOneShip::Vector{Float64}
    DistOneShip::Vector{Float64}
    CostOneShip::Vector{Float64}

    TimeAllShip::Vector{Float64}
    DistAllShip::Vector{Float64}
    CostAllShip::Vector{Float64}

    ReversedTimeAllShip::Vector{Float64}
    ReversedDistAllShip::Vector{Float64}
    ReversedCostAllShip::Vector{Float64}

    RateConstrained::Vector{Float64}

    PropToRemove::Vector{Float64}
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
    Reversed::String
    TacticAllBoats::String
    ReversedTacticAllBoats::String
    TacticLocalSearch::String
    IndexOneShip::Int64
    IndexAllShip::Int64
    IndexReversedAllShip::Int64
    IndexRateConstrained::Int64
    IndexPropToRemove::Int64
end


mutable struct FixedParameters
    OneBoat::String
    Alphaoneboat::Vector{Float64}
    Reversed::String
    AllBoat::String
    Alphaallboat::Vector{Float64}
    ReversedAllBoat::String
    Alphareversedallboat::Vector{Float64}
    LocalSearch::String
    LocalSearchOne::String
    LocalSearchAll::String
    LocalSearchRandom::Float64
    LocalSearchBoat::Float64
    Alpharateconstrained::Vector{Float64}
    Alphapropremove::Vector{Float64}
    WindowSize::Float64
    PushOnlyConstrained::Bool
    lookforconstrained::Bool
    restartParams::Int64
    GreedymaxNoImprove::Int64
    Until::Int64
    FocusRemoveUntil::Int64
    NbFocusRemove::Int64
    RateImproveReconstructOrPath::Float64
    windowLocalSearch::Float64
    pathRelinking::String
    MaxTimeRelinking::Float64
    LengthElite::Int64
    RemovePathRelinking::Float64
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


function initializeParam(fixedparam::FixedParameters)

    alphaoneship = fixedparam.Alphaoneboat
    alphallship = fixedparam.Alphaallboat
    alphareversedallship = fixedparam.Alphareversedallboat
    alpharateconstrained = fixedparam.Alpharateconstrained
    alphaprobtoremove = fixedparam.Alphapropremove

    averageoneship = fill(1,length(alphaoneship))
    averageallship = fill(1,length(alphallship))
    averageallshipreversed = fill(1,length(alphareversedallship))
    averagerateconstrained = fill(1,length(alpharateconstrained))
    averageprobtoremove = fill(1,length(alphaprobtoremove))

    nboneship = fill(0,length(alphaoneship))
    nballship = fill(0,length(alphallship))
    nbreversedallship = fill(0,length(alphareversedallship))
    nbconstrained = fill(0,length(alpharateconstrained))
    nbprobtoremove = fill(0,length(alphaprobtoremove))

    qoneship = fill(1,length(alphaoneship))
    qallship = fill(1,length(alphallship))
    qreversedallship = fill(1,length(alphareversedallship))
    qconstrained = fill(1,length(alpharateconstrained))
    qprobtoremove = fill(1,length(alphaprobtoremove))

    proboneship = fill(1/length(alphaoneship),length(alphaoneship))
    proballship = fill(1/length(alphallship),length(alphallship))
    probreversedallship = fill(1/length(alphareversedallship),length(alphareversedallship))
    probconstrained = fill(1/length(alpharateconstrained),length(alpharateconstrained))
    probprobtoremove = fill(1/length(alphaprobtoremove),length(alphaprobtoremove))


    AverageTacticOneBoat = [1.0,1.0,1.0]
    AverageTacticAllBoats = [1.0,1.0,1.0]
    AverageTacticLocalSearch = [1.0,1.0,1.0]
    AverageReversed = [1.0,1.0]
    AverageReversedTacticAllBoats = [1.0,1.0,1.0]

    NbTacticOneBoat = [0,0,0]
    NbTacticAllBoats = [0,0,0]
    NbTacticLocalSearch = [0,0,0]
    NbReversed = [0,0]
    NbReversedTacticAllBoats = [0,0,0]

    QTacticOneBoat = [1.0,1.0,1.0]
    QTacticAllBoats = [1.0,1.0,1.0]
    QTacticLocalSearch = [1.0,1.0,1.0]
    QReversed = [1.0,1.0]
    QReversedTacticAllBoats = [1.0,1.0,1.0]

    ProbaTacticOneBoat = [1/3,1/3,1/3]
    ProbaTacticAllBoats = [1/3,1/3,1/3]
    ProbaTacticLocalSearch = [1/3,1/3,1/3]
    ProbaReversed = [1/2,1/2]
    ProbaReversedTacticAllBoats = [1/3,1/3,1/3]

    AverageCostAll = AverageSolCost(AverageTacticOneBoat,AverageReversed, AverageTacticAllBoats,AverageReversedTacticAllBoats, AverageTacticLocalSearch,
    deepcopy(averageoneship),deepcopy(averageoneship),deepcopy(averageoneship),
    deepcopy(averageallship),deepcopy(averageallship),deepcopy(averageallship),
    deepcopy(averageallshipreversed),deepcopy(averageallshipreversed),deepcopy(averageallshipreversed),
    deepcopy(averagerateconstrained),deepcopy(averageprobtoremove))
    NbAll = NbExperiments(NbTacticOneBoat,NbReversed, NbTacticAllBoats,NbReversedTacticAllBoats,NbTacticLocalSearch, 
    deepcopy(nboneship),deepcopy(nboneship),deepcopy(nboneship),
    deepcopy(nballship),deepcopy(nballship),deepcopy(nballship),
    deepcopy(nbreversedallship),deepcopy(nbreversedallship),deepcopy(nbreversedallship),
    deepcopy(nbconstrained),deepcopy(nbprobtoremove))
    QAll = QParameters(QTacticOneBoat,QReversed, QTacticAllBoats,QReversedTacticAllBoats, QTacticLocalSearch,
    deepcopy(qoneship),deepcopy(qoneship),deepcopy(qoneship),
    deepcopy(qallship),deepcopy(qallship),deepcopy(qallship),
    deepcopy(qreversedallship),deepcopy(qreversedallship),deepcopy(qreversedallship),
    deepcopy(qconstrained),deepcopy(qprobtoremove))
    ProbaAll = Probabilities(ProbaTacticOneBoat,ProbaReversed, ProbaTacticAllBoats,ProbaReversedTacticAllBoats, ProbaTacticLocalSearch,
    deepcopy(proboneship),deepcopy(proboneship),deepcopy(proboneship),
    deepcopy(proballship),deepcopy(proballship),deepcopy(proballship),
    deepcopy(probreversedallship),deepcopy(probreversedallship),deepcopy(probreversedallship),
    deepcopy(probconstrained),deepcopy(probprobtoremove))
    AlphaAll = AlphaParameters(deepcopy(alphaoneship),deepcopy(alphaoneship),deepcopy(alphaoneship),
    deepcopy(alphallship),deepcopy(alphallship),deepcopy(alphallship),
    deepcopy(alphareversedallship),deepcopy(alphareversedallship),deepcopy(alphareversedallship),
    deepcopy(alpharateconstrained),deepcopy(alphaprobtoremove))
    
    All = AllParameters(ProbaAll,AverageCostAll,NbAll,QAll,AlphaAll)
    
    return All
end


function ChooseAllParam(allparam::AllParameters, paramfixed::FixedParameters)
    @unpack Proba, AverageCost, Nbexp, Q = allparam
    
    if paramfixed.OneBoat=="all"
        tacticoneship = sample(["cost","dist","time"], Weights(Proba.TacticOneBoat))
    end
    if paramfixed.OneBoat=="cost"
        tacticoneship="cost"
    end
    if paramfixed.OneBoat=="dist"
        tacticoneship="dist"
    end
    if paramfixed.OneBoat=="time"
        tacticoneship="time"
    end
    
    if paramfixed.Reversed=="all"
        reversedtactic = sample(["yes","no"], Weights(Proba.Reversed))
    end
    if paramfixed.Reversed=="yes"
        reversedtactic="yes"
    end
    if paramfixed.Reversed=="no"
        reversedtactic="no"
    end

    if paramfixed.AllBoat=="all"
        tacticallship = sample(["cost","dist","time"], Weights(Proba.TacticAllBoats))
    end
    if paramfixed.AllBoat=="cost"
        tacticallship="cost"
    end
    if paramfixed.AllBoat=="dist"
        tacticallship="dist"
    end
    if paramfixed.AllBoat=="time"
        tacticallship="time"
    end

    if paramfixed.ReversedAllBoat=="all"
        tacticreversedallship = sample(["cost","dist","time"], Weights(Proba.ReversedTacticAllBoats))
    end
    if paramfixed.ReversedAllBoat=="cost"
        tacticreversedallship="cost"
    end
    if paramfixed.ReversedAllBoat=="dist"
        tacticreversedallship="dist"
    end
    if paramfixed.ReversedAllBoat=="time"
        tacticreversedallship="time"
    end

    if paramfixed.LocalSearch=="all"
        tacticlocal = sample(["cost","dist","time"], Weights(Proba.TacticLocalSearch))
    end
    if paramfixed.LocalSearch=="cost"
        tacticlocal="cost"
    end
    if paramfixed.LocalSearch=="dist"
        tacticlocal="dist"
    end
    if paramfixed.LocalSearch=="time"
        tacticlocal="time"
    end
    
    indexone=0
   
    if tacticoneship=="cost"
        indexone = sample(1:length(allparam.Nbexp.CostOneShip), Weights(Proba.CostOneShip))
    end
    if tacticoneship=="dist"
        indexone = sample(1:length(allparam.Nbexp.DistOneShip), Weights(Proba.DistOneShip))
    end
    if tacticoneship=="time"
        indexone = sample(1:length(allparam.Nbexp.TimeOneShip), Weights(Proba.TimeOneShip))
    end
    
    #print('\n')
    #print(1:length(allparam.Nbexp.TimeAllShip))
    indexall=0
    if tacticallship=="cost"
        indexall = sample(1:length(allparam.Nbexp.CostAllShip), Weights(Proba.CostAllShip))
    end
    if tacticallship=="dist"
        indexall = sample(1:length(allparam.Nbexp.DistAllShip), Weights(Proba.DistAllShip))
    end
    if tacticallship=="time"
        indexall = sample(1:length(allparam.Nbexp.TimeAllShip), Weights(Proba.TimeAllShip))
    end
    #print('\n')
    #print(indexall)
    #print("<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>")
    indexreversedall=0

    if tacticreversedallship=="cost"
        indexreversedall = sample(1:length(allparam.Nbexp.ReversedCostAllShip), Weights(Proba.ReversedCostAllShip))
    end
    if tacticreversedallship=="dist"
        indexreversedall = sample(1:length(allparam.Nbexp.ReversedDistAllShip), Weights(Proba.ReversedDistAllShip))
    end
    if tacticreversedallship=="time"
        indexreversedall = sample(1:length(allparam.Nbexp.ReversedTimeAllShip), Weights(Proba.ReversedTimeAllShip))
    end

    indexconstraint = sample(1:length(allparam.Nbexp.RateConstrained), Weights(Proba.RateConstrained))
    indexprobtoremove = sample(1:length(allparam.Nbexp.PropToRemove), Weights(Proba.PropToRemove))
    print('\n')
    print(tacticoneship)
    return ChosenParameters(tacticoneship, reversedtactic, tacticallship, tacticreversedallship, tacticlocal, indexone, indexall, indexreversedall, indexconstraint, indexprobtoremove)
end


function UpdateAllParameters(paramchosen::ChosenParameters, allparam::AllParameters, cost, new_cost, constrained)
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
    if paramchosen.Reversed=="yes"
        index=1
    end
    if paramchosen.Reversed=="no"
        index=2
    end
    nbexp = allparam.Nbexp.Reversed[index]
    newparam.AverageCost.Reversed[index] = deepcopy((allparam.AverageCost.Reversed[index]*nbexp+new_cost)/(nbexp+1))
    newparam.Nbexp.Reversed[index]=deepcopy(nbexp+1)
    for i in 1:2
        newav = deepcopy(newparam.AverageCost.Reversed[i])
        newparam.Q.Reversed[i]=deepcopy(cost/newav)
    end
    for i in 1:2
        newq=deepcopy(newparam.Q.Reversed[i])
        newparam.Proba.Reversed[i] = deepcopy(newq/sum(newparam.Q.Reversed))
    end


    if paramchosen.Reversed=="no"
        index=0
        if paramchosen.TacticAllBoats=="cost"
            index=1
        end
        if paramchosen.TacticAllBoats=="dist"
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
    else
        index=0
        if paramchosen.ReversedTacticAllBoats=="cost"
            index=1
        end
        if paramchosen.ReversedTacticAllBoats=="dist"
            index=2
        end
        if paramchosen.ReversedTacticAllBoats=="time"
            index=3
        end
        nbexp = deepcopy(allparam.Nbexp.ReversedTacticAllBoats[index])
        newparam.AverageCost.ReversedTacticAllBoats[index] = deepcopy((allparam.AverageCost.ReversedTacticAllBoats[index]*nbexp+new_cost)/(nbexp+1))
        newparam.Nbexp.ReversedTacticAllBoats[index]=deepcopy(nbexp+1)
        for i in 1:3
            newav = deepcopy(newparam.AverageCost.ReversedTacticAllBoats[i])
            newparam.Q.ReversedTacticAllBoats[i]=deepcopy(cost/newav)
        end
        for i in 1:3
            newq=deepcopy(newparam.Q.ReversedTacticAllBoats[i])
            newparam.Proba.ReversedTacticAllBoats[i] = deepcopy(newq/sum(newparam.Q.ReversedTacticAllBoats))
        end
    end

    index=0
    if paramchosen.TacticLocalSearch=="cost"
        index=1
    end
    if paramchosen.TacticLocalSearch=="dist"
        index=2
    end
    if paramchosen.TacticLocalSearch=="time"
        index=3
    end
    nbexp = deepcopy(allparam.Nbexp.TacticLocalSearch[index])
    newparam.AverageCost.TacticLocalSearch[index] = deepcopy((allparam.AverageCost.TacticLocalSearch[index]*nbexp+new_cost)/(nbexp+1))
    newparam.Nbexp.TacticLocalSearch[index]=deepcopy(nbexp+1)
    for i in 1:3
        newav = deepcopy(newparam.AverageCost.TacticLocalSearch[i])
        newparam.Q.TacticLocalSearch[i]=deepcopy(cost/newav)
    end
    for i in 1:3
        newq=deepcopy(newparam.Q.TacticLocalSearch[i])
        newparam.Proba.TacticLocalSearch[i] = deepcopy(newq/sum(newparam.Q.TacticLocalSearch))
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


    if paramchosen.Reversed=="no"
        if paramchosen.TacticAllBoats=="cost"
            nbexp=allparam.Nbexp.CostAllShip[paramchosen.IndexAllShip]
            newparam.AverageCost.CostAllShip[paramchosen.IndexAllShip] = (allparam.AverageCost.CostAllShip[paramchosen.IndexAllShip]*nbexp+new_cost)/(nbexp+1)
            newparam.Nbexp.CostAllShip[paramchosen.IndexAllShip]=nbexp+1
            for i in 1:length(allparam.Nbexp.CostAllShip)
                newav = deepcopy(newparam.AverageCost.CostAllShip[i])
                newparam.Q.CostAllShip[i]=cost/newav
            end
            for i in 1:length(allparam.Nbexp.CostAllShip)
                newq=deepcopy(newparam.Q.CostAllShip[i])
                newparam.Proba.CostAllShip[i] = newq/sum(newparam.Q.CostAllShip)
            end
        end

        if paramchosen.TacticAllBoats=="dist"
            nbexp=deepcopy(allparam.Nbexp.DistAllShip[paramchosen.IndexAllShip])
            newparam.AverageCost.DistAllShip[paramchosen.IndexAllShip] = deepcopy((allparam.AverageCost.DistAllShip[paramchosen.IndexAllShip]*nbexp+new_cost)/(nbexp+1))
            newparam.Nbexp.DistAllShip[paramchosen.IndexAllShip]=deepcopy(nbexp+1)
            for i in 1:length(allparam.Nbexp.CostAllShip)
                newav = deepcopy(newparam.AverageCost.DistAllShip[i])
                newparam.Q.DistAllShip[i]=deepcopy(cost/newav)
            end
            for i in 1:length(allparam.Nbexp.CostAllShip)
                newq=deepcopy(newparam.Q.DistAllShip[i])
                newparam.Proba.DistAllShip[i] = deepcopy(newq/sum(newparam.Q.DistAllShip))
            end
        end

        if paramchosen.TacticAllBoats=="time"
            nbexp=deepcopy(allparam.Nbexp.TimeAllShip[paramchosen.IndexAllShip])
            newparam.AverageCost.TimeAllShip[paramchosen.IndexAllShip] = deepcopy((allparam.AverageCost.TimeAllShip[paramchosen.IndexAllShip]*nbexp+new_cost)/(nbexp+1))
            newparam.Nbexp.TimeAllShip[paramchosen.IndexAllShip]=deepcopy(nbexp+1)
            for i in 1:length(allparam.Nbexp.CostAllShip)
                newav = deepcopy(newparam.AverageCost.TimeAllShip[i])
                newparam.Q.TimeAllShip[i]=deepcopy(cost/newav)
            end
            for i in 1:length(allparam.Nbexp.CostAllShip)
                newq=deepcopy(newparam.Q.TimeAllShip[i])
                newparam.Proba.TimeAllShip[i] = newq/sum(newparam.Q.TimeAllShip)
            end
        end
    else
        if paramchosen.ReversedTacticAllBoats=="cost"
            nbexp=deepcopy(allparam.Nbexp.ReversedCostAllShip[paramchosen.IndexReversedAllShip])
            newparam.AverageCost.ReversedCostAllShip[paramchosen.IndexAllShip] = deepcopy((allparam.AverageCost.ReversedCostAllShip[paramchosen.IndexReversedAllShip]*nbexp+new_cost)/(nbexp+1))
            newparam.Nbexp.ReversedCostAllShip[paramchosen.IndexReversedAllShip]=deepcopy(nbexp+1)
            for i in 1:length(allparam.Nbexp.ReversedCostAllShip)
                newav = deepcopy(newparam.AverageCost.ReversedCostAllShip[i])
                newparam.Q.ReversedCostAllShip[i]=deepcopy(cost/newav)
            end
            for i in 1:length(allparam.Nbexp.ReversedCostAllShip)
                newq=deepcopy(newparam.Q.ReversedCostAllShip[i])
                newparam.Proba.ReversedCostAllShip[i] = deepcopy(newq/sum(newparam.Q.ReversedCostAllShip))
            end
        end
        if paramchosen.ReversedTacticAllBoats=="dist"
            nbexp=deepcopy(allparam.Nbexp.ReversedDistAllShip[paramchosen.IndexReversedAllShip])
            newparam.AverageCost.ReversedDistAllShip[paramchosen.IndexReversedAllShip] = deepcopy((allparam.AverageCost.ReversedDistAllShip[paramchosen.IndexReversedAllShip]*nbexp+new_cost)/(nbexp+1))
            newparam.Nbexp.ReversedDistAllShip[paramchosen.IndexReversedAllShip]=deepcopy(nbexp+1)
            for i in 1:length(allparam.Nbexp.ReversedCostAllShip)
                newav = deepcopy(newparam.AverageCost.ReversedDistAllShip[i])
                newparam.Q.ReversedDistAllShip[i]=deepcopy(cost/newav)
            end
            for i in 1:length(allparam.Nbexp.ReversedCostAllShip)
                newq=deepcopy(newparam.Q.ReversedDistAllShip[i])
                newparam.Proba.ReversedDistAllShip[i] = deepcopy(newq/sum(newparam.Q.ReversedDistAllShip))
            end
        end
        if paramchosen.ReversedTacticAllBoats=="time"
            nbexp=deepcopy(allparam.Nbexp.ReversedTimeAllShip[paramchosen.IndexReversedAllShip])
            newparam.AverageCost.ReversedTimeAllShip[paramchosen.IndexReversedAllShip] = deepcopy((allparam.AverageCost.ReversedTimeAllShip[paramchosen.IndexReversedAllShip]*nbexp+new_cost)/(nbexp+1))
            newparam.Nbexp.ReversedTimeAllShip[paramchosen.IndexReversedAllShip]=deepcopy(nbexp+1)
            for i in 1:length(allparam.Nbexp.ReversedCostAllShip)
                newav = deepcopy(newparam.AverageCost.ReversedTimeAllShip[i])
                newparam.Q.ReversedTimeAllShip[i]=deepcopy(cost/newav)
            end
            for i in 1:length(allparam.Nbexp.ReversedCostAllShip)
                newq=deepcopy(newparam.Q.ReversedTimeAllShip[i])
                newparam.Proba.ReversedTimeAllShip[i] = deepcopy(newq/sum(newparam.Q.ReversedTimeAllShip))
            end
        end
    end

    if constrained
        nbexp=deepcopy(allparam.Nbexp.RateConstrained[paramchosen.IndexRateConstrained])
        newparam.AverageCost.RateConstrained[paramchosen.IndexRateConstrained] = deepcopy((allparam.AverageCost.RateConstrained[paramchosen.IndexRateConstrained]*nbexp+new_cost)/(nbexp+1))
        newparam.Nbexp.RateConstrained[paramchosen.IndexRateConstrained]=deepcopy(nbexp+1)
        for i in 1:length(allparam.Nbexp.RateConstrained)
            newav = deepcopy(newparam.AverageCost.RateConstrained[i])
            newparam.Q.RateConstrained[i]=deepcopy(cost/newav)
        end
        for i in 1:length(allparam.Nbexp.RateConstrained)
            newq=deepcopy(newparam.Q.RateConstrained[i])
            newparam.Proba.RateConstrained[i] = deepcopy(newq/sum(newparam.Q.RateConstrained))
        end
    end

    nbexp=deepcopy(allparam.Nbexp.PropToRemove[paramchosen.IndexPropToRemove])
    newparam.AverageCost.PropToRemove[paramchosen.IndexPropToRemove] = deepcopy((allparam.AverageCost.PropToRemove[paramchosen.IndexPropToRemove]*nbexp+new_cost)/(nbexp+1))
    newparam.Nbexp.PropToRemove[paramchosen.IndexPropToRemove]=deepcopy(nbexp+1)
    for i in 1:length(allparam.Nbexp.PropToRemove)
        newav = deepcopy(newparam.AverageCost.PropToRemove[i])
        newparam.Q.PropToRemove[i]=deepcopy(cost/newav)
    end
    for i in 1:length(allparam.Nbexp.PropToRemove)
        newq=deepcopy(newparam.Q.PropToRemove[i])
        newparam.Proba.PropToRemove[i] = deepcopy(newq/sum(newparam.Q.PropToRemove))
    end

    return newparam
end



function ChooseAfterReconstruct(allparam::AllParameters, paramchosen::ChosenParameters, paramfixed::FixedParameters)
    @unpack Proba, AverageCost, Nbexp, Q = allparam
    new_paramchosen = deepcopy(paramchosen)
    
    if paramfixed.OneBoat=="all"
        tacticoneship = sample(["cost","dist","time"], Weights(Proba.TacticOneBoat))
    end
    if paramfixed.OneBoat=="cost"
        tacticoneship="cost"
    end
    if paramfixed.OneBoat=="dist"
        tacticoneship="dist"
    end
    if paramfixed.OneBoat=="time"
        tacticoneship="time"
    end
    
    if paramfixed.Reversed=="all"
        reversedtactic = sample(["yes","no"], Weights(Proba.Reversed))
    end
    if paramfixed.Reversed=="yes"
        reversedtactic="yes"
    end
    if paramfixed.Reversed=="no"
        reversedtactic="no"
    end

    if paramfixed.AllBoat=="all"
        tacticallship = sample(["cost","dist","time"], Weights(Proba.TacticAllBoats))
    end
    if paramfixed.AllBoat=="cost"
        tacticallship="cost"
    end
    if paramfixed.AllBoat=="dist"
        tacticallship="dist"
    end
    if paramfixed.AllBoat=="time"
        tacticallship="time"
    end

    if paramfixed.ReversedAllBoat=="all"
        tacticreversedallship = sample(["cost","dist","time"], Weights(Proba.ReversedTacticAllBoats))
    end
    if paramfixed.ReversedAllBoat=="cost"
        tacticreversedallship="cost"
    end
    if paramfixed.ReversedAllBoat=="dist"
        tacticreversedallship="dist"
    end
    if paramfixed.ReversedAllBoat=="time"
        tacticreversedallship="time"
    end
    
    indexone=0
   
    if tacticoneship=="cost"
        indexone = sample(1:length(allparam.Nbexp.CostOneShip), Weights(Proba.CostOneShip))
    end
    if tacticoneship=="dist"
        indexone = sample(1:length(allparam.Nbexp.DistOneShip), Weights(Proba.DistOneShip))
    end
    if tacticoneship=="time"
        indexone = sample(1:length(allparam.Nbexp.TimeOneShip), Weights(Proba.TimeOneShip))
    end

    indexall=0
    if tacticallship=="cost"
        indexall = sample(1:length(allparam.Nbexp.CostAllShip), Weights(Proba.CostAllShip))
    end
    if tacticallship=="dist"
        indexall = sample(1:length(allparam.Nbexp.DistAllShip), Weights(Proba.DistAllShip))
    end
    if tacticallship=="time"
        indexall = sample(1:length(allparam.Nbexp.TimeAllShip), Weights(Proba.TimeAllShip))
    end

    indexreversedall=0
    if tacticreversedallship=="cost"
        indexreversedall = sample(1:length(allparam.Nbexp.ReversedCostAllShip), Weights(Proba.ReversedCostAllShip))
    end
    if tacticreversedallship=="dist"
        indexreversedall = sample(1:length(allparam.Nbexp.ReversedCostAllShip), Weights(Proba.ReversedDistAllShip))
    end
    if tacticreversedallship=="time"
        indexreversedall = sample(1:length(allparam.Nbexp.ReversedTimeAllShip), Weights(Proba.ReversedTimeAllShip))
    end

    indexconstraint = sample(1:length(allparam.Nbexp.RateConstrained), Weights(Proba.RateConstrained))
    indexprobtoremove = sample(1:length(allparam.Nbexp.PropToRemove), Weights(Proba.PropToRemove))

    new_paramchosen.TacticOneShip = tacticoneship
    new_paramchosen.Reversed = reversedtactic
    new_paramchosen.TacticAllShip = tacticallship
    new_paramchosen.TacticReversedAllShip = tacticreversedallship
    new_paramchosen.IndexOneShip = indexone
    new_paramchosen.IndexReversedAllShip = indexreversedall
    new_paramchosen.IndexAllShip = indexall
    new_paramchosen.IndexConstraint = indexconstraint
    new_paramchosen.IndexProbToRemove = indexprobtoremove
    
    return new_paramchosen
end



function ChooseAfterHeur(allparam::AllParameters, paramchosen::ChosenParameters, paramfixed::FixedParameters)
    @unpack Proba, AverageCost, Nbexp, Q = allparam
    new_paramchosen = deepcopy(paramchosen)
    
    if paramfixed.OneBoat=="all"
        tacticoneship = sample(["cost","dist","time"], Weights(Proba.TacticOneBoat))
    end
    if paramfixed.OneBoat=="cost"
        tacticoneship="cost"
    end
    if paramfixed.OneBoat=="dist"
        tacticoneship="dist"
    end
    if paramfixed.OneBoat=="time"
        tacticoneship="time"
    end
    
    if paramfixed.Reversed=="all"
        reversedtactic = sample(["yes","no"], Weights(Proba.Reversed))
    end
    if paramfixed.Reversed=="yes"
        reversedtactic="yes"
    end
    if paramfixed.Reversed=="no"
        reversedtactic="no"
    end

    if paramfixed.AllBoat=="all"
        tacticallship = sample(["cost","dist","time"], Weights(Proba.TacticAllBoats))
    end
    if paramfixed.AllBoat=="cost"
        tacticallship="cost"
    end
    if paramfixed.AllBoat=="dist"
        tacticallship="dist"
    end
    if paramfixed.AllBoat=="time"
        tacticallship="time"
    end

    if paramfixed.ReversedAllBoat=="all"
        tacticreversedallship = sample(["cost","dist","time"], Weights(Proba.ReversedTacticAllBoats))
    end
    if paramfixed.ReversedAllBoat=="cost"
        tacticreversedallship="cost"
    end
    if paramfixed.ReversedAllBoat=="dist"
        tacticreversedallship="dist"
    end
    if paramfixed.ReversedAllBoat=="time"
        tacticreversedallship="time"
    end

    indexone=0
   
    if tacticoneship=="cost"
        indexone = sample(1:length(allparam.Nbexp.CostOneShip), Weights(Proba.CostOneShip))
    end
    if tacticoneship=="dist"
        indexone = sample(1:length(allparam.Nbexp.DistOneShip), Weights(Proba.DistOneShip))
    end
    if tacticoneship=="time"
        indexone = sample(1:length(allparam.Nbexp.TimeOneShip), Weights(Proba.TimeOneShip))
    end

    indexall=0
    if tacticallship=="cost"
        indexall = sample(1:length(allparam.Nbexp.CostAllShip), Weights(Proba.CostAllShip))
    end
    if tacticallship=="dist"
        indexall = sample(1:length(allparam.Nbexp.DistAllShip), Weights(Proba.DistAllShip))
    end
    if tacticallship=="time"
        indexall = sample(1:length(allparam.Nbexp.TimeAllShip), Weights(Proba.TimeAllShip))
    end

    indexreversedall=0
    if tacticreversedallship=="cost"
        indexreversedall = sample(1:length(allparam.Nbexp.ReversedCostAllShip), Weights(Proba.ReversedCostAllShip))
    end
    if tacticreversedallship=="dist"
        indexreversedall = sample(1:length(allparam.Nbexp.ReversedCostAllShip), Weights(Proba.ReversedDistAllShip))
    end
    if tacticreversedallship=="time"
        indexreversedall = sample(1:length(allparam.Nbexp.ReversedTimeAllShip), Weights(Proba.ReversedTimeAllShip))
    end

    indexconstraint = sample(1:length(allparam.Nbexp.RateConstrained), Weights(Proba.RateConstrained))

    new_paramchosen.TacticOneShip = tacticoneship
    new_paramchosen.Reversed = reversedtactic
    new_paramchosen.TacticAllShip = tacticallship
    new_paramchosen.TacticReversedAllShip = tacticreversedallship
    new_paramchosen.IndexOneShip = indexone
    new_paramchosen.IndexReversedAllShip = indexreversedall
    new_paramchosen.IndexAllShip = indexall
    new_paramchosen.IndexConstraint = indexconstraint
    
    return new_paramchosen
end


function ChooseAfterLocal(allparam::AllParameters, paramchosen::ChosenParameters, paramfixed::FixedParameters)
    @unpack Proba, AverageCost, Nbexp, Q = allparam
    new_paramchosen = deepcopy(paramchosen)
    
    if paramfixed.LocalSearch=="all"
        tacticlocal = sample(["cost","dist","time"], Weights(Proba.TacticLocalSearch))
    end
    if paramfixed.LocalSearch=="cost"
        tacticlocal="cost"
    end
    if paramfixed.LocalSearch=="dist"
        tacticlocal="dist"
    end
    if paramfixed.LocalSearch=="time"
        tacticlocal="time"
    end

    new_paramchosen.TacticLocalSearch = tacticlocal
    
    return new_paramchosen
end


function UpdateAfterReconstructParameters(paramchosen::ChosenParameters, allparam::AllParameters, cost, new_cost, constrained)
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
    if paramchosen.Reversed=="yes"
        index=1
    end
    if paramchosen.Reversed=="no"
        index=2
    end
    nbexp = allparam.Nbexp.Reversed[index]
    newparam.AverageCost.Reversed[index] = deepcopy((allparam.AverageCost.Reversed[index]*nbexp+new_cost)/(nbexp+1))
    newparam.Nbexp.Reversed[index]=deepcopy(nbexp+1)
    for i in 1:2
        newav = deepcopy(newparam.AverageCost.Reversed[i])
        newparam.Q.Reversed[i]=deepcopy(cost/newav)
    end
    for i in 1:2
        newq=deepcopy(newparam.Q.Reversed[i])
        newparam.Proba.Reversed[i] = deepcopy(newq/sum(newparam.Q.Reversed))
    end


    if paramchosen.Reversed=="no"
        index=0
        if paramchosen.TacticAllBoats=="cost"
            index=1
        end
        if paramchosen.TacticAllBoats=="dist"
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
    else
        index=0
        if paramchosen.ReversedTacticAllBoats=="cost"
            index=1
        end
        if paramchosen.ReversedTacticAllBoats=="dist"
            index=2
        end
        if paramchosen.ReversedTacticAllBoats=="time"
            index=3
        end
        nbexp = deepcopy(allparam.Nbexp.ReversedTacticAllBoats[index])
        newparam.AverageCost.ReversedTacticAllBoats[index] = deepcopy((allparam.AverageCost.ReversedTacticAllBoats[index]*nbexp+new_cost)/(nbexp+1))
        newparam.Nbexp.ReversedTacticAllBoats[index]=deepcopy(nbexp+1)
        for i in 1:3
            newav = deepcopy(newparam.AverageCost.ReversedTacticAllBoats[i])
            newparam.Q.ReversedTacticAllBoats[i]=deepcopy(cost/newav)
        end
        for i in 1:3
            newq=deepcopy(newparam.Q.ReversedTacticAllBoats[i])
            newparam.Proba.ReversedTacticAllBoats[i] = deepcopy(newq/sum(newparam.Q.ReversedTacticAllBoats))
        end
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


    if paramchosen.Reversed=="no"
        if paramchosen.TacticAllBoats=="cost"
            nbexp=allparam.Nbexp.CostAllShip[paramchosen.IndexAllShip]
            newparam.AverageCost.CostAllShip[paramchosen.IndexAllShip] = (allparam.AverageCost.CostAllShip[paramchosen.IndexAllShip]*nbexp+new_cost)/(nbexp+1)
            newparam.Nbexp.CostAllShip[paramchosen.IndexAllShip]=nbexp+1
            for i in 1:length(allparam.Nbexp.CostAllShip)
                newav = deepcopy(newparam.AverageCost.CostAllShip[i])
                newparam.Q.CostAllShip[i]=cost/newav
            end
            for i in 1:length(allparam.Nbexp.CostAllShip)
                newq=deepcopy(newparam.Q.CostAllShip[i])
                newparam.Proba.CostAllShip[i] = newq/sum(newparam.Q.CostAllShip)
            end
        end

        if paramchosen.TacticAllBoats=="dist"
            nbexp=deepcopy(allparam.Nbexp.DistAllShip[paramchosen.IndexAllShip])
            newparam.AverageCost.DistAllShip[paramchosen.IndexAllShip] = deepcopy((allparam.AverageCost.DistAllShip[paramchosen.IndexAllShip]*nbexp+new_cost)/(nbexp+1))
            newparam.Nbexp.DistAllShip[paramchosen.IndexAllShip]=deepcopy(nbexp+1)
            for i in 1:length(allparam.Nbexp.CostAllShip)
                newav = deepcopy(newparam.AverageCost.DistAllShip[i])
                newparam.Q.DistAllShip[i]=deepcopy(cost/newav)
            end
            for i in 1:length(allparam.Nbexp.CostAllShip)
                newq=deepcopy(newparam.Q.DistAllShip[i])
                newparam.Proba.DistAllShip[i] = deepcopy(newq/sum(newparam.Q.DistAllShip))
            end
        end

        if paramchosen.TacticAllBoats=="time"
            nbexp=deepcopy(allparam.Nbexp.TimeAllShip[paramchosen.IndexAllShip])
            newparam.AverageCost.TimeAllShip[paramchosen.IndexAllShip] = deepcopy((allparam.AverageCost.TimeAllShip[paramchosen.IndexAllShip]*nbexp+new_cost)/(nbexp+1))
            newparam.Nbexp.TimeAllShip[paramchosen.IndexAllShip]=deepcopy(nbexp+1)
            for i in 1:length(allparam.Nbexp.CostAllShip)
                newav = deepcopy(newparam.AverageCost.TimeAllShip[i])
                newparam.Q.TimeAllShip[i]=deepcopy(cost/newav)
            end
            for i in 1:length(allparam.Nbexp.CostAllShip)
                newq=deepcopy(newparam.Q.TimeAllShip[i])
                newparam.Proba.TimeAllShip[i] = newq/sum(newparam.Q.TimeAllShip)
            end
        end
    else
        if paramchosen.ReversedTacticAllBoats=="cost"
            nbexp=deepcopy(allparam.Nbexp.ReversedCostAllShip[paramchosen.IndexReversedAllShip])
            newparam.AverageCost.ReversedCostAllShip[paramchosen.IndexAllShip] = deepcopy((allparam.AverageCost.ReversedCostAllShip[paramchosen.IndexReversedAllShip]*nbexp+new_cost)/(nbexp+1))
            newparam.Nbexp.ReversedCostAllShip[paramchosen.IndexReversedAllShip]=deepcopy(nbexp+1)
            for i in 1:length(allparam.Nbexp.ReversedCostAllShip)
                newav = deepcopy(newparam.AverageCost.ReversedCostAllShip[i])
                newparam.Q.ReversedCostAllShip[i]=deepcopy(cost/newav)
            end
            for i in 1:length(allparam.Nbexp.ReversedCostAllShip)
                newq=deepcopy(newparam.Q.ReversedCostAllShip[i])
                newparam.Proba.ReversedCostAllShip[i] = deepcopy(newq/sum(newparam.Q.ReversedCostAllShip))
            end
        end
        if paramchosen.ReversedTacticAllBoats=="dist"
            nbexp=deepcopy(allparam.Nbexp.ReversedDistAllShip[paramchosen.IndexReversedAllShip])
            newparam.AverageCost.ReversedDistAllShip[paramchosen.IndexReversedAllShip] = deepcopy((allparam.AverageCost.ReversedDistAllShip[paramchosen.IndexReversedAllShip]*nbexp+new_cost)/(nbexp+1))
            newparam.Nbexp.ReversedDistAllShip[paramchosen.IndexReversedAllShip]=deepcopy(nbexp+1)
            for i in 1:length(allparam.Nbexp.ReversedCostAllShip)
                newav = deepcopy(newparam.AverageCost.ReversedDistAllShip[i])
                newparam.Q.ReversedDistAllShip[i]=deepcopy(cost/newav)
            end
            for i in 1:length(allparam.Nbexp.ReversedCostAllShip)
                newq=deepcopy(newparam.Q.ReversedDistAllShip[i])
                newparam.Proba.ReversedDistAllShip[i] = deepcopy(newq/sum(newparam.Q.ReversedDistAllShip))
            end
        end
        if paramchosen.ReversedTacticAllBoats=="time"
            nbexp=deepcopy(allparam.Nbexp.ReversedTimeAllShip[paramchosen.IndexReversedAllShip])
            newparam.AverageCost.ReversedTimeAllShip[paramchosen.IndexReversedAllShip] = deepcopy((allparam.AverageCost.ReversedTimeAllShip[paramchosen.IndexReversedAllShip]*nbexp+new_cost)/(nbexp+1))
            newparam.Nbexp.ReversedTimeAllShip[paramchosen.IndexReversedAllShip]=deepcopy(nbexp+1)
            for i in 1:length(allparam.Nbexp.ReversedCostAllShip)
                newav = deepcopy(newparam.AverageCost.ReversedTimeAllShip[i])
                newparam.Q.ReversedTimeAllShip[i]=deepcopy(cost/newav)
            end
            for i in 1:length(allparam.Nbexp.ReversedCostAllShip)
                newq=deepcopy(newparam.Q.ReversedTimeAllShip[i])
                newparam.Proba.ReversedTimeAllShip[i] = deepcopy(newq/sum(newparam.Q.ReversedTimeAllShip))
            end
        end
    end

    if constrained
        nbexp=deepcopy(allparam.Nbexp.RateConstrained[paramchosen.IndexRateConstrained])
        newparam.AverageCost.RateConstrained[paramchosen.IndexRateConstrained] = deepcopy((allparam.AverageCost.RateConstrained[paramchosen.IndexRateConstrained]*nbexp+new_cost)/(nbexp+1))
        newparam.Nbexp.RateConstrained[paramchosen.IndexRateConstrained]=deepcopy(nbexp+1)
        for i in 1:length(allparam.Nbexp.RateConstrained)
            newav = deepcopy(newparam.AverageCost.RateConstrained[i])
            newparam.Q.RateConstrained[i]=deepcopy(cost/newav)
        end
        for i in 1:length(allparam.Nbexp.RateConstrained)
            newq=deepcopy(newparam.Q.RateConstrained[i])
            newparam.Proba.RateConstrained[i] = deepcopy(newq/sum(newparam.Q.RateConstrained))
        end
    end

    nbexp=deepcopy(allparam.Nbexp.PropToRemove[paramchosen.IndexPropToRemove])
    newparam.AverageCost.PropToRemove[paramchosen.IndexPropToRemove] = deepcopy((allparam.AverageCost.PropToRemove[paramchosen.IndexPropToRemove]*nbexp+new_cost)/(nbexp+1))
    newparam.Nbexp.PropToRemove[paramchosen.IndexPropToRemove]=deepcopy(nbexp+1)
    for i in 1:length(allparam.Nbexp.PropToRemove)
        newav = deepcopy(newparam.AverageCost.PropToRemove[i])
        newparam.Q.PropToRemove[i]=deepcopy(cost/newav)
    end
    for i in 1:length(allparam.Nbexp.PropToRemove)
        newq=deepcopy(newparam.Q.PropToRemove[i])
        newparam.Proba.PropToRemove[i] = deepcopy(newq/sum(newparam.Q.PropToRemove))
    end

    return newparam
end


function UpdateAfterHeurParameters(paramchosen::ChosenParameters, allparam::AllParameters, cost, new_cost, constrained)
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
    if paramchosen.Reversed=="yes"
        index=1
    end
    if paramchosen.Reversed=="no"
        index=2
    end
    nbexp = allparam.Nbexp.Reversed[index]
    newparam.AverageCost.Reversed[index] = deepcopy((allparam.AverageCost.Reversed[index]*nbexp+new_cost)/(nbexp+1))
    newparam.Nbexp.Reversed[index]=deepcopy(nbexp+1)
    for i in 1:2
        newav = deepcopy(newparam.AverageCost.Reversed[i])
        newparam.Q.Reversed[i]=deepcopy(cost/newav)
    end
    for i in 1:2
        newq=deepcopy(newparam.Q.Reversed[i])
        newparam.Proba.Reversed[i] = deepcopy(newq/sum(newparam.Q.Reversed))
    end


    if paramchosen.Reversed=="no"
        index=0
        if paramchosen.TacticAllBoats=="cost"
            index=1
        end
        if paramchosen.TacticAllBoats=="dist"
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
    else
        index=0
        if paramchosen.ReversedTacticAllBoats=="cost"
            index=1
        end
        if paramchosen.ReversedTacticAllBoats=="dist"
            index=2
        end
        if paramchosen.ReversedTacticAllBoats=="time"
            index=3
        end
        nbexp = deepcopy(allparam.Nbexp.ReversedTacticAllBoats[index])
        newparam.AverageCost.ReversedTacticAllBoats[index] = deepcopy((allparam.AverageCost.ReversedTacticAllBoats[index]*nbexp+new_cost)/(nbexp+1))
        newparam.Nbexp.ReversedTacticAllBoats[index]=deepcopy(nbexp+1)
        for i in 1:3
            newav = deepcopy(newparam.AverageCost.ReversedTacticAllBoats[i])
            newparam.Q.ReversedTacticAllBoats[i]=deepcopy(cost/newav)
        end
        for i in 1:3
            newq=deepcopy(newparam.Q.ReversedTacticAllBoats[i])
            newparam.Proba.ReversedTacticAllBoats[i] = deepcopy(newq/sum(newparam.Q.ReversedTacticAllBoats))
        end
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


    if paramchosen.Reversed=="no"
        if paramchosen.TacticAllBoats=="cost"
            nbexp=allparam.Nbexp.CostAllShip[paramchosen.IndexAllShip]
            newparam.AverageCost.CostAllShip[paramchosen.IndexAllShip] = (allparam.AverageCost.CostAllShip[paramchosen.IndexAllShip]*nbexp+new_cost)/(nbexp+1)
            newparam.Nbexp.CostAllShip[paramchosen.IndexAllShip]=nbexp+1
            for i in 1:length(allparam.Nbexp.CostAllShip)
                newav = deepcopy(newparam.AverageCost.CostAllShip[i])
                newparam.Q.CostAllShip[i]=cost/newav
            end
            for i in 1:length(allparam.Nbexp.CostAllShip)
                newq=deepcopy(newparam.Q.CostAllShip[i])
                newparam.Proba.CostAllShip[i] = newq/sum(newparam.Q.CostAllShip)
            end
        end

        if paramchosen.TacticAllBoats=="dist"
            nbexp=deepcopy(allparam.Nbexp.DistAllShip[paramchosen.IndexAllShip])
            newparam.AverageCost.DistAllShip[paramchosen.IndexAllShip] = deepcopy((allparam.AverageCost.DistAllShip[paramchosen.IndexAllShip]*nbexp+new_cost)/(nbexp+1))
            newparam.Nbexp.DistAllShip[paramchosen.IndexAllShip]=deepcopy(nbexp+1)
            for i in 1:length(allparam.Nbexp.CostAllShip)
                newav = deepcopy(newparam.AverageCost.DistAllShip[i])
                newparam.Q.DistAllShip[i]=deepcopy(cost/newav)
            end
            for i in 1:length(allparam.Nbexp.CostAllShip)
                newq=deepcopy(newparam.Q.DistAllShip[i])
                newparam.Proba.DistAllShip[i] = deepcopy(newq/sum(newparam.Q.DistAllShip))
            end
        end

        if paramchosen.TacticAllBoats=="time"
            nbexp=deepcopy(allparam.Nbexp.TimeAllShip[paramchosen.IndexAllShip])
            newparam.AverageCost.TimeAllShip[paramchosen.IndexAllShip] = deepcopy((allparam.AverageCost.TimeAllShip[paramchosen.IndexAllShip]*nbexp+new_cost)/(nbexp+1))
            newparam.Nbexp.TimeAllShip[paramchosen.IndexAllShip]=deepcopy(nbexp+1)
            for i in 1:length(allparam.Nbexp.CostAllShip)
                newav = deepcopy(newparam.AverageCost.TimeAllShip[i])
                newparam.Q.TimeAllShip[i]=deepcopy(cost/newav)
            end
            for i in 1:length(allparam.Nbexp.CostAllShip)
                newq=deepcopy(newparam.Q.TimeAllShip[i])
                newparam.Proba.TimeAllShip[i] = newq/sum(newparam.Q.TimeAllShip)
            end
        end
    else
        if paramchosen.ReversedTacticAllBoats=="cost"
            nbexp=deepcopy(allparam.Nbexp.ReversedCostAllShip[paramchosen.IndexReversedAllShip])
            newparam.AverageCost.ReversedCostAllShip[paramchosen.IndexAllShip] = deepcopy((allparam.AverageCost.ReversedCostAllShip[paramchosen.IndexReversedAllShip]*nbexp+new_cost)/(nbexp+1))
            newparam.Nbexp.ReversedCostAllShip[paramchosen.IndexReversedAllShip]=deepcopy(nbexp+1)
            for i in 1:length(allparam.Nbexp.ReversedCostAllShip)
                newav = deepcopy(newparam.AverageCost.ReversedCostAllShip[i])
                newparam.Q.ReversedCostAllShip[i]=deepcopy(cost/newav)
            end
            for i in 1:length(allparam.Nbexp.ReversedCostAllShip)
                newq=deepcopy(newparam.Q.ReversedCostAllShip[i])
                newparam.Proba.ReversedCostAllShip[i] = deepcopy(newq/sum(newparam.Q.ReversedCostAllShip))
            end
        end
        if paramchosen.ReversedTacticAllBoats=="dist"
            nbexp=deepcopy(allparam.Nbexp.ReversedDistAllShip[paramchosen.IndexReversedAllShip])
            newparam.AverageCost.ReversedDistAllShip[paramchosen.IndexReversedAllShip] = deepcopy((allparam.AverageCost.ReversedDistAllShip[paramchosen.IndexReversedAllShip]*nbexp+new_cost)/(nbexp+1))
            newparam.Nbexp.ReversedDistAllShip[paramchosen.IndexReversedAllShip]=deepcopy(nbexp+1)
            for i in 1:length(allparam.Nbexp.ReversedCostAllShip)
                newav = deepcopy(newparam.AverageCost.ReversedDistAllShip[i])
                newparam.Q.ReversedDistAllShip[i]=deepcopy(cost/newav)
            end
            for i in 1:length(allparam.Nbexp.ReversedCostAllShip)
                newq=deepcopy(newparam.Q.ReversedDistAllShip[i])
                newparam.Proba.ReversedDistAllShip[i] = deepcopy(newq/sum(newparam.Q.ReversedDistAllShip))
            end
        end
        if paramchosen.ReversedTacticAllBoats=="time"
            nbexp=deepcopy(allparam.Nbexp.ReversedTimeAllShip[paramchosen.IndexReversedAllShip])
            newparam.AverageCost.ReversedTimeAllShip[paramchosen.IndexReversedAllShip] = deepcopy((allparam.AverageCost.ReversedTimeAllShip[paramchosen.IndexReversedAllShip]*nbexp+new_cost)/(nbexp+1))
            newparam.Nbexp.ReversedTimeAllShip[paramchosen.IndexReversedAllShip]=deepcopy(nbexp+1)
            for i in 1:length(allparam.Nbexp.ReversedCostAllShip)
                newav = deepcopy(newparam.AverageCost.ReversedTimeAllShip[i])
                newparam.Q.ReversedTimeAllShip[i]=deepcopy(cost/newav)
            end
            for i in 1:length(allparam.Nbexp.ReversedCostAllShip)
                newq=deepcopy(newparam.Q.ReversedTimeAllShip[i])
                newparam.Proba.ReversedTimeAllShip[i] = deepcopy(newq/sum(newparam.Q.ReversedTimeAllShip))
            end
        end
    end

    if constrained
        nbexp=deepcopy(allparam.Nbexp.RateConstrained[paramchosen.IndexRateConstrained])
        newparam.AverageCost.RateConstrained[paramchosen.IndexRateConstrained] = deepcopy((allparam.AverageCost.RateConstrained[paramchosen.IndexRateConstrained]*nbexp+new_cost)/(nbexp+1))
        newparam.Nbexp.RateConstrained[paramchosen.IndexRateConstrained]=deepcopy(nbexp+1)
        for i in 1:length(allparam.Nbexp.RateConstrained)
            newav = deepcopy(newparam.AverageCost.RateConstrained[i])
            newparam.Q.RateConstrained[i]=deepcopy(cost/newav)
        end
        for i in 1:length(allparam.Nbexp.RateConstrained)
            newq=deepcopy(newparam.Q.RateConstrained[i])
            newparam.Proba.RateConstrained[i] = deepcopy(newq/sum(newparam.Q.RateConstrained))
        end
    end
    
    return newparam
end


function UpdateAfterLocalParameters(paramchosen::ChosenParameters, allparam::AllParameters, cost, new_cost)
    newparam=deepcopy(allparam)
    index=0
    if paramchosen.TacticLocalSearch=="cost"
        index=1
    end
    if paramchosen.TacticLocalSearch=="dist"
        index=2
    end
    if paramchosen.TacticLocalSearch=="time"
        index=3
    end
    nbexp = deepcopy(allparam.Nbexp.TacticLocalSearch[index])
    newparam.AverageCost.TacticLocalSearch[index] = deepcopy((allparam.AverageCost.TacticLocalSearch[index]*nbexp+new_cost)/(nbexp+1))
    newparam.Nbexp.TacticLocalSearch[index]=deepcopy(nbexp+1)
    for i in 1:3
        newav = deepcopy(newparam.AverageCost.TacticLocalSearch[i])
        newparam.Q.TacticLocalSearch[i]=deepcopy(cost/newav)
    end
    for i in 1:3
        newq=deepcopy(newparam.Q.TacticLocalSearch[i])
        newparam.Proba.TacticLocalSearch[i] = deepcopy(newq/sum(newparam.Q.TacticLocalSearch))
    end



    return newparam
end


function RestartParamNb(allparam)
    newparam=deepcopy(allparam)
    newparam.Nbexp.Reversed=fill(1,length(newparam.Nbexp.Reversed))
    newparam.Nbexp.TacticAllBoats=fill(1,length(newparam.Nbexp.TacticAllBoats))
    newparam.Nbexp.TacticOneBoat=fill(1,length(newparam.Nbexp.TacticOneBoat))
    newparam.Nbexp.ReversedTacticAllBoats=fill(1,length(newparam.Nbexp.ReversedTacticAllBoats))
    newparam.Nbexp.TacticLocalSearch=fill(1,length(newparam.Nbexp.TacticLocalSearch))

    newparam.Nbexp.TimeOneShip=fill(1,length(newparam.Nbexp.TimeOneShip))
    newparam.Nbexp.DistOneShip=fill(1,length(newparam.Nbexp.DistOneShip))
    newparam.Nbexp.CostOneShip=fill(1,length(newparam.Nbexp.CostOneShip))

    newparam.Nbexp.TimeAllShip=fill(1,length(newparam.Nbexp.TimeAllShip))
    newparam.Nbexp.DistAllShip=fill(1,length(newparam.Nbexp.DistAllShip))
    newparam.Nbexp.CostAllShip=fill(1,length(newparam.Nbexp.CostAllShip))

    newparam.Nbexp.ReversedTimeAllShip=fill(1,length(newparam.Nbexp.ReversedTimeAllShip))
    newparam.Nbexp.ReversedDistAllShip=fill(1,length(newparam.Nbexp.ReversedDistAllShip))
    newparam.Nbexp.ReversedCostAllShip=fill(1,length(newparam.Nbexp.ReversedCostAllShip))

    newparam.Nbexp.RateConstrained=fill(1,length(newparam.Nbexp.RateConstrained))

    newparam.Nbexp.PropToRemove=fill(1,length(newparam.Nbexp.PropToRemove))
    return newparam
end
