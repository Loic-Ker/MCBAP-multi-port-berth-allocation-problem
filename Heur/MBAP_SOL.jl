include("../MBAP_INST.jl")


mutable struct NewVisit
    n::Int64    # ship number
    c::Int64    # visit number
    b::Int64  #postion
    t::Int64    # port
    cost::Float64    
    distance::Float64
    time::Float64   
    constrained::Bool
    store::ToStoreVisit 
end

mutable struct Probabilities
    TacticOneBoat::Vector{Float64}
    TacticAllBoats::Vector{Float64}
    TacticLocalSearch::Vector{Float64}

    TimeOneShip::Vector{Float64}
    DistOneShip::Vector{Float64}
    CostOneShip::Vector{Float64}

    TimeAllShip::Vector{Float64}
    DistAllShip::Vector{Float64}
    CostAllShip::Vector{Float64}

    RateConstrained::Vector{Float64}

    PropToRemove::Vector{Float64}
end

mutable struct SplitCosts
    all::Int64
    delay::Int64
    waiting::Int64
    penalty::Int64
    handling::Int64
    fuel::Int64    
end

mutable struct ToStoreVisit
    cost::SplitCosts
    when::Int64
    tacticBoat::String
    tacticAll::String
end


mutable struct ToStoreSol
    costHeur::SplitCosts
    costLocal::SplitCosts
    timeHeur::Float64
    timeLocalSearch::Float64
    parameters::Probabilities
end

mutable struct Visit
    # n::Int64    # ship number
    # c::Int64    # visit number
    planned::Bool # planned visit?
    p::Int64    # port
    b::Int64    # berth
    t::Int64    # time
    minT::Int64 # minimum berthing time given previous or next visit
    maxT::Int64 # maximum berthing time given previous or next visit
    store::ToStoreVisit
    failed::Int64
end

Visit() = Visit(false, -1,-1,-1,0,typemax(Int64), ToStoreVisit(SplitCosts(0,0,0,0,0,0),0,"",""), 0)

mutable struct Sol
    obj::Float64
    # schedules::Vector{Schedule}
    M::Vector{Vector{Array{Bool, 2}}}
    # unassigned::Vector{Tuple{Int64,Int64}}
    visits::Vector{Vector{Visit}}
    store::ToStoreSol
    failed::Int64
    reconstruct::Int64
    relinking::Int64
    usedLocalSearch::Int64
    better::Int64
    average_cost_elite::Float64
    average_dist_elite::Float64
    pushimprove::Int64
    countlocal::Int64
    pathcost
    countpath
end



Sol(inst::Instance) = Sol(0., [[ones(Bool, inst.Bp[inst.Pi[n][c]]-ceil(Int, inst.shipsIn[n].l/inst.qli)+1, ceil(Int, inst.maxT+1)) for c in 1:length(inst.Pi[n])] for n in 1:inst.N], [[Visit() for c in 1:length(inst.Pi[n])] for n in 1:inst.N], ToStoreSol(SplitCosts(0,0,0,0,0,0),SplitCosts(0,0,0,0,0,0),0.0,0.0,Probabilities([1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0],[1.0], [1.0], [1.0], [1.0])),0,0,0,0,0,0,0,0,0,0,0)

