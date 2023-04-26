
mutable struct Visit
    # n::Int64    # ship number
    # c::Int64    # visit number
    planned::Bool # planned visit?
    p::Int64    # port
    b::Int64    # berth
    t::Int64    # time
    minT::Int64 # minimum berthing time given previous or next visit
    maxT::Int64 # maximum berthing time given previous or next visit
end
Visit() = Visit(false, -1,-1,-1,0,typemax(Int64))

mutable struct Sol
    obj::Float64
    # schedules::Vector{Schedule}
    M::Vector{Vector{Array{Bool, 2}}}
    # unassigned::Vector{Tuple{Int64,Int64}}
    visits::Vector{Vector{Visit}}
end
Sol(inst::Instance) = Sol(0., [[ones(Bool, inst.Bp[inst.Pi[n][c]]-ceil(Int, inst.shipsIn[n].l/inst.qli)+1, ceil(Int, inst.maxT+1)) for c in 1:length(inst.Pi[n])] for n in 1:inst.N], [[Visit() for c in 1:length(inst.Pi[n])] for n in 1:inst.N])

## A solution is a matrix of one of length of : number of discretized berthing positions minus size of ship * maximum time limit for all boat and ports,
## Then we have this for each visit of each boat (a bit overkilled no ?).
## We also keep a list of the visit per port.
instance = readInstFromFile("D:/DTU-Courses/DTU-Thesis/berth_allocation/data_small/CP2_Inst_1_4_3_10.txt")

