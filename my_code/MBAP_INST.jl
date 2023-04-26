using UnPack 

# internal ship (to be optimized)
mutable struct Ship
    type::Int64 # (1:small, 2:medium, 3:large). This affects the fuel consumption parameter and the minimum handling time at each port
    l::Int64    # length (in meters)
    Bi::Vector{Int64} # desired position in quay
    sT::Vector{Int64} # min start time at each visit
    eT::Vector{Int64} # expected finish time (at each visit)
end
# external ship (fixed)
mutable struct ExternalShip
    time::Int64 # berthing time
    berth::Int64 # berthing position
    hand::Int64 # handling time (Previously Float64)
    length::Int64 # length (in meters)
end
# terminal specs
mutable struct Term
    # ships::Vector{Int64}
    Ql::Int64   # length of quay (in meters)
    minH::Vector{Float64}             # minimum handling time per ship type
end
mutable struct Terminal
    id::Int64
    name::String
    coords::Tuple{Float64,Float64}
    quays::Vector{Float64}
end
# Instance object
mutable struct Instance
    N::Int64                          # number of (internal) ships
    Ntot::Int64                       # all ships (in + out)
    shipsIn::Vector{Ship}             # internal ships specs
    shipsOut::Vector{Vector{ExternalShip}}  # one vector of external ships per port
    P::Int64                          # total number of ports
    ports::Vector{Term}               # specs of each port (terminal)
    visits::Vector{Int64}             # number of port visits for each ship (in and out)
    Pi::Vector{Vector{Int64}}         # array of ports visited by each ship (in visit order)
    Nl::Vector{Int64}                 # length (in meters) of each ship
    qli::Int64                        # discretization size of the quay (multiples of 10m)
    beta::Float64                     # relative increase in handling time from desired position (Bi) in meters
    dist::Array{Int64,2}              # distance between ports #TODO: have another parameter for initial dist from port
    Ic::Int64                         # idleness cost ($/hour)
    Hc::Int64                         # handling cost ($/hour)
    Dc::Int64                         # delay cost    ($/hour)
    Fc::Int64                         # fuel cost     ($/ton)
    Pc::Int64                         # penalty cost after LFS ($/hour)
    T::Array{Int64,3}                 # first and last time instant at each port visit for each ship
    # Needed for MIP
    S::Int64                          # Set of speeds
    delta::Vector{Float64}            # travel time per unit of distance based on speed
    gamma::Array{Float64,2}           # fuel consumed per unit of distance based on speed and ship
    # Used for graph (can be integrated in Term and Ships objects)
    Bp::Vector{Int64}                 # number of berthing positions at each port (taking discretization into account)
    # Bx::Vector{Int64}
    h::Vector{Vector{Vector{Float64}}}# handling time for each ship, port visit and berthing position
    xP::Array{Int64,2}                # number of times that ship i visits port p
    maxT::Int64                       # maximum plannign horizon
    maxC::Int64                       # maximum number of visits made by a ship
    minB::Vector{Vector{Int}}         # leftmost feasible position for ship n at visit c (already in qli segments)
    maxB::Vector{Vector{Int}}         # rightmost feasible position for ship n at visit c (already in qli segments)
end

# default constructor (add seed if randomized)
Instance(N::Int64, P::Int64, S::Int64, qli::Int64) = Instance(
                                                                N,
                                                                0,
                                                                Vector{Ship}(),
                                                                [Vector{ExternalShip}() for p in 1:P],
                                                                P,
                                                                Vector{Term}(),
                                                                Vector{Int64}(),
                                                                Vector{Vector{Int64}}(),
                                                                Vector{Int64}(),
                                                                qli,
                                                                0.,
                                                                zeros(Int64, P, P),
                                                                200,
                                                                200,
                                                                300,
                                                                500,
                                                                2000,
                                                                zeros(Int64, N, 4, 2), # we assume max 4 visits per route
                                                                S, # hard-coded
                                                                zeros(Float64, S),
                                                                zeros(Float64, N, S),
                                                                zeros(Int64, P),
                                                                [Vector{Vector{Float64}}() for i in 1:N],
                                                                zeros(Int64, N, P),
                                                                0,
                                                                0,
                                                                [Vector{Vector{Int}}() for i in 1:N],
                                                                [Vector{Vector{Int}}() for i in 1:N]
)
# function that tightens the possible time windows of each ship given the earliest and latest times
function updateTvals(inst::Instance)
    @unpack N, Pi, shipsIn, h, dist, delta, Bp, T = inst
    for n in 1:N
        @unpack type, l, sT, eT, Bi = shipsIn[n]
        te = sT[1]
        for (i,p) in enumerate(Pi[n])
            # tE[i] = te
            # inst.T[n,i,1] = te > inst.T[n,i,1] ? te : inst.T[n,i,1]
            inst.T[n,i,1] = max(te, sT[i], T[n,i,1])
            te = inst.T[n,i,1]
            minH = minimum(h[n][i])
            if i < length(Pi[n])
                te += Int(ceil(minH + delta[end]*dist[p, Pi[n][i+1]]))
            end
        end
        # now backwards (this can be merged with the previous for loop, faster)
        tf = inst.T[n,length(Pi[n]),2]
        for c in 1:length(Pi[n])
            i = length(Pi[n])-c + 1
            inst.T[n,i,2] = tf < inst.T[n,i,2] ? tf : inst.T[n,i,2]
            tf = inst.T[n,i,2]
            p = Pi[n][i]
            minH = minimum(h[n][i])
            tf -= Int(ceil(minH))
            # tF[i] = tf
            if i > 1
                tf -= Int(ceil(delta[end]*dist[Pi[n][i-1], p]))
            end
        end
    end
end
# compute the bounds for each ship at each visit based on their idela position and the parameter DEVx
function computeBerthBounds(inst::Instance, DEVx::Float64)
    @unpack N, Pi, shipsIn, qli, ports = inst
    minB = [zeros(Int, length(Pi[n])) for n in 1:N]
    maxB = [zeros(Int, length(Pi[n])) for n in 1:N]
    for n in 1:N
        toMove = shipsIn[n].l*DEVx
        for (c,p) in enumerate(Pi[n])
            leftB = max(0, shipsIn[n].Bi[c] - toMove)
            rightB = min(ports[p].Ql, shipsIn[n].Bi[c] + shipsIn[n].l + toMove)
            minB[n][c] = round(Int, leftB/qli)+1
            maxB[n][c] = floor(Int, rightB/qli) #TODO this does not seem completely accurate
        end
    end
    return minB, maxB
end
# reads instance from a txt file (Ex: Instances/Medium/CP2_Inst_1_25_3_10.txt )
function readInstFromFile(filename::String)
    f = open(filename,"r")
    lines = readlines(f)
    line_count = 1
    vals = parse.(Int64,split(strip(lines[line_count]),"\t"))
    N, Nout, qli, P, S = vals
    inst = Instance(N,P,S,qli)
    line_count += 1
    inst.beta = parse(Float64, lines[line_count])
    line_count += 1
    for n in 1:N
        # 1. line: type, l
        vals = parse.(Int64,split(strip(lines[line_count]),"\t"))
        type, l = vals
        line_count += 1
        # 2. line: Bi
        Bi = parse.(Int64,split(strip(lines[line_count]),"\t"))
        line_count += 1
        # 3. line: sT
        sT = parse.(Int64,split(strip(lines[line_count]),"\t"))
        line_count += 1
        # 4. line: eT
        eT = parse.(Int64,split(strip(lines[line_count]),"\t"))
        line_count += 1
        # 5. line: Pi
        Pi = parse.(Int64,split(strip(lines[line_count]),"\t"))
        line_count += 1
        ship = Ship(type, l, Bi, sT, eT)
        push!(inst.shipsIn, ship)
        push!(inst.Pi, Pi)
        push!(inst.visits, length(Pi))
        push!(inst.Nl, ship.l)
        for el in Pi
            inst.xP[n,el] += 1
        end
    end
    inst.Ntot = N
    for p in 1:P
        # 1. line: Ql and num of outShips
        vals = parse.(Int64,split(strip(lines[line_count]),"\t"))
        Ql, nOP = vals
        inst.Ntot += nOP
        line_count += 1
        # 2. line: min H
        minH = parse.(Float64,split(strip(lines[line_count]),"\t"))
        line_count += 1
        term = Term(Ql, minH)
        push!(inst.ports, term)
        # next nOP lines: each out ship
        for n in 1:nOP
            vals = parse.(Int64,split(strip(lines[line_count]),"\t"))
            time, berth, hand, length = vals
            line_count += 1
            outShip = ExternalShip(time, berth, hand, length)
            push!(inst.shipsOut[p], outShip)
            push!(inst.visits, 1)
            push!(inst.Pi, [p])
            push!(inst.Nl, outShip.length)
        end
    end
    for p in 1:P
        vals = parse.(Int64,split(strip(lines[line_count]),"\t"))
        line_count += 1
        for p2 in 1:P
            inst.dist[p,p2] = vals[p2]
        end
    end
    vals = parse.(Int64,split(strip(lines[line_count]),"\t"))
    Ic, Hc, Dc, Fc = vals
    line_count += 1
    for n in 1:N
        vals = parse.(Int64,split(strip(lines[line_count]),"\t"))
        line_count += 1
        for c in 1:4
            inst.T[n,c,1] = vals[c]
        end
    end
    for n in 1:N
        vals = parse.(Int64,split(strip(lines[line_count]),"\t"))
        line_count += 1
        for c in 1:4
            #NOTE: Temp
            # inst.T[n,c,2] = ceil(Int, 2.5*vals[c])
            inst.T[n,c,2] = vals[c]
        end
    end
    inst.delta = parse.(Float64,split(strip(lines[line_count]),"\t"))
    line_count += 1
    for type in 1:3
        inst.gamma[type,:] = parse.(Float64,split(strip(lines[line_count]),"\t"))
        line_count += 1
    end
    close(f)

    # This is for the graph (check H same)
    ## Compute the number of berthing position
    for p=1:P
        inst.Bp[p] = Int(floor(inst.ports[p].Ql/qli))
    end
    #h based on c or p (at first I'd say c, we also have Bi)
    for i=1:N, c in 1:inst.visits[i]
        p = inst.Pi[i][c]
        vecH = Float64[]
        idealPos = (inst.shipsIn[i].Bi[c]/qli)+1 # Bi is in meters b=1 is 0-10m
        type = inst.shipsIn[i].type
        for b in 1:inst.Bp[p]
            push!(vecH, (1 + inst.beta*abs(b-idealPos)*qli)*inst.ports[p].minH[type])
            # push!(inst.h[i,p], (1 + inst.beta*abs(b-idealPos)*qli)*inst.ports[p].minH[i])
        end
        push!(inst.h[i], vecH)
        #NOTE: Temp
        #inst.T[i,c,2] = ceil(Int, 2.5*inst.T[i,c,2])
    end
    updateTvals(inst)
    #NOTE: Temp
    inst.maxT = ceil(Int, maximum(inst.T[:,:,2]))
    # inst.maxT = ceil(Int, maximum(inst.T[:,:,2])*1.5) # same as in BP method
    inst.maxC = maximum(length.(inst.Pi))
    return inst
end

