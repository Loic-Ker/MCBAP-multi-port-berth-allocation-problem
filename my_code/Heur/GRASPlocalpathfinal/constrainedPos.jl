using UnPack 
using Statistics
import XLSX
using CSV, Tables
using DataFrames
using Random
using StatsBase
using Distributions
using JuMP, CPLEX
include("../../MBAP_INST.jl")
include("../toolsMatrixTimes.jl")
include("../check_solution.jl")
include("../utilInit.jl")


function findConstrainedPos(inst::Instance, sol::Sol, port::Int, t1::Int, t2::Int, n::Int, c::Int, l::Int)
    # find ship visits during t1-t2
    pos = Tuple{Int, Int}[]
    for n_ in 1:inst.N
        for (c_,p) in enumerate(inst.Pi[n_])
            if p == port
                if sol.visits[n_][c_].planned
                    @unpack b,t = sol.visits[n_][c_]
                    hand = ceil(Int, inst.h[n_][c_][b])
                    if doOverlap(t1, t2, t, t + hand)
                        len = ceil(Int, inst.shipsIn[n_].l/inst.qli)
                        # check 12 possible positions (16-4)
                        # this positions may still be infeasible
                        X = [b - l, b, b + len - l, b + len]
                        # Y = [t - h, t, t + hand - h, t + hand]
                        for x in b-l:b+len
                            if 1 <= x && x + l <= inst.Bp[p]
                                h = ceil(Int, inst.h[n][c][x])
                                for y in [t - h, t + hand]
                                    if t1 <= y && y + h <= t2
                                        
                                        push!(pos, (x, y))
                                    end
                                end
                            end
                        end
                        for x in [b - l, b - l-1, b + len, b+len+1]
                            if 1 <= x && x + l <= inst.Bp[p]
                                h = ceil(Int, inst.h[n][c][x])
                                for y in t-h:t+hand
                                    if t1 <= y && y + h <= t2
                                        push!(pos, (x, y))
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    # also consider external ships and quay limits
    for extS in inst.shipsOut[port]
        @unpack time, berth, hand, length = extS
        if doOverlap(t1, t2, time, time + hand)
            # check 12 possible positions (16-4)
            # this positions may still be infeasible
            X = [berth - l, berth, berth + length - l, berth + length]
            # Y = [t - h, t, t + hand - h, t + hand]
            for (ix, x) in enumerate(X)
                if 1 <= x && x + l <= inst.Bp[port]
                    h = ceil(Int, inst.h[n][c][x])
                    if ix in [2,3]
                        if t1 <= time - h && time <= t2
                            push!(pos, (x, time - h))
                        end
                        if t1 <= time + hand && time + hand + h <= t2
                            push!(pos, (x, time + hand))
                        end
                    end
                    if ix in [1,4]
                        if t1 <= time && time + h <= t2
                            push!(pos, (x, time))
                        end
                        if t1 <= time + hand - h && time + hand <= t2
                            push!(pos, (x, time + hand - h))
                        end
                    end
                end
            end
        end
    end
    return pos
end