using UnPack 
include("MBAP_SOL.jl")


# check if two periods overlap
function doOverlap(t, tF, t_, tF_)
    if t_ < t < tF_
        return true
    elseif t_ < tF < tF_
        return true
    elseif t < t_ < tF
        return true
    elseif t < tF_ < tF
        return true
    elseif t == t_ || tF == tF_
        return true
    end
    return false
end

# check if two rectangles overlap : for the final check function
function doOverlapRectangle(x1l, x1r, y1d, y1u, x2l, x2r, y2d, y2u)
    # check if one element is a line
    if x1l == x1r || y1u == y1d || x2l == x2r || y2u == y2d
        return false
    end
    # If one rectangle is on left side of other
    if x1l >= x2r || x2l >= x1r
        return false
    end
    # If one rectangle is above other
    if y1d >= y2u || y2d >= y1u
        return false
    end
    return true
end

# compute the lowest feasible speed to satisfy travel time t
function findLowestSpeed(t, delta::Vector{Float64}, dist::Int64)
    S = length(delta)
    for s in 1:S
        tt = delta[s]*dist
        if tt < t
            return s
        end
    end
    return -1
end

# generate matrix M of occupied positions
function generateOccupiedMx(inst::Instance, sol) #, unassigned::Vector{Tuple{Int64, Int64}})
    @unpack N, P, Pi, T, Bp, h, dist, qli, shipsIn, shipsOut, delta, maxT, maxC = inst
    ## The matrix we create contains all the time/distance discretization
    M = [[zeros(Bool, Bp[Pi[n][c]]-ceil(Int, shipsIn[n].l/qli)+1, maxT+1) for c in 1:length(Pi[n])] for n in 1:N]
    ## The between our bounded time intervalls we set the positions as possible (true)
    for i in 1:N
        nq = ceil(Int,shipsIn[i].l/qli)
        for (c,p) in enumerate(Pi[i])
            t1 = max(shipsIn[i].sT[c], T[i,c,1])
            t2 = maxT 
            for b in 1:Bp[p]-nq+1
                hand = ceil(Int, h[i][c][b])
                M[i][c][b,t1+1:t2-hand+2] .= true
            end
        end
    end

    ## Here we remove the positions if the places are already taken by other boats
    ## We use the same tactic that the updateMatrix function 
    for port in 1:P
        for el in shipsOut[port]
            @unpack time, berth, hand, length = el
            b = round(Int, berth/qli) + 1
            nq = ceil(Int, length/qli)
            for n in 1:N, (c,p) in enumerate(Pi[n])
                ## If it's in the same port
                if port == p
                    l = ceil(Int, shipsIn[n].l/qli)
                    ## We can not not position our new boat between b_boat-length_newboat
                    bMin = b - l + 1 < 0.5 ? 1 : b - l + 1
                    ## And b_boat+length_boat if the have the same time line
                    bMax = b + nq - 1 > Bp[p]-l+1 ?  Bp[p]-l+1 : b + nq - 1
                    for b_ in bMin:bMax
                        ## Then here we check the time line
                        ## For each of these position we need the time to not be between :
                        handling = ceil(Int, h[n][c][b_])
                        ## t-handling_newboat
                        tMin = time - handling + 1 < -0.1 ? 0 : time - handling + 1
                        ## and t+handling_boat
                        tMax = time + hand - 1
                        ## Then we set everything to false
                        M[n][c][b_,tMin+1:tMax+1] .= false
                    end
                end
            end
        end
    end

    ## Here we can also create the matrix from a previous constructed solution
    ## This could be great to use for memory algorithm ?
    for n in 1:N
        for (c,val) in enumerate(sol[n])
            if val.planned # !((n,c) in unassigned)
                @unpack p,b,t = val
                # (p,b,t) = val
                nq = ceil(Int,shipsIn[n].l/qli)
                hand = ceil(Int, h[n][c][b])
                M = updateMpositions(inst, n, c, b, t, nq, hand, M)
            end
        end
    end
    return M
end



function updateTimesInitialization(inst::Instance, sol::Sol)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp = inst
    for n in 1:N
        l = ceil(Int, shipsIn[n].l/qli)
        len = length(Pi[n])
        if len>1
            for c in 0:len-2
                tf = sol.visits[n][len-c].maxT
                handling=maximum(h[n][len-c])
                besttime=tf-handling
                if besttime<0
                    besttime=0
                end
                for b in 1:Bp[Pi[n][len-c]]-l
                    newhandling = h[n][len-c][b]
                    time = ceil(Int,tf - newhandling)
                    if time<0
                        time=0
                    end
                    if sol.M[n][len-c][b,time+1]
                        if time+1>besttime+1
                            besttime=time
                        end
                    end
                end
                sol.visits[n][len-c-1].maxT = ceil(Int,minimum([besttime-Int(ceil(delta[end]*dist[Pi[n][len-c-1], Pi[n][len-c]])),sol.visits[n][len-c-1].maxT]))
            end
            
        end
    end
    return sol
end

# update matrix M of occupied positions, on the same port 
function updateMpositions(inst::Instance, n::Int64, c::Int64, b::Int64, t::Int64, nq::Int64, hand::Int64, M)
    @unpack N, Pi, h, shipsIn, Bp, qli = inst
    ## We look at all the other boats and all their visits
    for n_ in 1:N, (c_,p_) in enumerate(Pi[n_])
        ## If it's the same port, and we have a different boat or a different visit for the same boat
        if ((n_ == n && c != c_) || (n_ != n)) && p_ == Pi[n][c]
            l = ceil(Int, shipsIn[n_].l/qli)
            ## We can not not position our new boat between b_boat-length_newboat
            bMin = b - l + 1 < 0.5 ? 1 : b - l + 1
            ## And b_boat+length_boat if the have the same time line
            bMax = b + nq - 1 > Bp[p_]-l+1 ?  Bp[p_]-l+1 : b + nq - 1
            for b_ in bMin:bMax
                ## Then here we check the time line
                ## For each of these position we need the time to not be between :
                handling = ceil(Int, h[n_][c_][b_])
                ## t-handling_newboat
                tMin = t - handling + 1 < -0.1 ? 0 : t - handling + 1
                ## and t+handling_boat
                tMax = t + hand - 1
                ## Then we set everything to false
                M[n_][c_][b_,tMin+1:tMax+1] .= false
            end
        end
    end
    return M
end

# update matrix M of occupied times, regarding conflicts with the same boat
function updateMtimes(inst::Instance, visits::Vector{Vector{Visit}}, n::Int64, c::Int, M)
    @unpack N, Pi, h, shipsIn, Bp, qli = inst
    l = ceil(Int, shipsIn[n].l/qli)
    for (c_,p_) in enumerate(Pi[n])
        if c_!=c
            tmin = visits[n][c_].minT
            tmax = visits[n][c_].maxT
            for t in 1:tmin
                for b_ in 1:Bp[p_]-l
                    M[n][c_][b_,t] = false
                end
            end
            for b_ in 1:Bp[p_]-l
                hand = h[n][c_][b_]
                startime = tmax-hand
                if startime<1
                    startime=1
                end
                M[n][c_][b_,ceil(Int,startime):end] .= false
            end
            for t in tmin:tmax
                for b_ in 1:Bp[p_]-l
                    hand = h[n][c_][b_]
                    if t+1+hand>tmax
                        M[n][c_][b_,t+1] = false
                    end
                end
            end
        end
    end
    return M
end

function updateTimesAfterVisit(inst::Instance, visits::Vector{Vector{Visit}}, n::Int, c::Int, b::Int, t::Int)
    @unpack N, P, Pi, shipsIn, shipsOut, h, dist, delta, qli, T, Bp = inst
    l = ceil(Int, shipsIn[n].l/qli)
    len = length(Pi[n])
    if c>1
        visits[n][c-1].maxT = ceil(Int,minimum([t-Int(ceil(delta[end]*dist[Pi[n][c-1], Pi[n][c]])),visits[n][c-1].maxT]))
    end
    if c<len
        visits[n][c+1].minT = ceil(Int,maximum([t+h[n][c][b]+Int(ceil(delta[end]*dist[Pi[n][c+1], Pi[n][c]])),visits[n][c+1].minT]))
    end
    return visits
end