using UnPack 
using Statistics
include("../MBAP_INST.jl")
include("MBAP_SOL.jl")


# compute the cost of ship n at port visit c berthing at pos b and time t given solution sol for the solution
function computeCostPosSol(inst::Instance, n::Int64, c::Int64, b::Int64, t::Int64 ,sol::Sol)
    @unpack P, h, Pc, shipsIn, T, delta, gamma, dist, Hc, Dc, Fc, Ic = inst
    @unpack type, eT = shipsIn[n]
    p = inst.Pi[n][c]
    cost = Hc*h[n][c][b]
    handling_cost=Hc*h[n][c][b]
    delay = t + h[n][c][b] > eT[c] ? t + h[n][c][b] - eT[c] : 0.
    delay_cost=Dc*delay
    cost += Dc*delay
    penalty = t + h[n][c][b] > T[n,c,2] ? Pc*(t + h[n][c][b] - T[n,c,2]) : 0.
    cost += penalty
    fuel_cost=0
    waiting_cost=0
    #########
    ## Remove this part for later : do not take into account the fuel
    if c < length(sol.visits[n]) && sol.visits[n][c+1].planned
        pN = sol.visits[n][c+1].p
        bN = sol.visits[n][c+1].b
        tN = sol.visits[n][c+1].t
        # (pN,bN,tN) = sch[c+1]
        deltaT = tN - (t + h[n][c][b])
        s = findLowestSpeed(deltaT, delta, dist[p,pN])
        if s == -1
            # this should only happen if we quit the port too late
            s = length(delta)
            return 0,0,0,0,0,0, false
        end
        ## Change here, we do not want to take into account in the optimization the 
        ## differences in distances between the different ports
        cost += Fc*gamma[type,s]*dist[p,pN]
        fuel_cost=Fc*gamma[type,s]*dist[p,pN]
    end
    if c > 1 && sol.visits[n][c-1].planned
        pP = sol.visits[n][c-1].p
        bP = sol.visits[n][c-1].b
        tP = sol.visits[n][c-1].t
        hP = h[n][c-1][bP]
        deltaT = t - (tP + hP)
        s = findLowestSpeed(deltaT, delta, dist[pP,p])
        if s == -1
            # this should only happen if we quit the port too late
            s = length(delta)
            return 0,0,0,0,0,0, false
        end
        ## Change here, we do not want to take into account in the optimization the 
        ## differences in distances between the different ports
        cost += Fc*gamma[type,s]*dist[pP,p]
        fuel_cost=Fc*gamma[type,s]*dist[pP,p]
        wait = t - (tP + hP + delta[s]*dist[pP,p])
        waiting_cost=Ic*wait
        cost += Ic*wait
    end
    return cost, delay_cost, waiting_cost, penalty, handling_cost, fuel_cost, true
end


function checkSolutionCost(inst::Instance, sol::Sol)
    # returns the cost of a given solution
    @unpack Fc, Ic, Dc, Hc, Pc, h, dist, delta, gamma, shipsIn, T = inst
    cost = 0.
    delay_cost=0.
    penalty_cost=0.
    handling_cost=0.
    fuel_cost=0.
    waiting_cost=0.
    for (n, sch) in enumerate(sol.visits)
        @unpack type, eT = shipsIn[n]
        for (c, vis) in enumerate(sch)
            @unpack p,b,t = vis
            # (p,b,t) = vis
            handling_cost += Hc*h[n][c][b]
            cost += Hc*h[n][c][b]
            delay = t + h[n][c][b] > eT[c] ? t + h[n][c][b] - eT[c] : 0.
            delay_cost += Dc*delay
            cost += Dc*delay
            penalty = t + h[n][c][b] > T[n,c,2] ? Pc*(t + h[n][c][b] - T[n,c,2]) : 0.
            penalty_cost+=penalty
            cost += penalty
            if c < length(sch)
                pN = sch[c+1].p
                bN = sch[c+1].b
                tN = sch[c+1].t
                # (pN,bN,tN) = sch[c+1]
                deltaT = tN - (t + h[n][c][b])
                s = findLowestSpeed(deltaT, delta, dist[p,pN])
                #TODO: FIX!!
                if s == -1
                    println("Infeasible speed!!!!")
                    s = length(delta)
                else
                    fuel_cost+=Fc*gamma[type,s]*dist[p,pN]
                    cost += Fc*gamma[type,s]*dist[p,pN]
                    wait = tN - (t + h[n][c][b] + delta[s]*dist[p,pN]) > 0 ? tN - (t + h[n][c][b] + delta[s]*dist[p,pN]) : 0.
                    waiting_cost+=Ic*wait
                    cost += Ic*wait
                end
            end
        end
    end
    return cost, delay_cost, waiting_cost, penalty_cost, handling_cost, fuel_cost
end

function checkSolutionFeasability(inst::Instance, sol::Sol)
    @unpack N, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp = inst
    ## We are gonna check everything
    feasible = true


    ## No conflicts with within one boat schedule :
    for n in 1:N
        # The times :
        for (c,p) in enumerate(inst.Pi[n])
            t1 = sol.visits[n][c].t
            t2 = t1 +  h[n][c][sol.visits[n][c].b]

            # if previous visit is planned, check current visit can be reached with max speed
            if c > 1 
                pP = sol.visits[n][c-1].p
                bP = sol.visits[n][c-1].b
                tP = sol.visits[n][c-1].t
                # (pP,bP,tP) = sol.visits[n][c-1]
                hP = h[n][c-1][bP]
                tA = ceil(Int, tP + hP + delta[end]*dist[pP,p])
                if tA>t1
                    print('\n')
                    print(tP)
                    print('\n')
                    print(tA)
                    print('\n')
                    print(t1)
                    print('\n')
                    print("error1")
                    feasible=false
                end
            end
            
            # if next visit is planned, check next visit time can be reached with max speed
            if c < length(Pi[n])
                pN = sol.visits[n][c+1].p
                bN = sol.visits[n][c+1].b
                tN = sol.visits[n][c+1].t
                # (pN,bN,tN) = sol.solution[n][c+1]
                tD = floor(Int, tN - delta[end]*dist[p,pN])
                if t2>tD
                    print("error2")
                    feasible=false
                end
            end
        end
    end

    # Let's check each boat has correct boxs within it's own schedule (no conflicts with it's own + berth limits)
    for n in 1:N
        l = ceil(Int, shipsIn[n].l/qli)
        for (c,p) in enumerate(inst.Pi[n])
            t1 = sol.visits[n][c].t
            t2 = t1 +  h[n][c][sol.visits[n][c].b]
            x1 = sol.visits[n][c].b
            x2 = x1 + l
            
            if x1<0 || x2>Bp[p] || x2<0
                print('\n')
                print("error3")
                print('\n')
                print(Bp[p])
                print('\n')
                print(x2)
                feasible = false
            end

            for (c_,p_) in enumerate(inst.Pi[n])
                if p==p_ && c!=c_
                    t1_ = sol.visits[n][c_].t
                    t2_ = t1 +  h[n][c_][sol.visits[n][c_].b]
                    x1_ = sol.visits[n][c_].b
                    x2_ = x1_ + l
                    if doOverlapRectangle(t1, t2, x1, x2, t1_, t2_, x1_, x2_)
                        print("error4")
                        feasible=false
                    end
                end
            end
        end
    end


    ## Now let's compare with the other boats
    for n in 1:inst.N-1
        l = ceil(Int, shipsIn[n].l/qli)
        for n_ in n+1:inst.N
            l_ = ceil(Int, shipsIn[n_].l/qli)
            for (c,p) in enumerate(inst.Pi[n])
                t1 = sol.visits[n][c].t
                t2 = t1 +  h[n][c][sol.visits[n][c].b]
                x1 = sol.visits[n][c].b
                x2 = x1 + l
                for (c_,p_) in enumerate(inst.Pi[n_])
                    if p==p_
                        t1_ = sol.visits[n_][c_].t
                        t2_ = t1_ +  h[n_][c_][sol.visits[n_][c_].b]
                        x1_ = sol.visits[n_][c_].b
                        x2_ = x1_ + l_
                        if doOverlapRectangle(t1, t2, x1, x2, t1_, t2_, x1_, x2_)
                            print('\n')
                            print("error5")
                            print('\n')
                            print("The port $p")
                            print('\n')
                            print('\n')
                            print("The boat $n , $c and $n_ , $c_")
                            print('\n')
                            print("The first box : [$t1,$t2] and [$x1,$x2]")
                            print('\n')
                            print("The second box : [$t1_,$t2_] and [$x1_,$x2_]")

                            feasible=false
                        end
                    end
                end
            end
        end
    end

    ## With the external boats 
    for n in 1:N
        l = ceil(Int, shipsIn[n].l/qli)
        for (c,p) in enumerate(inst.Pi[n])
            t1 = sol.visits[n][c].t
            t2 = t1 +  h[n][c][sol.visits[n][c].b]
            x1 = sol.visits[n][c].b
            x2 = x1 + l
            for extS in inst.shipsOut[p]
                @unpack time, berth, hand, length = extS
                l_ = ceil(Int, length/qli)
                t1_ = time
                t2_= t1_ + hand
                x1_ = round(Int, berth/qli) + 1
                x2_ = x1_ + l_
                if doOverlapRectangle(t1, t2, x1, x2, t1_, t2_, x1_, x2_)
                    print("error6")
                    feasible=false
                end
            end
        end
    end

    return feasible
end
    
