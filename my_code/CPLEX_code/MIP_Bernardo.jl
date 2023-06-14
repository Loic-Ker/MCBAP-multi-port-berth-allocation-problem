#include("../MBAP_INST.jl")
include("../MBAP_INST_CPLEX.jl")

function MIPmodel(N,Nout,seed,qli, tl)
    inst = readInstFromFile(location*"MCBAP-multi-port-berth-allocation-problem/data_small/CP2_Inst_$seed"*"_$N"*"_$Nout"*"_$qli"*".txt")
    @unpack N, Ntot, shipsIn, shipsOut, P, ports, visits, Pi, Nl, S, T, Ic, Dc, Hc, Fc, Pc, delta, gamma, dist, beta, Bp, qli, maxT = inst
    
    m = Model(CPLEX.Optimizer)
    set_silent(m)
    #set_optimizer_attribute(m, "CPXPARAM_Threads", NumberThreads)
    set_optimizer_attribute(m, "CPXPARAM_TimeLimit", tl)

    @variable(m, x[i in 1:Ntot, c in 1:visits[i] ] >= 1, Int) # leftmost berthing position
    @variable(m, y[i in 1:Ntot, c in 1:visits[i] ] >= 0, Int) # berthing start time
    @variable(m, v[i in 1:N, c in 1:visits[i]-1, s in 1:S], Bin) # 1 if speed s is chosen to travel between visit c and c+1 for ship i
    @variable(m, σ[i in 1:Ntot, j in 1:Ntot, c in 1:visits[i], c2 in 1:visits[j]; i != j && Pi[i][c] == Pi[j][c2] ], Bin) # 1 if ship i is to the left of ship j
    @variable(m, δ[i in 1:Ntot, j in 1:Ntot, c in 1:visits[i], c2 in 1:visits[j]; i != j && Pi[i][c] == Pi[j][c2]], Bin)  # 1 if ship i below ship j
    @variable(m, a[i in 1:N, c in 1:visits[i] ] >= 0) # arrival time
    @variable(m, d[i in 1:N, c in 1:visits[i] ] >= 0) # delay time
    @variable(m, u[i in 1:N, c in 1:visits[i] ] >= 0)# Bin) # penalty for exceeding latest finish time
    @variable(m, h[i in 1:Ntot, c in 1:visits[i] ] >= 0) # handling time
    @variable(m, r[i in 1:N, c in 1:visits[i]] >= 0) # deviation from ideal berthign position


    # define the big-M values and adapt lengths to quay discretization
    maxC = maximum(length.(Pi))
    Tm = zeros(Int64, Ntot, maxC)
    nq = Int64[]
    bi = zeros(Float64, N, maxC)
    for n in 1:N
        @unpack l, Bi = shipsIn[n]
        push!(nq, ceil(Int, l/qli))
        for (i,p) in enumerate(Pi[n])
            Tm[n,i] = T[n,i,2]
            bi[n,i] =  (Bi[i]/qli)+1 #round(Int, Bi[i]/qli)+1
        end
    end
    # fix external ships:
    count = length(shipsIn)
    for p in 1:P
        for ship in shipsOut[p]
            # println(ship)
            @unpack time, berth, hand, length = ship
            count += 1
            b = round(Int, berth/qli) + 1
            fix(x[count,1], b, force=true)
            fix(y[count,1], time, force=true)
            fix(h[count,1], hand, force=true)
            Tm[count, 1] = time + hand
            push!(nq, ceil(Int, length/qli))
        end
    end
    # we cannot use Tmax as M value if we allow to exceed the LFT
    Tmax = zeros(Int64, Ntot, Ntot, maxC, maxC)
    for i in 1:Ntot, j in 1:Ntot, c in 1:visits[i], c2 in 1:visits[j]
        if i != j && Pi[i][c] == Pi[j][c2]
            Tmax[i,j,c,c2] = max(Tm[i,c], Tm[j,c2])
        end
    end

    @objective(m, Min,
    sum( sum( Ic*(y[i,c] - a[i,c]) + Hc*h[i,c] + Dc*d[i,c] + Pc*u[i,c] for c in 1:visits[i])
    + sum( Fc*( dist[p,Pi[i][idx+1]]*sum(gamma[shipsIn[i].type,s]*v[i,idx,s] for s in 1:S)) for (idx,p) in enumerate(Pi[i][1:end-1]))
    for i in 1:N))

    @constraint(m, c1[i in 1:N, c in 1:visits[i]], x[i,c] + nq[i] - 1 <= Bp[Pi[i][c]]) # out ships already must satisfy this
    @constraint(m, c2[i in 1:Ntot, j in 1:Ntot, c in 1:visits[i], c2 in 1:visits[j]; i != j && Pi[i][c] == Pi[j][c2]], x[i,c] + nq[i] <= x[j,c2] + (Bp[Pi[i][c]]+1)*(1 - σ[i,j,c,c2]))
    @constraint(m, c3[i in 1:Ntot, j in 1:Ntot, c in 1:visits[i], c2 in 1:visits[j]; i != j && Pi[i][c] == Pi[j][c2]], y[i,c] + h[i,c] <= y[j,c2] + maxT*(1 - δ[i,j,c,c2])) # Tmax[i,j,c,c2] max(T[i,c,2], T[j,c2,2])
    @constraint(m, c4a[i in 1:N, j in 1:Ntot, c in 1:visits[i], c2 in 1:visits[j]; i < j && Pi[i][c] == Pi[j][c2]], σ[i,j,c,c2] + σ[j,i,c2,c] + δ[i,j,c,c2] + δ[j,i,c2,c] >= 1)
    # Assumption: all shipsIn visit at least 2 ports (one travel time minimum)
    @constraint(m, c5[i in 1:N, c in 1:visits[i]-1], y[i,c] + h[i,c] + dist[Pi[i][c],Pi[i][c+1]]*sum(delta[s]*v[i,c,s] for s in 1:S) == a[i,c+1] )
    @constraint(m, c6[i in 1:N, c in 1:visits[i]], a[i,c] <= y[i,c])
    @constraint(m, c7[i in 1:N, c in 1:visits[i]], shipsIn[i].sT[c] <= y[i,c])
    @constraint(m, c8[i in 1:N, c in 1:visits[i]], y[i,c] + h[i,c] - shipsIn[i].eT[c] <= d[i,c])
    @constraint(m, c8b[i in 1:N, c in 1:visits[i]], y[i,c] + h[i,c] <= T[i,c,2] + u[i,c])
    @constraint(m, c9[i in 1:N, c in 1:visits[i]], (1 + beta*qli*(r[i,c]))*ports[Pi[i][c]].minH[shipsIn[i].type] == h[i,c])
    @constraint(m, c10[i in 1:N, c in 1:visits[i]], x[i,c] - bi[i,c] <= r[i,c])
    @constraint(m, c11[i in 1:N, c in 1:visits[i]], bi[i,c] - x[i,c] <= r[i,c])
    @constraint(m, c12[i in 1:N, c in 1:visits[i]-1], sum(v[i,c,s] for s in 1:S) == 1)

    optimize!(m)
    status = termination_status(m)
    print('\n')
    print("The status :")
    print('\n')
    print(status)
    return JuMP.objective_value(m)
end