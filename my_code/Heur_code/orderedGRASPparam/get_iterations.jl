function prepareSolIterSoft(seed,N,Nout,qli, nb_iter, inst, sol, cost, allparam, paramfixed, expname)
    @unpack N, Ntot, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT, Nl = inst
    d = Dict{Any,Any}()
    
    ### Here add: 
        #the type of the boats with every visit
        #The cost divided in total and for each visit added
        #The position of each visit added
        # Parameters at the end of the solution make
        # Which tactic for each boat

    d["cost_solheur"] = sol.store.costHeur.all
    d["delay_cost_solheur"] = sol.store.costHeur.delay
    d["waiting_cost_solheur"] = sol.store.costHeur.waiting
    d["penalty_solheur"] = sol.store.costHeur.penalty
    d["handling_cost_solheur"] = sol.store.costHeur.handling
    d["fuel_cost_solheur"] = sol.store.costHeur.fuel

    d["cost_sollocal"] = sol.store.costLocal.all
    d["delay_cost_sollocal"] = sol.store.costLocal.delay
    d["waiting_cost_sollocal"] = sol.store.costLocal.waiting
    d["penalty_sollocal"] = sol.store.costLocal.penalty
    d["handling_cost_sollocal"] = sol.store.costLocal.handling
    d["fuel_cost_sollocal"] = sol.store.costLocal.fuel

    d["timeheur"]=sol.store.timeHeur
    d["timelocal"]=sol.store.timeLocalSearch

    list_tacticboat=sol.store.parameters.TacticOneBoat
    list_tacticall = sol.store.parameters.TacticAllBoats
    list_tacticlocalsearch = sol.store.parameters.TacticLocalSearch
    
    for el in 1:length(list_tacticboat)
        list_tacticboat[el] = round(list_tacticboat[el], digits=4)
    end
    for el in 1:length(list_tacticall)
        list_tacticall[el] = round(list_tacticall[el], digits=4)
    end
    for el in 1:length(list_tacticlocalsearch)  
        list_tacticlocalsearch[el] = round(list_tacticlocalsearch[el], digits=4)
    end

    d["proba_tacticboat"] = sol.store.parameters.TacticOneBoat
    d["proba_tacticall"] = sol.store.parameters.TacticAllBoats
    d["proba_tacticlocalsearch"] = sol.store.parameters.TacticLocalSearch

    d["oneboatdistance"]=allparam.Proba.DistOneShip
    d["oneboatcost"]=allparam.Proba.CostOneShip
    d["oneboattime"]=allparam.Proba.TimeOneShip
    d["allboatscost"]=allparam.Proba.CostAllShip
    d["allboatscount"]=allparam.Proba.CountAllShip
    d["allboatstime"]=allparam.Proba.TimeAllShip
    d["localsearchboat"]=paramfixed.LocalSearchBoat
    d["localsearchrandom"]=paramfixed.LocalSearchRandom

    d["rateconstrained"]=allparam.Proba.RateConstrained
    
    d["failed"] = sol.failed
    d["better"] = sol.better
    d["reconstruct"]=sol.reconstruct
    d["usedCPLEX"]=sol.usedCPLEX

    ## Get the instances 
    d["inst"]= inst.Pi

    #print("Some proba :")
    #print('\n')
    #print("Tactic boat (should not evolve) :")
    #print(list_tacticboat)
    #print('\n')
    #print("Time proba one ship (should evolve) :")
    #print(allparam.Proba.TimeOneShip)
    #print('\n')
    
    ## Get the port calls
    T_vec=Vector{Vector{Vector{Int64}}}()
    for n in 1:N
        push!(T_vec,Vector{Vector{Int64}}(undef, 0))
        for c in 1:length(inst.Pi[n])
            push!(T_vec[n],T[n,c,:])
        end
    end
    d["calls"]=T_vec

    ## Get the boats length :
    Nlen=Vector{Int64}()
        for n in Nl
            push!(Nlen,ceil(Int, n/qli))
        end
    d["length_boats"]= Nlen
    
    cost_visit = Vector{Vector{Any}}(undef, 0)
    delay_cost_visit = Vector{Vector{Any}}(undef, 0)
    waiting_cost_visit = Vector{Vector{Any}}(undef, 0)
    penalty_visit = Vector{Vector{Any}}(undef, 0)
    handling_cost_visit = Vector{Vector{Any}}(undef, 0)
    fuel_cost_visit = Vector{Vector{Any}}(undef, 0)
    tacticboat_chosen = Vector{Vector{Any}}(undef, 0)
    tacticall_chosen = Vector{Vector{Any}}(undef, 0)
    when = Vector{Vector{Any}}(undef, 0)
    for n in 1:N
        cost_thisboat = Vector{Any}(undef, 0)
        delay_cost_thisboat = Vector{Any}(undef, 0)
        waiting_cost_thisboat = Vector{Any}(undef, 0)
        penalty_thisboat = Vector{Any}(undef, 0)
        handling_cost_thisboat = Vector{Any}(undef, 0)
        fuel_cost_thisboat = Vector{Any}(undef, 0)
        tacticboat_thisboat = Vector{Any}(undef, 0)
        tacticall_thisboat = Vector{Any}(undef, 0)
        when_thisboat = Vector{Any}(undef, 0)
        for (c,p) in enumerate(inst.Pi[n])
            push!(cost_thisboat, sol.visits[n][c].store.cost.all)
            push!(delay_cost_thisboat, sol.visits[n][c].store.cost.delay)
            push!(waiting_cost_thisboat, sol.visits[n][c].store.cost.waiting)
            push!(penalty_thisboat, sol.visits[n][c].store.cost.penalty)
            push!(handling_cost_thisboat, sol.visits[n][c].store.cost.handling)
            push!(fuel_cost_thisboat, sol.visits[n][c].store.cost.fuel)
            push!(tacticboat_thisboat, sol.visits[n][c].store.tacticBoat)
            push!(tacticall_thisboat, sol.visits[n][c].store.tacticAll)
            push!(when_thisboat, sol.visits[n][c].store.when)
        end
        push!(cost_visit,cost_thisboat)
        push!(delay_cost_visit,delay_cost_thisboat)
        push!(waiting_cost_visit,waiting_cost_thisboat)
        push!(penalty_visit,penalty_thisboat)
        push!(fuel_cost_visit,fuel_cost_thisboat)
        push!(handling_cost_visit,handling_cost_thisboat)
        push!(tacticboat_chosen, tacticboat_thisboat)
        push!(tacticall_chosen, tacticall_thisboat)
        push!(when, when_thisboat)
    end

    d["objectif"]=cost
    d["cost_visit"]=cost_visit
    d["delay_cost_visit"]=delay_cost_visit
    d["waiting_cost_visit"]=waiting_cost_visit
    d["handling_cost_visit"]=handling_cost_visit
    d["fuel_cost_visit"]=fuel_cost_visit
    d["penalty_visit"]=penalty_visit
    d["tacticboat_chosen"]=tacticboat_chosen
    d["tacticall_chosen"]=tacticall_chosen
    d["when"]=when

    return d
end


function prepareSolIter(seed,N,Nout,qli, nb_iter, inst, sol, cost, allparam, paramfixed, expname)
    @unpack N, Ntot, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT, Nl = inst
    d = Dict{Any,Any}()
    
    ### Here add: 
        #the type of the boats with every visit
        #The cost divided in total and for each visit added
        #The position of each visit added
        # Parameters at the end of the solution make
        # Which tactic for each boat

    d["cost_solheur"] = sol.store.costHeur.all
    d["delay_cost_solheur"] = sol.store.costHeur.delay
    d["waiting_cost_solheur"] = sol.store.costHeur.waiting
    d["penalty_solheur"] = sol.store.costHeur.penalty
    d["handling_cost_solheur"] = sol.store.costHeur.handling
    d["fuel_cost_solheur"] = sol.store.costHeur.fuel

    d["cost_sollocal"] = sol.store.costLocal.all
    d["delay_cost_sollocal"] = sol.store.costLocal.delay
    d["waiting_cost_sollocal"] = sol.store.costLocal.waiting
    d["penalty_sollocal"] = sol.store.costLocal.penalty
    d["handling_cost_sollocal"] = sol.store.costLocal.handling
    d["fuel_cost_sollocal"] = sol.store.costLocal.fuel

    d["timeheur"]=sol.store.timeHeur
    d["timelocal"]=sol.store.timeLocalSearch

    list_tacticboat=sol.store.parameters.TacticOneBoat
    list_tacticall = sol.store.parameters.TacticAllBoats
    list_tacticlocalsearch = sol.store.parameters.TacticLocalSearch
    
    for el in 1:length(list_tacticboat)
        list_tacticboat[el] = round(list_tacticboat[el], digits=4)
    end
    for el in 1:length(list_tacticall)
        list_tacticall[el] = round(list_tacticall[el], digits=4)
    end
    for el in 1:length(list_tacticlocalsearch)  
        list_tacticlocalsearch[el] = round(list_tacticlocalsearch[el], digits=4)
    end

    d["proba_tacticboat"] = sol.store.parameters.TacticOneBoat
    d["proba_tacticall"] = sol.store.parameters.TacticAllBoats
    d["proba_tacticlocalsearch"] = sol.store.parameters.TacticLocalSearch

    d["oneboatdistance"]=allparam.Proba.DistOneShip
    d["oneboatcost"]=allparam.Proba.CostOneShip
    d["oneboattime"]=allparam.Proba.TimeOneShip
    d["allboatscost"]=allparam.Proba.CostAllShip
    d["allboatscount"]=allparam.Proba.CountAllShip
    d["allboatstime"]=allparam.Proba.TimeAllShip
    d["localsearchboat"]=paramfixed.LocalSearchBoat
    d["localsearchrandom"]=paramfixed.LocalSearchRandom

    d["rateconstrained"]=allparam.Proba.RateConstrained
    
    d["failed"] = sol.failed
    d["better"] = 0
    d["reconstruct"]=sol.reconstruct
    d["usedCPLEX"]=sol.usedCPLEX

    ## Get the instances 
    d["inst"]= inst.Pi

    #print("Some proba :")
    #print('\n')
    #print("Tactic boat (should not evolve) :")
    #print(list_tacticboat)
    #print('\n')
    #print("Time proba one ship (should evolve) :")
    #print(allparam.Proba.TimeOneShip)
    #print('\n')
    
    ## Get the port calls
    T_vec=Vector{Vector{Vector{Int64}}}()
    for n in 1:N
        push!(T_vec,Vector{Vector{Int64}}(undef, 0))
        for c in 1:length(inst.Pi[n])
            push!(T_vec[n],T[n,c,:])
        end
    end
    d["calls"]=T_vec

    ## Get the boats length :
    Nlen=Vector{Int64}()
        for n in Nl
            push!(Nlen,ceil(Int, n/qli))
        end
    d["length_boats"]= Nlen
    
    x = Vector{Vector{Any}}(undef, 0)
    y = Vector{Vector{Any}}(undef, 0)
    hand = Vector{Vector{Any}}(undef, 0)
    v = Vector{Vector{Any}}(undef, 0)
    a = Vector{Vector{Any}}(undef, 0)
    cost_visit = Vector{Vector{Any}}(undef, 0)
    delay_cost_visit = Vector{Vector{Any}}(undef, 0)
    waiting_cost_visit = Vector{Vector{Any}}(undef, 0)
    penalty_visit = Vector{Vector{Any}}(undef, 0)
    handling_cost_visit = Vector{Vector{Any}}(undef, 0)
    fuel_cost_visit = Vector{Vector{Any}}(undef, 0)
    tacticboat_chosen = Vector{Vector{Any}}(undef, 0)
    tacticall_chosen = Vector{Vector{Any}}(undef, 0)
    when = Vector{Vector{Any}}(undef, 0)
    for n in 1:N
        x_thisboat= Vector{Any}(undef, 0)
        y_thisboat= Vector{Any}(undef, 0)
        hand_thisboat= Vector{Any}(undef, 0)
        v_thisboat= Vector{Any}(undef, 0)
        a_thisboat= Vector{Any}(undef, 0)
        cost_thisboat = Vector{Any}(undef, 0)
        delay_cost_thisboat = Vector{Any}(undef, 0)
        waiting_cost_thisboat = Vector{Any}(undef, 0)
        penalty_thisboat = Vector{Any}(undef, 0)
        handling_cost_thisboat = Vector{Any}(undef, 0)
        fuel_cost_thisboat = Vector{Any}(undef, 0)
        tacticboat_thisboat = Vector{Any}(undef, 0)
        tacticall_thisboat = Vector{Any}(undef, 0)
        when_thisboat = Vector{Any}(undef, 0)
        for (c,p) in enumerate(inst.Pi[n])
            push!(x_thisboat, sol.visits[n][c].b)
            push!(y_thisboat, sol.visits[n][c].t)
            push!(hand_thisboat, h[n][c][sol.visits[n][c].b])
            push!(cost_thisboat, sol.visits[n][c].store.cost.all)
            push!(delay_cost_thisboat, sol.visits[n][c].store.cost.delay)
            push!(waiting_cost_thisboat, sol.visits[n][c].store.cost.waiting)
            push!(penalty_thisboat, sol.visits[n][c].store.cost.penalty)
            push!(handling_cost_thisboat, sol.visits[n][c].store.cost.handling)
            push!(fuel_cost_thisboat, sol.visits[n][c].store.cost.fuel)
            push!(tacticboat_thisboat, sol.visits[n][c].store.tacticBoat)
            push!(tacticall_thisboat, sol.visits[n][c].store.tacticAll)
            push!(when_thisboat, sol.visits[n][c].store.when)
            if c < length(inst.Pi[n])
                pN = sol.visits[n][c+1].p
                bN = sol.visits[n][c+1].b
                tN = sol.visits[n][c+1].t
                # (pN,bN,tN) = sch[c+1]
                deltaT = tN - (sol.visits[n][c].t +  h[n][c][sol.visits[n][c].b])
                s = findLowestSpeed(deltaT, delta, dist[sol.visits[n][c].p,pN])
                push!(v_thisboat, ceil(Int, s))
            end
            if c>1
                pN = sol.visits[n][c-1].p
                bN = sol.visits[n][c-1].b
                tN = sol.visits[n][c-1].t
                # (pN,bN,tN) = sch[c+1]
                deltaT = sol.visits[n][c].t - (tN +  h[n][c-1][bN])
                s = findLowestSpeed(deltaT, delta, dist[pN,sol.visits[n][c].p])
                if s!=-1
                    push!(a_thisboat, tN + h[n][c-1][bN] + delta[s]*dist[pN,sol.visits[n][c].p])
                end
            end
            if c==1
                push!(a_thisboat,sol.visits[n][c].t)
            end
        end
        push!(x,x_thisboat)
        push!(y,y_thisboat)
        push!(hand,hand_thisboat)
        push!(v,v_thisboat)
        push!(a,a_thisboat)
        push!(cost_visit,cost_thisboat)
        push!(delay_cost_visit,delay_cost_thisboat)
        push!(waiting_cost_visit,waiting_cost_thisboat)
        push!(penalty_visit,penalty_thisboat)
        push!(fuel_cost_visit,fuel_cost_thisboat)
        push!(handling_cost_visit,handling_cost_thisboat)
        push!(tacticboat_chosen, tacticboat_thisboat)
        push!(tacticall_chosen, tacticall_thisboat)
        push!(when, when_thisboat)
    end


    shipsOutFlatt=collect(Iterators.flatten(shipsOut))
    for n in N+1:Ntot
        x_thisboat= Vector{Any}(undef, 0)
        y_thisboat= Vector{Any}(undef, 0)
        hand_thisboat= Vector{Any}(undef, 0)
        for (c,p) in enumerate(inst.Pi[n])
            push!(x_thisboat, round(shipsOutFlatt[n-N].berth/qli)+1)
            push!(y_thisboat, shipsOutFlatt[n-N].time)
            push!(hand_thisboat, shipsOutFlatt[n-N].hand)
        end
        push!(x,x_thisboat)
        push!(y,y_thisboat)
        push!(hand,hand_thisboat)
    end          
    d["x"]=x
    d["y"]=y
    d["a"]=a
    d["v"]=v
    d["hand"]=hand
    d["objectif"]=cost
    d["cost_visit"]=cost_visit
    d["delay_cost_visit"]=delay_cost_visit
    d["waiting_cost_visit"]=waiting_cost_visit
    d["handling_cost_visit"]=handling_cost_visit
    d["fuel_cost_visit"]=fuel_cost_visit
    d["penalty_visit"]=penalty_visit
    d["tacticboat_chosen"]=tacticboat_chosen
    d["tacticall_chosen"]=tacticall_chosen
    d["when"]=when

    return d
end

