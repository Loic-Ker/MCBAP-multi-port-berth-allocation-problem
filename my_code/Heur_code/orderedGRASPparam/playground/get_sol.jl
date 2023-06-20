include("parametric_free.jl")


#instance = readInstFromFile("/zhome/c3/6/164957/code_git/MCBAP-multi-port-berth-allocation-problem/data_small/CP2_Inst_2_4_3_10.txt")
#sol = greedyrandomizedconstruction(instance, false, 5)


function prepareSol(inst, sol, cost)
    @unpack N, Ntot, P, Pi, visits, shipsIn, shipsOut, h, dist, delta, qli, T, Bp, maxT, Nl = inst
    d = Dict{Any,Any}()
    
    ### Here add: 
        #the type of the boats with every visit
        #The cost divided in total and for each visit added
        #The position of each visit added
        # Parameters at the end of the solution make
        # Which tactic for each boat


    ## Get the instances 
    d["inst"]= inst.Pi

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
    for n in 1:N
        x_thisboat= Vector{Any}(undef, 0)
        y_thisboat= Vector{Any}(undef, 0)
        hand_thisboat= Vector{Any}(undef, 0)
        v_thisboat= Vector{Any}(undef, 0)
        a_thisboat= Vector{Any}(undef, 0)
        for (c,p) in enumerate(inst.Pi[n])
            push!(x_thisboat, sol.visits[n][c].b)
            push!(y_thisboat, sol.visits[n][c].t)
            push!(hand_thisboat, h[n][c][sol.visits[n][c].b])
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
    return d
end



newbenchmark = [Any["Seed","N","Nout","qli","OldLB","OldUB","OldTime","NewCost","MyTime"]]

import XLSX


function makeExpText(type1, type2, type3, paramfixed, list_paramvisit, list_paramconstrained, time_local, max_time_heur, max_time, expname, location)
    
    filename = location*"results_jobs/benchmarks_HEUR/orderedGRASPparam/$expname"*"/explanations.txt"
    types = "There are no fiwed types but three possible parameters"
    txttype1 = "The tactic for one boat is $type1"
    txttype2 = "The tactic for all boats is $type2"
    txttype3 = "The tactic for the ocal search is $type3"

    window=paramfixed.WindowSize
    txtwindow = "The window size for the visits to look at is $window"

    lookforconstraint=paramfixed.lookforconstrained
    txtlookforconstraint = "Should we look for constrained at the beginning : $lookforconstraint"

    alphaboat=paramfixed.LocalSearchBoat
    alpharandom=paramfixed.LocalSearchRandom
    txtalphaboat = "For the local search with the boats we take $alphaboat boats"
    txtalpharandom = "For the local search with random visits we take $alpharandom visits"

    proptoremove=paramfixed.PropToRemove
    txtproptoremove = "The prop of ships to remove for the reconstruction : $proptoremove"

    pushatconstraint=paramfixed.PushOnlyConstrained
    txtpushatconstraint = "Should we push at constraints position : $pushatconstraint"

    
    txtoneship = "The alpha parameters for one ship : "
    for el in list_paramvisit
        txtoneship = txtoneship*"$el "
    end
    txtallship = "The alpha parameters for all ships : "
    for el in list_paramconstrained
        txtallship = txtallship*"$el "
    end


    txttimelocal = "The maximum time for the local search : $time_local "
    txttimeheur = "The maximum time for the heuristic : $max_time_heur "
    txttimemax = "The time of the exp : $max_time"


    open(filename, "w") do file
        write(file,'\n')
        write(file, txtlookforconstraint)
        write(file,'\n')
        write(file, txtwindow)
        write(file,'\n')
        write(file, txtproptoremove)
        write(file,'\n')
        write(file, txtpushatconstraint)
        write(file,'\n')
        write(file, txtalphaboat)
        write(file,'\n')
        write(file, txtalpharandom)
        write(file,'\n')
        write(file, txttype1)
        write(file,'\n')
        write(file, txttype2)
        write(file,'\n')
        write(file, txttype3)
        write(file,'\n')
        write(file, txtoneship)
        write(file,'\n')
        write(file, txtallship)
        write(file,'\n')
        write(file, txttimeheur)
        write(file,'\n')
        write(file, txttimelocal)
        write(file,'\n')
        write(file, txttimemax)
        write(file,'\n')
    end
end


function makeSolHeur(type1, type2, type3, paramfixed, temperature, time_local, max_time_heur, max_time, expname, location, seed_chosen, Nchosen, Noutchosen, qlichosen)
    xf = CSV.read(location*"MCBAP-multi-port-berth-allocation-problem/Small_Inst_Res.csv", DataFrame)
    newbenchmark = DataFrame(Seed= [0],N= [0],Nout= [0],qli= [0],HeurCost= [0])
    all_instances = readdir(location*"MCBAP-multi-port-berth-allocation-problem/Large")
    for instance_name in all_instances
	split_instance = split(instance_name,"_")
	seed=parse(Int64,split_instance[3])
	N=parse(Int64,split_instance[4])
	Nout=parse(Int64,split_instance[5])
	qli=parse(Int64,split(split_instance[6],".")[1])
	if seed ==seed_chosen && N==Nchosen && Nout==Noutchosen && qli==qlichosen
		inst = readInstFromFile(location*"MCBAP-multi-port-berth-allocation-problem/Large/CP2_Inst_$seed"*"_$N"*"_$Nout"*"_$qli"*".txt")
		print("The instance : $seed"*"_$N"*"_$Nout"*"_$qli")
		if isdir(location*"results_jobs/benchmarks_HEUR/orderedGRASPparam/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli")==false
		    mkdir(location*"results_jobs/benchmarks_HEUR/orderedGRASPparam/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli")
		end
		if isdir(location*"results_jobs/benchmarks_HEUR/orderedGRASPparam/$expname"*"/iterations_before_local/sol_$seed"*"_$N"*"_$Nout"*"_$qli")==false
		    mkdir(location*"results_jobs/benchmarks_HEUR/orderedGRASPparam/$expname"*"/iterations_before_local/sol_$seed"*"_$N"*"_$Nout"*"_$qli")
		end
		sol, cost, allparam = GRASP_reactive(seed,N,Nout,qli, type1, type2, type3, paramfixed, temperature, time_local, max_time_heur, max_time, expname, location)
		#print('\n')
		#print("The solution :")
		#print('\n')
		#print(sol.visits)
		print('\n')
		print("And the cost is ")
		print('\n')
		print(cost)
		feasible=true
		for n in 1:N
		    # The times :
		    for (c,p) in enumerate(inst.Pi[n])
		        if sol.visits[n][c].planned == false
		            feasible = false              
		        end
		    end
		end
		if feasible && checkSolutionFeasability(inst, sol)
		    d=prepareSol(inst, sol, cost)
		    CSV.write(location*"results_jobs/benchmarks_HEUR/orderedGRASPparam/$expname"*"/final_sols/sol_$seed"*"_$N"*"_$Nout"*"_$qli"*".csv", d)
		end
	    
		this_benchmark=DataFrame(Seed= [seed],N= [N],Nout= [Nout],qli= [qli],HeurCost= [ ceil(Int, cost)])
		newbenchmark=append!(newbenchmark,this_benchmark)
	end
    end
    return newbenchmark
end

location = "D:/DTU-Courses/DTU-Thesis/berth_allocation/"
#location="/zhome/c3/6/164957/code_git/"


# The parameters of the experiment :

# The experience name :
expname="playground"

# The tactic types :
type1="time" 
type2="cost" 
type3="random"

# The window size for the visits to look at :
window=0.1

# The number of boat to remove for the local search :
alphaboat=5
alpharandom=23

# The prop of ships to remove for the reconstruction :
proptoremove=0.001

# Should we push at constraints position :
pushatconstraint=true

# Look for constrained at the beginning :
lookforconstraint=true

# All the parameters :
paramfixed=FixedParameters(alpharandom,alphaboat,proptoremove,window,pushatconstraint, lookforconstraint)

# Maximum time for the local search :
time_local=5
#time_local= parse(Int64,ARGS[1])

# Maximum time for the heuristic :
max_time_heur=30

# maximum time for the experiment :
max_time=2400
#max_time = parse(Int64,ARGS[2])

# the temperature parameter :
temperature=0.93

# look for a specific seed
#seedchosen = parse(Int64,ARGS[3])
#Nchosen=parse(Int64,ARGS[4])
#Noutchosen=parse(Int64,ARGS[5])
seedchosen=2
Nchosen=50
Noutchosen=10
qlichosen=10

allparam = initializeParam()
list_paramvisit = allparam.Alpha.CostOneShip
list_paramconstrained = allparam.Alpha.RateConstrained
makeExpText(type1, type2, type3, paramfixed, list_paramvisit, list_paramconstrained, time_local, max_time_heur, max_time, expname, location)
newbenchmark = makeSolHeur(type1, type2, type3, paramfixed, temperature, time_local, max_time_heur, max_time, expname, location, seedchosen, Nchosen, Noutchosen, qlichosen)
CSV.write(location*"results_jobs/benchmarks_HEUR/orderedGRASPparam/$expname"*"/NLarge_playground_$seedchosen"*"n_$Nchosen"*"_$Noutchosen.csv", newbenchmark)
newbenchmark


