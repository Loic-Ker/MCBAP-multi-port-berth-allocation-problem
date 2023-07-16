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

function makeExpText(temperature, paramfixed, time_local, max_time_heur, max_time, expname, location)
    
    filename = location*"results_jobs/benchmarks_HEUR/finalGRASP/$expname"*"/explanations.txt"
    type1 = paramfixed.OneBoat
    type2 = paramfixed.AllBoat
    type3 = paramfixed.ReversedAllBoat
    type4 = paramfixed.LocalSearchOne
    type7 = paramfixed.LocalSearchAll
    txttype1 = "The tactic for one boat is $type1"
    thisoneboat=paramfixed.Alphaoneboat
    txtalphatype1 = "The alpha parameters possible are $thisoneboat"
    window=paramfixed.WindowSize
    txtwindow = "The window size for the visits to look at is $window"
    txttype2 = "The tactic for all boats is $type2"
    thisalphaboat = paramfixed.Alphaallboat
    txtalphatype2 = "The alpha parameters possible are $thisalphaboat"
    thisreverse=paramfixed.Reversed
    txttype3bis = "Do we reverse the order of the visits : $thisreverse (all is both)"
    txttype3 = "The tactic for all boats in reversed is $type3"
    thisallboats=paramfixed.Alphareversedallboat
    txtalphatype3 = "The alpha parameters possible are $thisallboats"
    txttype4 = "The tactic for the local search one boat : $type4"
    txttype7 = "The tactic for the local search all boats : $type7"

    alphaboatmin=paramfixed.LocalSearchBoat
    alpharandommin=paramfixed.LocalSearchRandom
    # Which sol to take from the heuristic for the local search :
    windowlocalsearch=paramfixed.windowLocalSearch
    txtwindowlocalsearch = "The window cost size for the local search (which one to take) : $windowlocalsearch"
    txtalphaboat = "For the local search we take a proportion of $alphaboatmin worst ships by cost to removre"
    txtalpharandom = "For the local search we take a proportion of $alpharandommin random boats to remove"
    pushatconstraint=paramfixed.PushOnlyConstrained
    txtpushatconstraint = "Should we push only at constraints position during the time push of local search : $pushatconstraint"

    lookforconstraint=paramfixed.lookforconstrained
    txtlookforconstraint = "Should we look for constrained at the beginning : $lookforconstraint"
    type5 = paramfixed.Alpharateconstrained
    txttype5 = "The alpha parameters for the proportion of constrained are $type5"

    # When do we start from zero for the reconstruction :
    restartparams=paramfixed.restartParams
    txtrestartparams = "After n iterations we restart the probabilities : $restartparams"

    txttimelocal = "The maximum time for the local search : $time_local "
    txttimeheur = "The maximum time for the heuristic : $max_time_heur "
    timerelinking = paramfixed.MaxTimeRelinking
    txttimerelinking = "The maximum time for the relinking : $timerelinking "
    txttimemax = "The time of the exp : $max_time"



    open(filename, "w") do file
        write(file,'\n')
        write(file,txttype1)
        write(file,'\n')
        write(file,txtalphatype1)
        write(file,'\n')
        write(file,txtwindow)
        write(file,'\n')
        write(file,txttype2)
        write(file,'\n')
        write(file,txtalphatype2)
        write(file,'\n')
        write(file,txttype3bis)
        write(file,'\n')
        write(file,txttype3)
        write(file,'\n')
        write(file,txtalphatype3)
        write(file,'\n')
        write(file,txttype4)
        write(file,'\n')
        write(file,txttype7)
        write(file,'\n')
        write(file,txtwindowlocalsearch)
        write(file,'\n')
        write(file,txtalphaboat)
        write(file,'\n')
        write(file,txtalpharandom)
        write(file,'\n')
        write(file,txtpushatconstraint)
        write(file,'\n')
        write(file,txtlookforconstraint)
        write(file,'\n')
        write(file,txttype5)
        write(file,'\n')
        write(file,txtrestartparams)
        write(file,'\n')
        write(file,txttimelocal)
        write(file,'\n')
        write(file,txttimeheur)
        write(file,'\n')
        write(file,txttimerelinking)
        write(file,'\n')
        write(file,txttimemax)
        write(file,'\n')
    end
end



function makeSolHeur(paramfixed, temperature, time_local, max_time_heur, max_time, expname, location, seed_chosen, Nchosen, Noutchosen, qlichosen)
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
		if isdir(location*"results_jobs/benchmarks_HEUR/finalGRASP/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli")==false
		    mkdir(location*"results_jobs/benchmarks_HEUR/finalGRASP/$expname"*"/iterations/sol_$seed"*"_$N"*"_$Nout"*"_$qli")
		end
		if isdir(location*"results_jobs/benchmarks_HEUR/finalGRASP/$expname"*"/iterations_before_local/sol_$seed"*"_$N"*"_$Nout"*"_$qli")==false
		    mkdir(location*"results_jobs/benchmarks_HEUR/finalGRASP/$expname"*"/iterations_before_local/sol_$seed"*"_$N"*"_$Nout"*"_$qli")
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
		    CSV.write(location*"results_jobs/benchmarks_HEUR/finalGRASP/$expname"*"/final_sols/sol_$seed"*"_$N"*"_$Nout"*"_$qli"*".csv", d)
		end
	    
		this_benchmark=DataFrame(Seed= [seed],N= [N],Nout= [Nout],qli= [qli],HeurCost= [ ceil(Int, cost)])
		newbenchmark=append!(newbenchmark,this_benchmark)
	end
    end
    return newbenchmark
end

#location = "D:/DTU-Courses/DTU-Thesis/berth_allocation/"
location="/zhome/c3/6/164957/code_git/"

# The tactic types :
type1="time" 
type2="cost" 
type3="random"

# The parameters of the experiment :

# The experience name :
expname="GRASPonly4local21_01_2"

# The window size for the visits to look at :
window=0.2

# The prop of ships to remove for the reconstruction :
proptoremove=0.01

# Should we push at constraints position :
pushatconstraint=true

# Look for constrained at the beginning :
lookforconstraint=false

# The max proportion of boats to remove for the local search :
alphaboat=0.2
alpharandom=0.1

# When do we start from zero the parameters (every n iterations) :
restartparams=50

# Which sol to take from the heuristic for the local search :
windowlocalsearch=0.1

# One boat tactic :
oneboat="cost"
onboatvec = [0.001, 0.1, 0.2, 0.3]

# All boat tactic :
allboat="all"
allboatvec = [0.001, 0.1, 0.3, 0.5]

# Reversed all boat tactic :
reversedallboat="all"
reversedallboatvec = [0.001, 0.1, 0.3, 0.5, 0.7]


# Local search tactics :
localsearch="all"
localsearchone="cost"
localsearchall="cost"

# Rate constrained :
alpharateconstrained = [0.2,0.4,0.6,0.8]

# Prop to remove :
alphapropremove = [0.001,0.05,0.15]

# Number of non improvement greedy algo :
greedymaxnoimprove=120000

# Make the heuristic without reconstruct until :
until=1000000

# Dont focus on removal without recontrusct improvement until :
focusremoveuntil=1000000

# Number of remove in a row :
nbfocusremove = 1000000

# Needed rate improvement to accept the reconstruction or the pathrelinking (after local) :
rateimprovereconstruct = 0.0015


# Make path relinking instead of ALNS :
pathrelinking="no"
# Max time for the relinking
maxtimerelinking=2
# The length of the elite set :
lengthelite=8
removepathrelinking= 0.2

reversed="no"
# All the parameters :
paramfixed = FixedParameters(oneboat, onboatvec, reversed, allboat, allboatvec, reversedallboat, reversedallboatvec, localsearch, localsearchone, localsearchall, 
alpharandom, alphaboat, 
alpharateconstrained, alphapropremove, 
window, pushatconstraint, lookforconstraint,
restartparams, greedymaxnoimprove, until, focusremoveuntil, nbfocusremove, 
rateimprovereconstruct, windowlocalsearch, 
pathrelinking, maxtimerelinking, lengthelite, removepathrelinking)

# Maximum time for the local search :
time_local=2
#time_local= parse(Int64,ARGS[1])

# Maximum time for the heuristic :
max_time_heur=30

# maximum time for the experiment :
max_time=2400
#max_time = parse(Int64,ARGS[2])
#max_time=100

# the temperature parameter :
temperature=0.93

# look for a specific seed
seedchosen = parse(Int64,ARGS[3])
seedchosen = 2
Nchosen=parse(Int64,ARGS[4])
Nchosen=30
#Noutchosen=parse(Int64,ARGS[5])
Noutchosen=5
#qlichosen=parse(Int64,ARGS[6])
qlichosen=10

allparam = initializeParam(paramfixed)
#list_paramvisit = allparam.Alpha.CostOneShip
#list_paramconstrained = allparam.Alpha.RateConstrained
makeExpText(temperature, paramfixed, time_local, max_time_heur, max_time, expname, location)
newbenchmark = makeSolHeur(paramfixed, temperature, time_local, max_time_heur, max_time, expname, location, seedchosen, Nchosen, Noutchosen, qlichosen)
CSV.write(location*"results_jobs/benchmarks_HEUR/finalGRASP/$expname"*"/NLarge_playground_test_$seedchosen"*"n_$Nchosen"*"_$Noutchosen_$qli.csv", newbenchmark)
newbenchmark


