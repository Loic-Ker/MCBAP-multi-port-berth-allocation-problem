import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ast
from os import listdir
from os.path import isfile, join
import pandas as pd

def getsolfromfileCPLEX(Nin,Nout,seed,qli, type, time_limit):
    file=''
    if type=='small':
        if time_limit==5:
            file = "D:/DTU-Courses/DTU-Thesis/berth_allocation/MCBAP-multi-port-berth-allocation-problem/results_jobs/benchmarks_CPLEX/sols_5min/CPLEX_sol_{}_{}_{}_{}.csv".format(seed,Nin,Nout,qli)
        else:
            file = "D:/DTU-Courses/DTU-Thesis/berth_allocation/MCBAP-multi-port-berth-allocation-problem/results_jobs/benchmarks_CPLEX/sols_10min/CPLEX_sol_{}_{}_{}_{}.csv".format(seed,Nin,Nout,qli)
    if type == 'large':
        if time_limit==5:
            file = "D:/DTU-Courses/DTU-Thesis/berth_allocation/results_jobs/benchmarks_CPLEX/sols_2400s/CPLEX_sol_{}_{}_{}_{}.csv".format(seed,Nin,Nout,qli)
        else:
            file = "D:/DTU-Courses/DTU-Thesis/berth_allocation/results_jobs/benchmarks_CPLEX/sols_2400s/CPLEX_sol_{}_{}_{}_{}.csv".format(seed,Nin,Nout,qli)
    results_sol = pd.read_csv(file)
    dictresults=dict()

    listinst = list(results_sol[results_sol['first']=='inst']['second'])
    for i in range(0,len(listinst)):
        listinst[i]=ast.literal_eval(listinst[i])
    listinst=listinst[0]
    dictinst=dict()
    for n in range(0,len(listinst)):
        dictinst[n+1]= listinst[n]

    
    listlength=list(results_sol[results_sol['first']=='length_boats']['second'])
    for i in range(0,len(listlength)):
        listlength[i]=ast.literal_eval(listlength[i])
    listlength=listlength[0]
    dictlength=dict()
    for n in dictinst.keys():
        dictlength[n]=listlength[n-1]

    listcalls = ast.literal_eval(list(results_sol[results_sol['first']=='calls']['second'])[0])
    dictcalls=dict()
    for n in range(0,Nin):
        dictcalls[n+1]=dict()
        for c in range(0,len(listinst[n])):
            dictcalls[n+1][c+1]=listcalls[n][c]



    listx=list(results_sol[results_sol['first']=='x']['second'])[0].replace(' ','').split('\n')
    listy=list(results_sol[results_sol['first']=='y']['second'])[0].replace(' ','').split('\n')
    listh=list(results_sol[results_sol['first']=='hand']['second'])[0].replace(' ','').split('\n')
    for i in range(0,len(listx)):
        listx[i]=listx[i].split('=')
        listx[i][0]=ast.literal_eval(listx[i][0])
        listx[i][1]=float(listx[i][1])
        listy[i]=listy[i].split('=')
        listy[i][0]=ast.literal_eval(listy[i][0])
        listy[i][1]=float(listy[i][1])
        listh[i]=listh[i].split('=')
        listh[i][0]=ast.literal_eval(listh[i][0])
        listh[i][1]=float(listh[i][1])
    dictx=dict()
    dicty=dict()
    dicth=dict()
    for i in range(0,len(listh)):
        n=listh[i][0][0]
        if n not in dicth.keys():
            dicth[n]=dict()
        c=listh[i][0][1]
        dicth[n][c]=listh[i][1]
    for i in range(0,len(listx)):
        n=listx[i][0][0]
        if n not in dictx.keys():
            dictx[n]=dict()
        c=listx[i][0][1]
        dictx[n][c]=listx[i][1]
    for i in range(0,len(listy)):
        n=listy[i][0][0]
        if n not in dicty.keys():
            dicty[n]=dict()
        c=listy[i][0][1]
        dicty[n][c]=listy[i][1]
    dictresults['x']=dictx
    dictresults['y']=dicty
    dictresults['hand']=dicth
    dictresults['length']=dictlength
    dictresults['calls']=dictcalls
    dictresults['inst']=listinst
    dictresults['objectif']=float(list(results_sol[results_sol['first']=='objectif']['second'])[0])
    return dictresults

def getsolfromfileHEUR(Nin,Nout,seed,qli,algo_folder,exp):
    file = "D:/DTU-Courses/DTU-Thesis/berth_allocation/results_jobs/benchmarks_HEUR/{}/{}/final_sols/sol_{}_{}_{}_{}.csv".format(algo_folder,exp,seed,Nin,Nout,qli)
    results_sol = pd.read_csv(file)
    dictresults=dict()

    listinst = list(results_sol[results_sol['first']=='inst']['second'])
    for i in range(0,len(listinst)):
        listinst[i]=ast.literal_eval(listinst[i])
    listinst=listinst[0]
    dictinst=dict()
    for n in range(0,len(listinst)):
        dictinst[n+1]= listinst[n]

    
    listlength=list(results_sol[results_sol['first']=='length_boats']['second'])
    for i in range(0,len(listlength)):
        listlength[i]=ast.literal_eval(listlength[i])
    listlength=listlength[0]
    dictlength=dict()
    for n in dictinst.keys():
        dictlength[n]=listlength[n-1]

    listcalls = ast.literal_eval(list(results_sol[results_sol['first']=='calls']['second'])[0])
    dictcalls=dict()
    for n in range(0,Nin):
        dictcalls[n+1]=dict()
        for c in range(0,len(listinst[n])):
            dictcalls[n+1][c+1]=listcalls[n][c]



    listx=ast.literal_eval(list(results_sol[results_sol['first']=='x']['second'])[0].replace("Vector{Any}",''))
    listy=ast.literal_eval(list(results_sol[results_sol['first']=='y']['second'])[0].replace("Vector{Any}",''))
    listh=ast.literal_eval(list(results_sol[results_sol['first']=='hand']['second'])[0].replace("Vector{Any}",''))
    dictx=dict()
    dicty=dict()
    dicth=dict()
    for n in range(1,len(listh)+1):
        if n not in dicth.keys():
            dicth[n]=dict()
        for c in range(1,len(listh[n-1])+1):
            dicth[n][c]=listh[n-1][c-1]
    for n in range(1,len(listx)+1):
        if n not in dictx.keys():
            dictx[n]=dict()
        for c in range(1,len(listx[n-1])+1):
            dictx[n][c]=listx[n-1][c-1]
    for n in range(1,len(listy)+1):
        if n not in dicty.keys():
            dicty[n]=dict()
        for c in range(1,len(listy[n-1])+1):
            dicty[n][c]=listy[n-1][c-1]

    dictresults['x']=dictx
    dictresults['y']=dicty
    dictresults['hand']=dicth
    dictresults['length']=dictlength
    dictresults['calls']=dictcalls
    dictresults['inst']=listinst
    dictresults['objectif']=float(list(results_sol[results_sol['first']=='objectif']['second'])[0])
    return dictresults

def getspecificiterfromfileHEUR(Nin,Nout,seed,qli,algo_folder,exp,iternb):
    file = "D:/DTU-Courses/DTU-Thesis/berth_allocation/results_jobs/benchmarks_HEUR/{}/{}/iterations/sol_{}_{}_{}_{}/iter_{}.csv".format(algo_folder,exp,seed,Nin,Nout,qli,iternb)
    results_sol = pd.read_csv(file)
    dictresults=dict()

    listinst = list(results_sol[results_sol['first']=='inst']['second'])
    for i in range(0,len(listinst)):
        listinst[i]=ast.literal_eval(listinst[i])
    listinst=listinst[0]
    dictinst=dict()
    for n in range(0,len(listinst)):
        dictinst[n+1]= listinst[n]

    
    listlength=list(results_sol[results_sol['first']=='length_boats']['second'])
    for i in range(0,len(listlength)):
        listlength[i]=ast.literal_eval(listlength[i])
    listlength=listlength[0]
    dictlength=dict()
    for n in dictinst.keys():
        dictlength[n]=listlength[n-1]

    listcalls = ast.literal_eval(list(results_sol[results_sol['first']=='calls']['second'])[0])
    dictcalls=dict()
    for n in range(0,Nin):
        dictcalls[n+1]=dict()
        for c in range(0,len(listinst[n])):
            dictcalls[n+1][c+1]=listcalls[n][c]



    listx=ast.literal_eval(list(results_sol[results_sol['first']=='x']['second'])[0].replace("Vector{Any}",''))
    listy=ast.literal_eval(list(results_sol[results_sol['first']=='y']['second'])[0].replace("Vector{Any}",''))
    listh=ast.literal_eval(list(results_sol[results_sol['first']=='hand']['second'])[0].replace("Vector{Any}",''))
    dictx=dict()
    dicty=dict()
    dicth=dict()
    for n in range(1,len(listh)+1):
        if n not in dicth.keys():
            dicth[n]=dict()
        for c in range(1,len(listh[n-1])+1):
            dicth[n][c]=listh[n-1][c-1]
    for n in range(1,len(listx)+1):
        if n not in dictx.keys():
            dictx[n]=dict()
        for c in range(1,len(listx[n-1])+1):
            dictx[n][c]=listx[n-1][c-1]
    for n in range(1,len(listy)+1):
        if n not in dicty.keys():
            dicty[n]=dict()
        for c in range(1,len(listy[n-1])+1):
            dicty[n][c]=listy[n-1][c-1]

    dictresults['x']=dictx
    dictresults['y']=dicty
    dictresults['hand']=dicth
    dictresults['length']=dictlength
    dictresults['calls']=dictcalls
    dictresults['inst']=listinst
    dictresults['objectif']=float(list(results_sol[results_sol['first']=='objectif']['second'])[0])
    return dictresults

def getiterfromfileHEUR(Nin,Nout,seed,qli,algo_folder,exp):
    location = "D:/DTU-Courses/DTU-Thesis/berth_allocation/results_jobs/benchmarks_HEUR/{}/{}/iterations/sol_{}_{}_{}_{}".format(algo_folder,exp,seed,Nin,Nout,qli)
    list_iter = [f for f in listdir(location) if isfile(join(location, f))]
    nb_iter=0
    dict_all_iters=dict()
    for file in list_iter:
        nb_iter+=1
        results_sol = pd.read_csv(location+'/'+file)

        dictresults=dict()

        

        listinst = list(results_sol[results_sol['first']=='inst']['second'])
        for i in range(0,len(listinst)):
            listinst[i]=ast.literal_eval(listinst[i])
        listinst=listinst[0]
        dictinst=dict()
        for n in range(0,len(listinst)):
            dictinst[n+1]= listinst[n]

        
        listlength=list(results_sol[results_sol['first']=='length_boats']['second'])
        for i in range(0,len(listlength)):
            listlength[i]=ast.literal_eval(listlength[i])
        listlength=listlength[0]
        dictlength=dict()
        for n in dictinst.keys():
            dictlength[n]=listlength[n-1]

        listcalls = ast.literal_eval(list(results_sol[results_sol['first']=='calls']['second'])[0])
        dictcalls=dict()
        for n in range(0,Nin):
            dictcalls[n+1]=dict()
            for c in range(0,len(listinst[n])):
                dictcalls[n+1][c+1]=listcalls[n][c]



        listx=ast.literal_eval(list(results_sol[results_sol['first']=='x']['second'])[0].replace("Vector{Any}",''))
        listy=ast.literal_eval(list(results_sol[results_sol['first']=='y']['second'])[0].replace("Vector{Any}",''))
        listh=ast.literal_eval(list(results_sol[results_sol['first']=='hand']['second'])[0].replace("Vector{Any}",''))
        list_cost_visit=ast.literal_eval(list(results_sol[results_sol['first']=='cost_visit']['second'])[0].replace("Vector{Any}",''))
        list_delay_cost_visit=ast.literal_eval(list(results_sol[results_sol['first']=='delay_cost_visit']['second'])[0].replace("Vector{Any}",''))
        list_when=ast.literal_eval(list(results_sol[results_sol['first']=='when']['second'])[0].replace("Vector{Any}",''))
        list_tacticall_chosen=ast.literal_eval(list(results_sol[results_sol['first']=='tacticall_chosen']['second'])[0].replace("Vector{Any}",''))
        list_tacticboat_chosen=ast.literal_eval(list(results_sol[results_sol['first']=='tacticboat_chosen']['second'])[0].replace("Vector{Any}",''))
        list_penalty_visit=ast.literal_eval(list(results_sol[results_sol['first']=='penalty_visit']['second'])[0].replace("Vector{Any}",''))
        list_fuel_cost_visit=ast.literal_eval(list(results_sol[results_sol['first']=='fuel_cost_visit']['second'])[0].replace("Vector{Any}",''))
        list_handling_cost_visit=ast.literal_eval(list(results_sol[results_sol['first']=='handling_cost_visit']['second'])[0].replace("Vector{Any}",''))
        list_waiting_cost_visit=ast.literal_eval(list(results_sol[results_sol['first']=='waiting_cost_visit']['second'])[0].replace("Vector{Any}",''))
        
        dictx=dict()
        dicty=dict()
        dicth=dict()
        dict_cost_visit=dict()
        dict_delay_cost_visit=dict()
        dict_when=dict()
        dict_tacticall_chosen=dict()
        dict_tacticboat_chosen=dict()
        dict_penalty_visit=dict()
        dict_fuel_cost_visit=dict()
        dict_handling_cost_visit=dict()
        dict_waiting_cost_visit=dict()
        for n in range(1,len(listh)+1):
            if n not in dicth.keys():
                dicth[n]=dict()
            for c in range(1,len(listh[n-1])+1):
                dicth[n][c]=listh[n-1][c-1]
        for n in range(1,len(listx)+1):
            if n not in dictx.keys():
                dictx[n]=dict()
            for c in range(1,len(listx[n-1])+1):
                dictx[n][c]=listx[n-1][c-1]
        for n in range(1,len(listy)+1):
            if n not in dicty.keys():
                dicty[n]=dict()
            for c in range(1,len(listy[n-1])+1):
                dicty[n][c]=listy[n-1][c-1]
        for n in range(1,len(list_cost_visit)+1):
            if n not in dict_cost_visit.keys():
                dict_cost_visit[n]=dict()
            for c in range(1,len(list_cost_visit[n-1])+1):
                dict_cost_visit[n][c]=list_cost_visit[n-1][c-1]
        for n in range(1,len(list_delay_cost_visit)+1):
            if n not in dict_delay_cost_visit.keys():
                dict_delay_cost_visit[n]=dict()
            for c in range(1,len(list_delay_cost_visit[n-1])+1):
                dict_delay_cost_visit[n][c]=list_delay_cost_visit[n-1][c-1]
        for n in range(1,len(list_when)+1):
            if n not in dict_when.keys():
                dict_when[n]=dict()
            for c in range(1,len(list_when[n-1])+1):
                dict_when[n][c]=list_when[n-1][c-1]
        for n in range(1,len(list_tacticall_chosen)+1):
            if n not in dict_tacticall_chosen.keys():
                dict_tacticall_chosen[n]=dict()
            for c in range(1,len(list_tacticall_chosen[n-1])+1):
                dict_tacticall_chosen[n][c]=list_tacticall_chosen[n-1][c-1]
        for n in range(1,len(list_tacticboat_chosen)+1):
            if n not in dict_tacticboat_chosen.keys():
                dict_tacticboat_chosen[n]=dict()
            for c in range(1,len(list_tacticboat_chosen[n-1])+1):
                dict_tacticboat_chosen[n][c]=list_tacticboat_chosen[n-1][c-1]
        for n in range(1,len(list_penalty_visit)+1):
            if n not in dict_penalty_visit.keys():
                dict_penalty_visit[n]=dict()
            for c in range(1,len(list_penalty_visit[n-1])+1):
                dict_penalty_visit[n][c]=list_penalty_visit[n-1][c-1]
        for n in range(1,len(list_fuel_cost_visit)+1):
            if n not in dict_fuel_cost_visit.keys():
                dict_fuel_cost_visit[n]=dict()
            for c in range(1,len(list_fuel_cost_visit[n-1])+1):
                dict_fuel_cost_visit[n][c]=list_fuel_cost_visit[n-1][c-1]
        for n in range(1,len(list_handling_cost_visit)+1):
            if n not in dict_handling_cost_visit.keys():
                dict_handling_cost_visit[n]=dict()
            for c in range(1,len(list_handling_cost_visit[n-1])+1):
                dict_handling_cost_visit[n][c]=list_handling_cost_visit[n-1][c-1]
        for n in range(1,len(list_waiting_cost_visit)+1):
            if n not in dict_waiting_cost_visit.keys():
                dict_waiting_cost_visit[n]=dict()
            for c in range(1,len(list_waiting_cost_visit[n-1])+1):
                dict_waiting_cost_visit[n][c]=list_waiting_cost_visit[n-1][c-1]

        dictresults['x']=dictx
        dictresults['y']=dicty
        dictresults['hand']=dicth
        dictresults['cost_visit']=dict_cost_visit
        dictresults['delay_cost_visit']=dict_delay_cost_visit
        dictresults['when']=dict_when
        dictresults['tacticall_chosen']=dict_tacticall_chosen
        dictresults['tacticboat_chosen']=dict_tacticboat_chosen
        dictresults['penalty_visit']=dict_penalty_visit
        dictresults['fuel_cost_visit']=dict_fuel_cost_visit
        dictresults['handling_cost_visit']=dict_handling_cost_visit
        dictresults['waiting_cost_visit']=dict_waiting_cost_visit
        dictresults['length']=dictlength
        dictresults['calls']=dictcalls
        dictresults['inst']=listinst
        dictresults['objectif']=float(list(results_sol[results_sol['first']=='objectif']['second'])[0])

        dictresults["cost_solheur"] = float(results_sol[results_sol['first']=="cost_solheur"]['second'])
        dictresults["delay_cost_solheur"] = float(results_sol[results_sol['first']=="delay_cost_solheur"]['second'])
        dictresults["waiting_cost_solheur"] = float(results_sol[results_sol['first']=="waiting_cost_solheur"]['second'])
        dictresults["penalty_solheur"] = float(results_sol[results_sol['first']=="penalty_solheur"]['second'])
        dictresults["handling_cost_solheur"] = float(results_sol[results_sol['first']=="handling_cost_solheur"]['second'])
        dictresults["fuel_cost_solheur"] = float(results_sol[results_sol['first']=="fuel_cost_solheur"]['second'])

        dictresults["cost_sollocal"] = float(results_sol[results_sol['first']=="cost_sollocal"]['second'])
        dictresults["delay_cost_sollocal"] = float(results_sol[results_sol['first']=="delay_cost_sollocal"]['second'])
        dictresults["waiting_cost_sollocal"] = float(results_sol[results_sol['first']=="waiting_cost_sollocal"]['second'])
        dictresults["penalty_sollocal"] = float(results_sol[results_sol['first']=="penalty_sollocal"]['second'])
        dictresults["handling_cost_sollocal"] = float(results_sol[results_sol['first']=="handling_cost_sollocal"]['second'])
        dictresults["fuel_cost_sollocal"] =float(results_sol[results_sol['first']=="fuel_cost_sollocal"]['second'])

        dictresults["timeheur"]= float(results_sol[results_sol['first']=="timeheur"]['second'])
        dictresults["timelocal"]= float(results_sol[results_sol['first']=="timelocal"]['second'])
        dictresults["oneboatdistance"] = ast.literal_eval(list(results_sol[results_sol['first']=="oneboatdistance"]['second'])[0])
        dictresults["oneboatcost"] = ast.literal_eval(list(results_sol[results_sol['first']=="oneboatcost"]['second'])[0])
        dictresults["oneboattime"] = ast.literal_eval(list(results_sol[results_sol['first']=="oneboattime"]['second'])[0])
        dictresults["allboatsdist"] = ast.literal_eval(list(results_sol[results_sol['first']=="allboatsdist"]['second'])[0])
        dictresults["allboatstime"] = ast.literal_eval(list(results_sol[results_sol['first']=="allboatstime"]['second'])[0])
        dictresults["allboatscost"] = ast.literal_eval(list(results_sol[results_sol['first']=="allboatscost"]['second'])[0])
        dictresults["rateconstrained"] = ast.literal_eval(list(results_sol[results_sol['first']=="rateconstrained"]['second'])[0])
        
        dictresults["proba_tacticboat"] = ast.literal_eval(list(results_sol[results_sol['first']=="proba_tacticboat"]['second'])[0])
        dictresults["proba_tacticall"] = ast.literal_eval(list(results_sol[results_sol['first']=="proba_tacticall"]['second'])[0])
        dictresults["proba_tacticlocalsearch"] = ast.literal_eval(list(results_sol[results_sol['first']=="proba_tacticlocalsearch"]['second'])[0])
        
        dictresults["failed"] = float(results_sol[results_sol['first']=="failed"]['second'])
        dictresults["better"] = float(results_sol[results_sol['first']=="better"]['second'])
        
        dictresults["reconstruct"] = float(list(results_sol[results_sol['first']=="reconstruct"]['second'])[0])
        dictresults["usedCPLEX"] = float(list(results_sol[results_sol['first']=="usedCPLEX"]['second'])[0])

        dict_all_iters[nb_iter]=dictresults

    return dict_all_iters


def getiterfromfileHEURSoft(Nin,Nout,seed,qli,algo_folder,exp):
    location = "D:/DTU-Courses/DTU-Thesis/berth_allocation/results_jobs/benchmarks_HEUR/{}/{}/iterations/sol_{}_{}_{}_{}".format(algo_folder,exp,seed,Nin,Nout,qli)
    list_iter = [f for f in listdir(location) if isfile(join(location, f))]
    nb_iter=0
    dict_all_iters=dict()
    for file in list_iter:
        nb_iter+=1
        results_sol = pd.read_csv(location+'/'+file)

        dictresults=dict()

        

        listinst = list(results_sol[results_sol['first']=='inst']['second'])
        for i in range(0,len(listinst)):
            listinst[i]=ast.literal_eval(listinst[i])
        listinst=listinst[0]
        dictinst=dict()
        for n in range(0,len(listinst)):
            dictinst[n+1]= listinst[n]

        
        listlength=list(results_sol[results_sol['first']=='length_boats']['second'])
        for i in range(0,len(listlength)):
            listlength[i]=ast.literal_eval(listlength[i])
        listlength=listlength[0]
        dictlength=dict()
        for n in dictinst.keys():
            dictlength[n]=listlength[n-1]

        listcalls = ast.literal_eval(list(results_sol[results_sol['first']=='calls']['second'])[0])
        dictcalls=dict()
        for n in range(0,Nin):
            dictcalls[n+1]=dict()
            for c in range(0,len(listinst[n])):
                dictcalls[n+1][c+1]=listcalls[n][c]

        dictresults['length']=dictlength
        dictresults['calls']=dictcalls
        dictresults['inst']=listinst
        dictresults['objectif']=float(list(results_sol[results_sol['first']=='objectif']['second'])[0])

        dictresults["chosen_tacticoneboat"] = results_sol[results_sol['first']=="chosen_tacticoneboat"]['second']
        dictresults["chosen_reversed"] = results_sol[results_sol['first']=="chosen_reversed"]['second']
        dictresults["chosen_tacticallboats"] = results_sol[results_sol['first']=="chosen_tacticallboats"]['second']
        dictresults["chosen_reversedtacticallboats"] = results_sol[results_sol['first']=="chosen_reversedtacticallboats"]['second']
        dictresults["chosen_tacticlocalsearch"] = results_sol[results_sol['first']=="chosen_tacticlocalsearch"]['second']
        dictresults["chosen_indexoneship"] = results_sol[results_sol['first']=="chosen_indexoneship"]['second']
        dictresults["chosen_indexallship"] = results_sol[results_sol['first']=="chosen_indexallship"]['second']
        dictresults["chosen_indexreversedallship"] = results_sol[results_sol['first']=="chosen_indexreversedallship"]['second']
        dictresults["chosen_indexrateconstrained"] = results_sol[results_sol['first']=="chosen_indexrateconstrained"]['second']
        dictresults["chosen_indexproptoremove"] = results_sol[results_sol['first']=="chosen_indexproptoremove"]['second']
    
        dictresults["cost_solheur"] = float(results_sol[results_sol['first']=="cost_solheur"]['second'])
        dictresults["delay_cost_solheur"] = float(results_sol[results_sol['first']=="delay_cost_solheur"]['second'])
        dictresults["waiting_cost_solheur"] = float(results_sol[results_sol['first']=="waiting_cost_solheur"]['second'])
        dictresults["penalty_solheur"] = float(results_sol[results_sol['first']=="penalty_solheur"]['second'])
        dictresults["handling_cost_solheur"] = float(results_sol[results_sol['first']=="handling_cost_solheur"]['second'])
        dictresults["fuel_cost_solheur"] = float(results_sol[results_sol['first']=="fuel_cost_solheur"]['second'])

        dictresults["cost_sollocal"] = float(results_sol[results_sol['first']=="cost_sollocal"]['second'])
        dictresults["delay_cost_sollocal"] = float(results_sol[results_sol['first']=="delay_cost_sollocal"]['second'])
        dictresults["waiting_cost_sollocal"] = float(results_sol[results_sol['first']=="waiting_cost_sollocal"]['second'])
        dictresults["penalty_sollocal"] = float(results_sol[results_sol['first']=="penalty_sollocal"]['second'])
        dictresults["handling_cost_sollocal"] = float(results_sol[results_sol['first']=="handling_cost_sollocal"]['second'])
        dictresults["fuel_cost_sollocal"] =float(results_sol[results_sol['first']=="fuel_cost_sollocal"]['second'])

        dictresults["timeheur"]= float(results_sol[results_sol['first']=="timeheur"]['second'])
        dictresults["timelocal"]= float(results_sol[results_sol['first']=="timelocal"]['second'])
        dictresults["oneboatdistance"] = ast.literal_eval(list(results_sol[results_sol['first']=="oneboatdistance"]['second'])[0])
        dictresults["oneboatcost"] = ast.literal_eval(list(results_sol[results_sol['first']=="oneboatcost"]['second'])[0])
        dictresults["oneboattime"] = ast.literal_eval(list(results_sol[results_sol['first']=="oneboattime"]['second'])[0])
        dictresults["allboatsdist"] = ast.literal_eval(list(results_sol[results_sol['first']=="allboatsdist"]['second'])[0])
        dictresults["allboatstime"] = ast.literal_eval(list(results_sol[results_sol['first']=="allboatstime"]['second'])[0])
        dictresults["allboatscost"] = ast.literal_eval(list(results_sol[results_sol['first']=="allboatscost"]['second'])[0])
        dictresults["rateconstrained"] = ast.literal_eval(list(results_sol[results_sol['first']=="rateconstrained"]['second'])[0])
        
        dictresults["proba_tacticboat"] = ast.literal_eval(list(results_sol[results_sol['first']=="proba_tacticboat"]['second'])[0])
        dictresults["proba_tacticall"] = ast.literal_eval(list(results_sol[results_sol['first']=="proba_tacticall"]['second'])[0])
        dictresults["proba_tacticlocalsearch"] = ast.literal_eval(list(results_sol[results_sol['first']=="proba_tacticlocalsearch"]['second'])[0])
        
        dictresults["failed"] = float(results_sol[results_sol['first']=="failed"]['second'])
        
        try:
            dictresults["better"] = float(results_sol[results_sol['first']=="better"]['second'])
        except:
            dictresults["better"] = "None"
        
        try:
            dictresults["reconstruct"] = float(list(results_sol[results_sol['first']=="reconstruct"]['second'])[0])
        except:
            dictresults["reconstruct"] = "None"
        
        try:
            dictresults["pathrelinking"] = ast.literal_eval(list(results_sol[results_sol['first']=="pathrelinking"]['second'])[0])
        except:
            dictresults["pathrelinking"] = "None"
            
        try:
            dictresults["usedLocalSearch"] = ast.literal_eval(list(results_sol[results_sol['first']=="usedLocalSearch"]['second'])[0])
        except:
            dictresults["usedLocalSearch"] = "None"
        
        try:
            dictresults["average_cost_elite"] = ast.literal_eval(list(results_sol[results_sol['first']=="average_cost_elite"]['second'])[0])
        except:
            dictresults["average_cost_elite"] = "None"
        
        try:
            dictresults["average_dist_elite"] = ast.literal_eval(list(results_sol[results_sol['first']=="average_dist_elite"]['second'])[0])
        except:
            dictresults["average_dist_elite"] = "None"
            
        try:
            dictresults["pushimprove"] = ast.literal_eval(list(results_sol[results_sol['first']=="pushimprove"]['second'])[0])
        except:
            dictresults["pushimprove"] = "None"

        dict_all_iters[nb_iter]=dictresults

    return dict_all_iters


def make_datasetiter(algo_folder, exp, seed, Nin,Nout,qli):
    #results = getiterfromfileHEURSoft(algo_folder=algo_folder, exp=exp, seed=seed,Nin=Nin,Nout=Nout,qli=qli)
    results = getiterfromfileHEUR(algo_folder=algo_folder, exp=exp, seed=seed,Nin=Nin,Nout=Nout,qli=qli)
    dataset = pd.DataFrame(columns=('iter','n','c','port',
            'cost_visit',
            'delay_cost_visit',
            'when',
            'tacticall_chosen',
            'tacticboat_chosen',
            'penalty_visit',
            'fuel_cost_visit',
            'handling_cost_visit',
            'waiting_cost_visit',
            'length',
            'calls',
            'inst',
            'objectif',
            "cost_solheur",
            "delay_cost_solheur",
            "waiting_cost_solheur",
            "penalty_solheur",
            "handling_cost_solheur",
            "fuel_cost_solheur",
            "cost_sollocal",
            "delay_cost_sollocal",
            "waiting_cost_sollocal",
            "penalty_sollocal",
            "handling_cost_sollocal",
            "fuel_cost_sollocal",
            "timeheur",
            "timelocal",
            "oneboatdistance",
            "oneboatcost",
            "oneboattime",
            "allboatsdist",
            "allboatstime",
            "allboatscost",
            "rateconstained",
            "failed",
            "reconstruct",
            "usedCPLEX"))

    for i in results.keys():
        for n in range(1,Nin+1):
            for c in range(1,len(results[i]['inst'][n-1])+1):
                dataset_row = {'iter':i,'n':n,'c':c,'port':results[i]['inst'][n-1][c-1],
'cost_visit':results[i]['cost_visit'][n][c],
'delay_cost_visit':results[i]['delay_cost_visit'][n][c],
'when':results[i]['when'][n][c],
'tacticall_chosen':results[i]['tacticall_chosen'][n][c],
'tacticboat_chosen':results[i]['tacticboat_chosen'][n][c],
'penalty_visit':results[i]['penalty_visit'][n][c],
'fuel_cost_visit':results[i]['fuel_cost_visit'][n][c],
'handling_cost_visit':results[i]['handling_cost_visit'][n][c],
'waiting_cost_visit':results[i]['waiting_cost_visit'][n][c],
'length':results[i]['length'][n],
'calls':results[i]['calls'][n],
'objectif':results[i]['objectif'],
"cost_solheur":results[i]['cost_solheur'],
"delay_cost_solheur":results[i]['delay_cost_solheur'],
"waiting_cost_solheur":results[i]['waiting_cost_solheur'],
"penalty_solheur":results[i]['penalty_solheur'],
"handling_cost_solheur":results[i]['handling_cost_solheur'],
"fuel_cost_solheur":results[i]['fuel_cost_solheur'],
"cost_sollocal":results[i]['cost_sollocal'],
"delay_cost_sollocal":results[i]['delay_cost_sollocal'],
"waiting_cost_sollocal":results[i]['waiting_cost_sollocal'],
"penalty_sollocal":results[i]['penalty_sollocal'],
"handling_cost_sollocal":results[i]['handling_cost_sollocal'],
"fuel_cost_sollocal":results[i]['fuel_cost_sollocal'],
"timeheur":results[i]['timeheur'],
"timelocal":results[i]['timelocal'],
"oneboatdistance":results[i]['oneboatdistance'],
"oneboatcost":results[i]['oneboatcost'],
"oneboattime":results[i]['oneboattime'],
"allboatsdist":results[i]['allboatsdist'],
"allboatstime":results[i]['allboatstime'],
"allboatscost":results[i]['allboatscost'],
"rateconstained":results[i]['rateconstrained'],
"failed":results[i]['failed'],
'reconstruct':results[i]['reconstruct'],
'usedCPLEX':results[i]['usedCPLEX']}
                dataset = dataset.append(dataset_row, ignore_index=True)
    return dataset


def make_datasetiterSoft(algo_folder, exp, seed, Nin,Nout,qli):
    results = getiterfromfileHEURSoft(algo_folder=algo_folder, exp=exp, seed=seed,Nin=Nin,Nout=Nout,qli=qli)
    #results = getiterfromfileHEUR(algo_folder=algo_folder, exp=exp, seed=seed,Nin=Nin,Nout=Nout,qli=qli)
    dataset = pd.DataFrame(columns=('iter',
            'inst',
            'objectif',
            "chosen_tacticoneboat",
            "chosen_reversed",
            "chosen_tacticallboats",
            "chosen_reversedtacticallboats",
            "chosen_tacticlocalsearch",
            "chosen_indexoneship",
            "chosen_indexallship",
            "chosen_indexreversedallship",
            "chosen_indexrateconstrained",
            "chosen_indexproptoremove",
            "cost_solheur",
            "delay_cost_solheur",
            "waiting_cost_solheur",
            "penalty_solheur",
            "handling_cost_solheur",
            "fuel_cost_solheur",
            "cost_sollocal",
            "delay_cost_sollocal",
            "waiting_cost_sollocal",
            "penalty_sollocal",
            "handling_cost_sollocal",
            "fuel_cost_sollocal",
            "timeheur",
            "timelocal",
            "proba_tacticboat",
            "proba_tacticall",
            "proba_tacticlocalsearch",
            "oneboatdistance",
            "oneboatcost",
            "oneboattime",
            "allboatsdist",
            "allboatstime",
            "allboatscost",
            "rateconstained",
            "failed",
            "better",
            "reconstruct",
            "pathrelinking",
            "usedLocalSearch",
            "average_cost_elite",
            "average_dist_elite",
            "pushimprove"
            ))

    for i in results.keys():
        dataset_row = {'iter':i,
'objectif':results[i]['objectif'],
"chosen_tacticoneboat":results[i]['chosen_tacticoneboat'],
"chosen_reversed":results[i]['chosen_reversed'],
"chosen_tacticallboats":results[i]['chosen_tacticallboats'],
"chosen_reversedtacticallboats":results[i]['chosen_reversedtacticallboats'],
"chosen_tacticlocalsearch":results[i]['chosen_tacticlocalsearch'],
"chosen_indexoneship":results[i]['chosen_indexoneship'],
"chosen_indexallship":results[i]['chosen_indexallship'],
"chosen_indexreversedallship":results[i]['chosen_indexreversedallship'],
"chosen_indexrateconstrained":results[i]['chosen_indexrateconstrained'],
"chosen_indexproptoremove":results[i]['chosen_indexproptoremove'],
"cost_solheur":results[i]['cost_solheur'],
"delay_cost_solheur":results[i]['delay_cost_solheur'],
"waiting_cost_solheur":results[i]['waiting_cost_solheur'],
"penalty_solheur":results[i]['penalty_solheur'],
"handling_cost_solheur":results[i]['handling_cost_solheur'],
"fuel_cost_solheur":results[i]['fuel_cost_solheur'],
"cost_sollocal":results[i]['cost_sollocal'],
"delay_cost_sollocal":results[i]['delay_cost_sollocal'],
"waiting_cost_sollocal":results[i]['waiting_cost_sollocal'],
"penalty_sollocal":results[i]['penalty_sollocal'],
"handling_cost_sollocal":results[i]['handling_cost_sollocal'],
"fuel_cost_sollocal":results[i]['fuel_cost_sollocal'],
"timeheur":results[i]['timeheur'],
"timelocal":results[i]['timelocal'],
"proba_tacticboat":results[i]['proba_tacticboat'],
"proba_tacticall":results[i]['proba_tacticall'],
"proba_tacticlocalsearch":results[i]['proba_tacticlocalsearch'],
"oneboatdistance":results[i]['oneboatdistance'],
"oneboatcost":results[i]['oneboatcost'],
"oneboattime":results[i]['oneboattime'],
"allboatsdist":results[i]['allboatsdist'],
"allboatstime":results[i]['allboatstime'],
"allboatscost":results[i]['allboatscost'],
"rateconstained":results[i]['rateconstrained'],
"failed":results[i]['failed'],
"better":results[i]['better'],
'reconstruct':results[i]['reconstruct'],
"pathrelinking":results[i]['pathrelinking'],
"usedLocalSearch":results[i]['usedLocalSearch'],
"average_cost_elite":results[i]['average_cost_elite'],
"average_dist_elite":results[i]['average_dist_elite'],
"pushimprove":results[i]['pushimprove']
}
        
        dataset = dataset.append(dataset_row, ignore_index=True)
    return dataset