import matplotlib.pyplot as plt
from getdict_results import *

def get_cmap(n, name='hsv'):
    return plt.cm.get_cmap(name, n)

def get_max_nested_dict_x(results,length, port):
    max=0
    for boat in results.keys():
        for visit in results["x"][boat]:
            new_length = results[boat][visit]+length[boat]
            if results[boat][visit]>max:
                max=results[boat][visit]
    return max

def getgraphs_port(algo_folder, exp, seed,Nin,Nout,qli):
    results = getsolfromfileHEUR(algo_folder=algo_folder, exp=exp, seed=seed,Nin=Nin,Nout=Nout,qli=qli)
    resultsCPLEX = getsolfromfileCPLEX(seed=seed,Nin=Nin,Nout=Nout,qli=qli)
    color =get_cmap(Nin)
    fig, axs = plt.subplots(3,2,figsize=(15, 15),sharex=True, sharey=True)
    max_time=0
    for port in [1,2,3]:
        max_pos=0
        ## Make the plot for each port 
        for n in range(0,len(results['inst'])):
            for c in range(0,len(results['inst'][n])):
                if results['inst'][n][c]==port:
                    if n+1<=Nin:
                        #print(n+1)
                        #print('Time :')
                        #print(results['y'][n+1][c+1])
                        #print(results['y'][n+1][c+1]+results['hand'][n+1][c+1])
                        #print('Postions :')
                        #print(results['x'][n+1][c+1])
                        #print(results['x'][n+1][c+1]+boatlen)
                        #print('#####')
                        boatlen = results['length'][n+1]
                        if results['y'][n+1][c+1]+results['hand'][n+1][c+1]>max_time:
                            max_time = results['y'][n+1][c+1]+results['hand'][n+1][c+1]
                        if results['x'][n+1][c+1]+boatlen>max_pos:
                            max_pos = results['x'][n+1][c+1]+boatlen
                            axs[port-1,0].set_ylim(0, max_pos)
                        rectangle = plt.Rectangle((results['y'][n+1][c+1],results['x'][n+1][c+1]), results['hand'][n+1][c+1],boatlen, fc=color(n),ec="black")
                        axs[port-1,0].add_patch(rectangle)
                        line = plt.Line2D((results['calls'][n+1][c+1][0], results['calls'][n+1][c+1][0]), (0, results['x'][n+1][c+1]+boatlen),color=color(n), lw=1.5)
                        axs[port-1,0].add_line(line)
                        centerx=results['y'][n+1][c+1]+(results['hand'][n+1][c+1]/20)
                        centery=results['x'][n+1][c+1]+(boatlen/2)
                        axs[port-1,0].text(centerx, centery,"b{}_v{}".format(n+1,c+1), color="black", fontsize=8)

                        if resultsCPLEX['y'][n+1][c+1]+resultsCPLEX['hand'][n+1][c+1]>max_time:
                            max_time = resultsCPLEX['y'][n+1][c+1]+resultsCPLEX['hand'][n+1][c+1]
                        if resultsCPLEX['x'][n+1][c+1]+boatlen>max_pos:
                            max_pos = resultsCPLEX['x'][n+1][c+1]+boatlen
                            axs[port-1,1].set_ylim(0, max_pos)
                        rectangle = plt.Rectangle((resultsCPLEX['y'][n+1][c+1],resultsCPLEX['x'][n+1][c+1]), resultsCPLEX['hand'][n+1][c+1],boatlen, fc=color(n),ec="black")
                        axs[port-1,1].add_patch(rectangle)
                        line = plt.Line2D((resultsCPLEX['calls'][n+1][c+1][0], resultsCPLEX['calls'][n+1][c+1][0]), (0, resultsCPLEX['x'][n+1][c+1]+boatlen),color=color(n), lw=1.5)
                        axs[port-1,1].add_line(line)
                        centerx=resultsCPLEX['y'][n+1][c+1]+(resultsCPLEX['hand'][n+1][c+1]/20)
                        centery=resultsCPLEX['x'][n+1][c+1]+(boatlen/2)
                        axs[port-1,1].text(centerx, centery,"b{}_v{}".format(n+1,c+1), color="black", fontsize=8)
                    else:
                        boatlen = results['length'][n+1]
                        if results['y'][n+1][c+1]+results['hand'][n+1][c+1]>max_time:
                            max_time = results['y'][n+1][c+1]+results['hand'][n+1][c+1]
                        if results['x'][n+1][c+1]+boatlen>max_pos:
                            max_pos = results['x'][n+1][c+1]+boatlen
                            axs[port-1,0].set_ylim(0, max_pos)
                        rectangle = plt.Rectangle((results['y'][n+1][c+1],results['x'][n+1][c+1]), results['hand'][n+1][c+1],boatlen, fc="black",ec="black")
                        axs[port-1,0].add_patch(rectangle)
                        centerx=results['y'][n+1][c+1]+(results['hand'][n+1][c+1]/20)
                        centery=results['x'][n+1][c+1]+(boatlen/2)
                        axs[port-1,0].text(centerx, centery,"b{}_v{}".format(n+1,c+1), color="white",fontsize=8)

                        if resultsCPLEX['y'][n+1][c+1]+resultsCPLEX['hand'][n+1][c+1]>max_time:
                            max_time = resultsCPLEX['y'][n+1][c+1]+resultsCPLEX['hand'][n+1][c+1]
                        if resultsCPLEX['x'][n+1][c+1]+boatlen>max_pos:
                            max_pos = resultsCPLEX['x'][n+1][c+1]+boatlen
                            axs[port-1,1].set_ylim(0, max_pos)
                        rectangle = plt.Rectangle((resultsCPLEX['y'][n+1][c+1],resultsCPLEX['x'][n+1][c+1]), resultsCPLEX['hand'][n+1][c+1],boatlen, fc="black",ec="black")
                        axs[port-1,1].add_patch(rectangle)
                        centerx=resultsCPLEX['y'][n+1][c+1]+(resultsCPLEX['hand'][n+1][c+1]/20)
                        centery=resultsCPLEX['x'][n+1][c+1]+(boatlen/2)
                        axs[port-1,1].text(centerx, centery,"b{}_v{}".format(n+1,c+1), color="white",fontsize=8)
    axs[0,0].set_title("Heuristics, cost : {} $".format(int(results["objectif"])))  
    axs[0,1].set_title("CPLEX, cost : {} $".format(int(resultsCPLEX["objectif"])))
    for port in [1,2,3]:
        axs[port-1,0].set_ylabel("Port {}".format(port))
        for i in [1,2]:
            axs[port-1, i-1].set_xlim(-1, max_time)  
    fig.tight_layout() 
    plt.show()         

