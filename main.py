import heapq
from search.algorithms import State
from search.map import Map
import getopt
import sys

def main():
    """
    Function for testing your implementation. Run it with a -help option to see the options available. 
    """
    optlist, _ = getopt.getopt(sys.argv[1:], 'h:m:r:', ['testinstances', 'plots', 'help'])

    plots = False
    for o, _ in optlist:
        if o in ("-help"):
            print("Examples of Usage:")
            print("Solve set of test instances: main.py --testinstances")
            print("Solve set of test instances and generate plots: main.py --testinstances --plots")
            exit()
        elif o in ("--plots"):
            plots = True
    test_instances = "test-instances/testinstances.txt"
    gridded_map = Map("dao-map/brc000d.map")
    
    nodes_expanded_biastar = []   
    nodes_expanded_astar = []   
    nodes_expanded_mm = []
    
    start_states = []
    goal_states = []
    solution_costs = []
       
    file = open(test_instances, "r")
    for instance_string in file:
        list_instance = instance_string.split(",")
        start_states.append(State(int(list_instance[0]), int(list_instance[1])))
        goal_states.append(State(int(list_instance[2]), int(list_instance[3])))
        
        solution_costs.append(float(list_instance[4]))
    file.close()
        
    for i in range(0, len(start_states)):   

        start = start_states[i]
        goal = goal_states[i]
    
        cost, expanded_astar = astar(start, goal, gridded_map)
        nodes_expanded_astar.append(expanded_astar)

        if cost != solution_costs[i]:
            print("There is a mismatch in the solution cost found by A* and what was expected for the problem:")
            print("Start state: ", start)
            print("Goal state: ", goal)
            print("Solution cost encountered: ", cost)
            print("Solution cost expected: ", solution_costs[i])
            print()
        
        cost, expanded_mm = mm(start, goal, gridded_map)
        nodes_expanded_mm.append(expanded_mm)
        
        if cost != solution_costs[i]:
            print("There is a mismatch in the solution cost found by MM and what was expected for the problem:")
            print("Start state: ", start)
            print("Goal state: ", goal)
            print("Solution cost encountered: ", cost)
            print("Solution cost expected: ", solution_costs[i])
            print()
        
        cost, expanded_biastar = biastar(start, goal, gridded_map)
        nodes_expanded_biastar.append(expanded_biastar)
        
        if cost != solution_costs[i]:
            print("There is a mismatch in the solution cost found by Bi-A* and what was expected for the problem:")
            print("Start state: ", start)
            print("Goal state: ", goal)
            print("Solution cost encountered: ", cost)
            print("Solution cost expected: ", solution_costs[i])
            print()
        
    print('Finished running all tests. The implementation of an algorithm is likely correct if you do not see mismatch messages for it.')

    if plots:
        from search.plot_results import PlotResults
        plotter = PlotResults()
        plotter.plot_results(nodes_expanded_mm, nodes_expanded_astar, "Nodes Expanded (MM)", "Nodes Expanded (A*)", "nodes_expanded_mm_astar")
        plotter.plot_results(nodes_expanded_mm, nodes_expanded_biastar, "Nodes Expanded (MM)", "Nodes Expanded (Bi-A*)", "nodes_expanded_mm_biastar")

def octileH(state1, state2):
    xdiff = abs(state1.get_x() - state2.get_x())
    ydiff = abs(state1.get_y() - state2.get_y())
    h = 1.5*min(xdiff, ydiff) + abs(xdiff-ydiff)
    return h

def astar(start, goal, gridded_map):
    #print(octileH(start, goal))
    cost = 0
    expanded_astar = 0
    # create open list (heap)
    openList = []
    # create closed list (python dictionary)
    closedList = {}
    start.set_cost(octileH(start, goal))
    heapq.heappush(openList,start) # insert start to open list
    closedList[start.state_hash()] = start # insert start to closed list
    while openList!=[]:
        n_exp = heapq.heappop(openList)
        expanded_astar +=1
        if n_exp==goal:
            cost = n_exp.get_cost()
            return cost, expanded_astar
        children = gridded_map.successors(n_exp)
        
        for n in children:
            hash_value = n.state_hash()
            if hash_value not in closedList:
                n.set_cost(n.get_g()+octileH(n, goal))
                heapq.heappush(openList,n)
                closedList[hash_value] = n
            newf = n.get_g()+ octileH(n, goal)
            oldf = closedList[hash_value].get_cost()
            if hash_value in closedList and newf < oldf:
                n.set_cost(newf)
                heapq.heappush(openList,n)
                closedList[hash_value].set_cost(newf) # need to update cost in closedList? Not in psudocode
                heapq.heapify(openList)
    return -1, expanded_astar

def biastar(start, goal, gridded_map):
    cost = 0
    expanded_biastar = 0
    openf = []
    openb = []
    closedf = {}
    closedb = {}

    start.set_cost(octileH(start, goal))
    goal.set_cost(octileH(goal, start))
    heapq.heappush(openf,start)
    heapq.heappush(openb,goal)
    cost = float('inf')
    while openf!=[] and openb!=[]:
        minf = min(openf).get_cost()
        minb = min(openb).get_cost()
        if cost <= min(minf,minb):
            return cost, expanded_biastar
        if minf < minb:
            # forward
            n_exp = heapq.heappop(openf)
            expanded_biastar+=1
            children = gridded_map.successors(n_exp)
            for n in children:
                hash_value = n.state_hash()
                if hash_value in closedb:
                    cost = min(cost, n.get_g() + closedb[hash_value].get_g())
                if hash_value not in closedf:
                    n.set_cost(n.get_g() + octileH(n, goal))
                    heapq.heappush(openf,n)
                    closedf[hash_value] = n
                newf = n.get_g()+ octileH(n, goal)
                oldf = closedf[hash_value].get_cost()
               
                if hash_value in closedf and newf < oldf:
                    n.set_cost(newf)
                    heapq.heappush(openf,n)
                    closedf[hash_value].set_cost(newf) # need to update cost in closedList? Not in psudocode
                    #heapq.heapify(openf)
        else:
            #backward
            n_exp = heapq.heappop(openb)
            expanded_biastar+=1
            children = gridded_map.successors(n_exp)
            for n in children:
                hash_value = n.state_hash()
                if hash_value in closedf:
                    cost = min(cost, n.get_g() + closedf[hash_value].get_g())
                if hash_value not in closedb:
                    n.set_cost(n.get_g() + octileH(n, start))
                    heapq.heappush(openb,n)
                    closedb[hash_value] = n
                newf = n.get_g()+ octileH(n, start)
                oldf = closedb[hash_value].get_cost()
                if hash_value in closedb and newf < oldf:
                    n.set_cost(newf)
                    heapq.heappush(openb,n)
                    closedb[hash_value].set_cost(newf) # need to update cost in closedList? Not in psudocode
                    #heapq.heapify(openb)
    return  -1, expanded_biastar

def mm(start, goal, gridded_map):
    cost = 0
    expanded_mm = 0
    openf = []
    openb = []
    closedf = {}
    closedb = {}

    start.set_cost(octileH(start, goal))
    goal.set_cost(octileH(goal, start))
    heapq.heappush(openf,start)
    heapq.heappush(openb,goal)
    cost = float('inf')
    while openf!=[] and openb!=[]:
        minf = min(openf).get_cost()
        minb = min(openb).get_cost()
        if cost <= min(minf,minb):
            return cost, expanded_mm
        if minf < minb:
            # forward
            n_exp = heapq.heappop(openf)
            expanded_mm+=1
            children = gridded_map.successors(n_exp)
            for n in children:
                hash_value = n.state_hash()
                if hash_value in closedb:
                    cost = min(cost, n.get_g() + closedb[hash_value].get_g())
                if hash_value not in closedf:
                    p = max(n.get_g() + octileH(n, goal), 2*n.get_g())
                    n.set_cost(p)
                    heapq.heappush(openf,n)
                    closedf[hash_value] = n
                
                newp = max(n.get_g()+ octileH(n, goal), 2*n.get_g())
                oldp = closedf[hash_value].get_cost()
               
                if hash_value in closedf and newp < oldp:
                    n.set_cost(newp)
                    heapq.heappush(openf,n)
                    closedf[hash_value].set_cost(newp) # need to update cost in closedList? Not in psudocode
                    #heapq.heapify(openf)
        else:
            #backward
            n_exp = heapq.heappop(openb)
            expanded_mm+=1
            children = gridded_map.successors(n_exp)
            for n in children:
                hash_value = n.state_hash()
                if hash_value in closedf:
                    cost = min(cost, n.get_g() + closedf[hash_value].get_g())
                if hash_value not in closedb:
                    p = max(n.get_g() + octileH(n, start), 2*n.get_g())
                    n.set_cost(p)
                    heapq.heappush(openb,n)
                    closedb[hash_value] = n
                newp = max(n.get_g()+ octileH(n, start), 2*n.get_g())
                oldp = closedb[hash_value].get_cost()
                if hash_value in closedb and newp < oldp:
                    n.set_cost(newp)
                    heapq.heappush(openb,n)
                    closedb[hash_value].set_cost(newp) # need to update cost in closedList? Not in psudocode
                    #heapq.heapify(openb)
    return  -1, expanded_mm

if __name__ == "__main__":
    main()