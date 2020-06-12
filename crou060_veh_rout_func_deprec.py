from pulp import *
import coinor.dippy as dippy
import random

from math import floor, ceil
import matplotlib.pyplot as plt
from veh_rout_prob import FIGSIZE, get_graphs, get_subtour

from tsp_prob import TSPProb, get_subtour_tsp
from tsp_func import formulate_tsp, solve_and_display_tsp, get_assignments_tsp

tol = pow(pow(2, -20), 2.0 / 3.0)
myopts = {
    "Tol": tol,
    "Interval": 1000,
    "Antisymmetry": "on",
    "Tight": "on",
    "Aggregate": "on",
    "Priority": "x",
    "Root": "Root TSP",
    # "Node": "Frac Fit",
    "Cuts": "CGL"
}


def formulate(vrp, options={}):
    prob = dippy.DipProblem("VRP",
                            # display_mode='matplotlib',
                            display_mode='None',
                            display_interval=10)

    assign_vars = LpVariable.dicts("y",
                                   [(i, j, k) for i in vrp.EXTLOCS
                                    for j in vrp.EXTLOCS
                                    for k in vrp.VEHS
                                    if i != j],
                                   cat=LpBinary)
    use_vars = LpVariable.dicts("x", vrp.VEHS, cat=LpBinary)

    # Objective function
    prob += lpSum(vrp.dist[i, j] * assign_vars[i, j, k]
                  for i in vrp.EXTLOCS
                  for j in vrp.EXTLOCS
                  for k in vrp.VEHS
                  if i != j), "min_dist"

    # Each node (excluding 'O') must have one arc entering from any other node (including 'O')
    for j in vrp.LOCS:
        prob += lpSum(assign_vars[i, j, k]
                      for i in vrp.EXTLOCS
                      for k in vrp.VEHS
                      if i != j) == 1

    # Each node (excluding 'O') must have one arc leaving to any other node (including 'O')
    for i in vrp.LOCS:
        prob += lpSum(assign_vars[i, j, k]
                      for j in vrp.EXTLOCS
                      for k in vrp.VEHS
                      if j != i) == 1

    for k in vrp.VEHS:

        # Specify all vehicles must leave the depot
        if vrp.allused:
            prob += lpSum(assign_vars['O', j, k] for j in vrp.LOCS) >= 1

        # Conservation of flows
        for j in vrp.LOCS:
            prob += lpSum(assign_vars[i_1, j, k] for i_1 in vrp.EXTLOCS
                          if i_1 != j) == lpSum(assign_vars[j, i_2, k]
                                                for i_2 in vrp.EXTLOCS
                                                if i_2 != j)

    for k in vrp.VEHS:
        prob += lpSum(assign_vars[i, 'O', k] for i in vrp.LOCS) == use_vars[k]
        prob += lpSum(assign_vars['O', j, k] for j in vrp.LOCS) == use_vars[k]

        # For each vehicle k, ensure that the maximum distance travelled is less than the distance
        # capacity and 0 if that vehicle is not used.
        if vrp.distcap is not None:
            prob += lpSum(vrp.dist[i, j] * assign_vars[i, j, k]
                          for i in vrp.EXTLOCS
                          for j in vrp.EXTLOCS
                          if i != j) <= vrp.distcap * use_vars[k]

        # Cardinality of arcs for vehicles in use
        else:
            prob += lpSum(assign_vars[i, j, k]
                          for i in vrp.EXTLOCS
                          for j in vrp.EXTLOCS
                          if i != j) <= len(vrp.EXTLOCS) * use_vars[k]

    if ('Antisymmetry' in options) and (options['Antisymmetry'] == 'on'):
        # Order the use of the vehicles
        for k in vrp.VEHS:
            if k != vrp.VEHS[-1]:
                prob += use_vars[k] >= use_vars[k + 1]

    if ('Tight' in options) and (options['Tight'] == 'on'):
        for k in vrp.VEHS:
            for i in vrp.EXTLOCS:
                for j in vrp.EXTLOCS:
                    if i != j:
                        prob += assign_vars[i, j, k] <= use_vars[k]

    # Attach the problem data and variable dictionaries to the DipProblem
    prob.vrp = vrp
    prob.assign_vars = assign_vars
    prob.use_vars = use_vars

    if "Tol" in options:
        prob.tol = options["Tol"]
    else:
        prob.tol = pow(pow(2, -24), 2.0 / 3.0)

    return prob


def solve(prob, options={}):
    prob.options = options
    prob.is_solution_feasible = is_solution_feasible
    prob.generate_cuts = generate_cuts
    prob.branch_method = my_branch

    if ("Root" in options) and (options["Root"] is not None):
        prob.is_root_node = True
        prob.root_heuristic = True
        prob.heuristics = my_heuristics
    else:
        prob.root_heuristic = False

    if ("Node" in options) and (options["Node"] is not None):
        prob.is_root_node = True
        prob.node_heuristic = True
        prob.heuristics = my_heuristics
    else:
        prob.node_heuristic = False

    dippyOpts = {
        # 'CutCGL': 1,  # <----- Cuts turned on
        'CutCGL': 0,  # <----- Cuts turned off
        # 'LogDumpModel': 5,
        # 'LogDebugLevel': 5,
    }
    if ("Cuts" in options) and (options["Cuts"] == "CGL"):
        dippyOpts['CutCGL'] = 1
    if "Interval" in options:
        prob.display_interval = options["Interval"]

    # plt.figure(figsize=FIGSIZE)
    status, message, primals, duals = dippy.Solve(prob, dippyOpts)

    if status == LpStatusOptimal:
        return dict((var, var.value()) for var in prob.variables())
    else:
        return None


def solve_and_display(prob, options={}):
    xopt = solve(prob, options)

    if xopt is not None:
        for var in prob.variables():
            if abs(xopt[var]) > options["Tol"]:
                print(var.name, "=", xopt[var])
    else:
        print("Dippy could not find and optimal solution")

    if prob.display_mode != 'off':
        tree_nodes = prob.Tree.get_node_list()
        numNodes = len(tree_nodes)
        print("Number of nodes =", numNodes)
        print("Tree depth =", max(prob.Tree.get_node(n).attr['level'] for n in tree_nodes) - 1)
        if prob.display_mode in ['pygame', 'xdot', 'matplotlib']:
            prob.Tree.display(pause=True, wait_for_click=False)

    return xopt


def generate_cuts(prob, sol):

    # No constraints added initially
    cons = []
    cons_added = 0

    # Get the assignment variables
    assign_vars = prob.assign_vars
    assign_vals = dict([((i, j, k), sol[assign_vars[i, j, k]]) for (i, j, k) in assign_vars.keys()])

    # Get the threshold for whether an arc should be considered
    # as part of the solution (almost = 1 by default)
    if "Tours" in prob.options:
        threshold = prob.options["Tours"]
    else:
        threshold = 1.0 - prob.tol  # Default is only consider integer arcs

    # Get the graph structure for each vehicle
    vehNodes, vehArcs = get_graphs(prob.vrp, assign_vals, threshold)

    # Loop over the vehicles that are used
    for k in prob.vrp.VEHS:

        # Extract the nodes and arcs for each vehicle
        if (sol[prob.use_vars[k]] > prob.tol):
            kNodes = vehNodes[k]
            kArcs = vehArcs[k]

            # Define the set of nodes that have not been put in a connected component
            not_connected = set(kNodes)

            # While any nodes are not in a connected component
            while not_connected:

                # Get an unconnected node
                start = not_connected.pop()

                # Find a subtour from that node
                nodes, arcs = get_subtour(kNodes, kArcs, start)

                # If it is a subtour (and not a complete tour), add a subtour elimination constraint
                if (len(nodes) == len(arcs)) and (len(nodes) < len(kNodes)):
                    cons_added += 1

                    # Option 1
                    cons.append(lpSum(assign_vars[i, j, k]
                                      for (i, j) in kArcs) <= len(kArcs) - 1)

                # Remove the subtour nodes as they are now connected
                not_connected -= set(nodes)

    if len(cons) > 0:
        return cons
    else:
        return None


def is_solution_feasible(prob, sol, tol):


    assign_vars = prob.assign_vars
    assign_vals = dict([((i, j, k), sol[assign_vars[i, j, k]])
                        for (i, j, k) in assign_vars.keys()])

    if "Tours" in prob.options:
        threshold = prob.options["Tours"]
    else:
        threshold = 1.0 - prob.tol  # Default is only consider integer arcs

    # Get the graphs for each vehicle
    vehNodes, vehArcs = get_graphs(prob.vrp, assign_vals, threshold)

    # Loop over the vehicles that are used
    for k in prob.vrp.VEHS:

        # Extract the nodes and arcs for each vehicle to pass into get_subtour()
        if (sol[prob.use_vars[k]] > prob.tol):
            kNodes = vehNodes[k]
            kArcs = vehArcs[k]

            if len(kNodes) > 0:

                # Look for subtour in each vehicles graph from the first nod
                nodes, arcs = get_subtour(kNodes, kArcs, kNodes[0])

                # If a subtour is found then the solution is not feasible, so will generate cuts
                if (len(nodes) == len(arcs)) and (len(nodes) < len(kNodes)):
                    return False

    # Else no subtours
    return True


def get_assignments(prob, sol, tol):
    assignments = {}
    for tup, var in prob.assign_vars.items():
        if sol[var] is not None:
            if sol[var] > tol:
                assignments[tup] = sol[var]

    return assignments


def my_branch(prob, sol):

    options = prob.options
    bounds = None

    if ("Aggregate" in options) and (options["Aggregate"] == "on"):
        bounds = symmetry(prob, sol)
    elif ("Priority" in options) and (options["Priority"] == "x"):
        bounds = most_frac_use(prob, sol)

    return bounds


def most_frac_use(prob, sol):

    vrp         = prob.vrp
    use_vars    = prob.use_vars
    tol         = prob.tol
    most        = float('-inf')

    Veh = None
    for k in vrp.VEHS:
        alpha = sol[use_vars[k]]
        up = ceil(alpha)  # Round up to next nearest integer
        down = floor(alpha)  # Round down
        frac = min(up - alpha, alpha - down)
        if frac > tol:  # Is fractional?
            if frac > most:
                most = frac
                Veh = k

    down_lbs = {}
    down_ubs = {}
    up_lbs = {}
    up_ubs = {}

    if Veh is not None:
        down_ubs[use_vars[Veh]] = 0.0
        up_lbs[use_vars[Veh]] = 1.0

        return down_lbs, down_ubs, up_lbs, up_ubs


def symmetry(prob, sol):

    vrp         = prob.vrp
    use_vars    = prob.use_vars
    tol         = prob.tol

    alpha = sum(sol[use_vars[k]] for k in vrp.VEHS)
    up   = int(ceil(alpha))  # Round up to next nearest integer
    down = int(floor(alpha))  # Round down
    frac = min(up - alpha, alpha - down)
    if frac > tol:  # Is fractional?

        down_lbs = {}
        down_ubs = {}
        up_lbs = {}
        up_ubs = {}
        for n in range(up - 1, len(vrp.VEHS)):
            down_ubs[use_vars[vrp.VEHS[n]]] = 0.0
        for n in range(up):
            up_lbs[use_vars[vrp.VEHS[n]]] = 1.0

        return down_lbs, down_ubs, up_lbs, up_ubs


def my_heuristics(prob, xhat, cost):

    options = prob.options
    sol = None

    if prob.is_root_node:
        prob.is_root_node = False
        if prob.root_heuristic:
            if ("Root" in options) and (options["Root"] == "Root TSP"):
                sol = root_tsp(prob)
    #
    # else:
    #     if prob.node_heuristic:  # Are we using a node heuristic?
    #         if ("Node" in options) and (options["Node"] == "Frac Fit"):
    #             sol = frac_fit(prob, xhat)

    if sol is not None:
        return [sol]


def initial_route(prob):

    # One vehicle, unless maxdist or all used specified
    # One vehicle: From depot to the closest node, then add the next closest node, repeat until you have a closed loop
    # Maxdist: Same as One Vehicle, however we move onto the next vehicle when the route is too long
    # Allused: Loop through each vehicle first and then insert a node into the remaining ones

    vrp = prob.vrp
    use_vars = prob.use_vars
    assign_vars = prob.assign_vars
    dist = vrp.dist

    savings = {}

    # Calculate the savings dict, only taking the top triangle
    for i in vrp.LOCS:
        for j in vrp.LOCS:
            if i < j:
                # savingDict[i, j] = dist['O', i] + dist[j, 'O'] - dist[i, j]
                savings[(i, j)] = dist['O', i] + dist[j, 'O'] - dist[i, j]

    # Order the savings dict
    savings = {key: val for key, val in sorted(savings.items(),
                                               key=lambda item: item[1],
                                               reverse=True)}

    print(savings)

    sol = {}

    # Must use all nodes
    for i in vrp.LOCS:
        sol[use_vars[i]] = 1

    # Set all arcs to be empty
    for i in vrp.EXTLOCS:
        for j in vrp.EXTLOCS:
            for k in vrp.VEHS:
                if i != j:
                    sol[assign_vars[i, j, k]] = 0

    # Arbitrarily assign each node to its own route
    k = 1
    for i in vrp.LOCS:
        sol[assign_vars['O', i, k]] = 1
        sol[assign_vars[i, 'O', k]] = 1
        k += 1

    # Iterate through the savings dict, and assign the largest value to a route,
    currentroute = 1          # <- Starting vehicle

    for key, val in savings.items():
        i = key[0]
        j = key[1]

        # If i and j are both adjacent to the depot:
        # if (sum(sol[assign_vars[i, 'O', k]] for k in vrp.VEHS) or sum(sol[assign_vars['O', i, k]] for k in vrp.VEHS)) \
        #         and (sum(sol[assign_vars[j, 'O', k]] for k in vrp.VEHS) or sum(sol[assign_vars['O', j, k]] for k in vrp.VEHS))




        # Are (i, j) in a route together already? If no, assign them to a new route
        # if sum(sol[assign_vars[i, j, k]] for k in vrp.VEHS) == 0:
        #     sol[assign_vars[i, j, currentroute]] = 1

        # If one of (i, j) are in a route already, and that point is adjacent to the depot
        # Else if i is in a route and is adjacent to the depot or j is in a route and adjacent to the depot,
        # add the new arc into the route
        # elif sum(sol[assign_vars[i, j_, k]] for j_ in vrp.EXTLOCS for k in vrp.VEHS)


        # Add the new arc into the route



    # Parallel Method: Allused

    # Start from the top of the savings (maximum saving) list and execute the following.
    # If next edge has a common node with the existing route and the common node is not an interior of the route then
    # connect that edge to the existing route, otherwise start a new route with next edge.
    # Repeat the above step until all nodes serviced or no edges left in sorted edges

    # Sequential Method: Maxdist

    # This is same as parallel method except for one change.In sequential method routes built sequentially.
    # That is, we can not start new routes until existing routes are filled.

    # Add the closest node to the




    list_in = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    n = 3

    random.shuffle(list_in)
    print([list_in[i::n] for i in range(n)])
    return [list_in[i::n] for i in range(n)]


def run_tsp(route, x, y):

    tsp_opts = {
        "Tol": tol,
        "Interval": 1000,
        #      "Tours": tol,
        #      "Tours": 0.5 - tol,
        "Tours": 1.0 - tol,
    }

    # Initialise TSP
    tsp = TSPProb(LOCS=route,
                  x=x,
                  y=y)

    # Formulate_TSP
    prob_tsp = formulate_tsp(tsp, options=tsp_opts)
    prob_tsp.display_mode = 'off'

    # Solve the TSP and display the B-&-B tree
    xopt = solve_and_display_tsp(prob_tsp, options=tsp_opts)

    # Extract the arcs assigned to the tour
    return get_assignments_tsp(prob_tsp, xopt, tol)


def root_tsp(prob):

    vrp         = prob.vrp
    use_vars    = prob.use_vars
    assign_vars = prob.assign_vars
    dist        = vrp.dist
    LOCS        = vrp.LOCS
    allused     = vrp.allused
    maxdist     = vrp.distcap
    x           = vrp.x
    y           = vrp.y
    tol         = prob.tol

    # Initialise sol
    sol = {}

    if allused:
        n = vrp.fixed
        for k in vrp.VEHS:
            sol[use_vars[k]] = 1                    # Must use all vehicles
    else:
        n = 1                                       # Just use one vehicle

    # Partition nodes into n sets
    random.shuffle(LOCS)
    partitions = [LOCS[i::n] for i in range(n)]

    # Set all arcs to be empty
    for i in vrp.EXTLOCS:
        for j in vrp.EXTLOCS:
            for k in vrp.VEHS:
                if i != j:
                    sol[assign_vars[i, j, k]] = 0

    if maxdist is None:                                         # Separate heuristic for maxdist problems
        # Route must include the depot
        for vehicle, route in enumerate(partitions):

            # # Initialise TSP
            # tsp = TSPProb(LOCS=route,
            #               x=x,
            #               y=y)
            #
            # # Formulate_TSP
            # prob_tsp = formulate_tsp(tsp, options=tsp_opts)
            # prob_tsp.display_mode = 'off'
            #
            # # Solve the TSP and display the B-&-B tree
            # xopt = solve_and_display_tsp(prob_tsp, options=tsp_opts)
            #
            # # Extract the arcs assigned to the tour
            # assignments = get_assignments_tsp(prob_tsp, xopt, tol)
            assignments = run_tsp(route, x, y)

            # Update sol with these assignments
            sol[use_vars[vehicle+1]] = 1
            for (i, j) in assignments.keys():
                sol[assign_vars[i, j, vehicle+1]] = 1

    elif maxdist is not None:
        # Start with a partition of two vehicles
        n = 2
        count = 0

        # Start with two partitions, calculate TSP for each, find the total cost of each tour, if any one of the tours
        # breaches the maximum distance, reshuffle the nodes and repeat. If after five reshufflings there still isn't a
        # feasible option, then increase the number of partitions and repeat

        while n < vrp.fixed:

            if count >= 5:
                n += 1
                count = 0

            # Partition LOCS into sets
            random.shuffle(LOCS)
            partitions = [LOCS[i::n] for i in range(n)]

            total = {}
            assignments = {}

            for vehicle, route in enumerate(partitions):

                # # Initialise TSP
                # tsp = TSPProb(LOCS=route,
                #               x=x,
                #               y=y)
                #
                # # Formulate_TSP
                # prob_tsp = formulate_tsp(tsp, options=tsp_opts)
                # prob_tsp.display_mode = 'off'
                #
                # # Solve the TSP and display the B-&-B tree
                # xopt = solve_and_display_tsp(prob_tsp, options=tsp_opts)
                #
                # # Extract the arcs assigned to the tour
                # assignments[vehicle] = get_assignments_tsp(prob_tsp, xopt, tol)
                assignments[vehicle] = run_tsp(route, x, y)

                # Calculate the total_distance of the tour
                total[vehicle] = 0
                for arc in assignments[vehicle]:
                    # print(" ", arc[0], arc[1], vrp.dist[arc[0], arc[1]])
                    total[vehicle] += vrp.dist[arc[0], arc[1]]

            # If all tours have a distance <= maxdist, return this solution
            if len({val for (key, val) in total.items() if val > maxdist}) == 0:

                # Add these to the solution
                for vehicle, route in enumerate(partitions):

                    print(assignments[vehicle].keys())

                    # Update sol with these assignments
                    sol[use_vars[vehicle + 1]] = 1
                    for (i, j) in assignments[vehicle].keys():
                        sol[assign_vars[i, j, vehicle + 1]] = 1

                return sol

    return sol


















