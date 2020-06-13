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
    # "Node": "Rounding",
    # "Cuts": "CGL",
    # # # "Cuts": "Lazy Tight"
    }


# Formulate the IP and necessary constraints
def formulate(vrp, options={}):
    prob = dippy.DipProblem("VRP",
                            # display_mode='matplotlib',
                            display_mode='none',
                            display_interval=10)

    assign_vars = LpVariable.dicts("y",
                                   [(i, j, k) for i in vrp.EXTLOCS
                                    for j in vrp.EXTLOCS
                                    for k in vrp.VEHS
                                    if i != j],
                                   cat=LpBinary)
    use_vars = LpVariable.dicts("x", vrp.VEHS, cat=LpBinary)

    # Objective function: minimise the distance between nodes * whether that arc is used by any vehicle.
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
        # Conservation of flows
        # If an arc enters a certain node j from any other node, then there must be
        # an arc leaving j to any other node.
        for j in vrp.LOCS:
            prob += lpSum(assign_vars[i_1, j, k]
                          for i_1 in vrp.EXTLOCS
                          if i_1 != j) == lpSum(assign_vars[j, i_2, k]
                                                for i_2 in vrp.EXTLOCS
                                                if i_2 != j)

        # If all ncurr vehicles specified in the veh_rout_cart[i].py are to be used
        if vrp.allused:

            # Specify that all vehicles must enter the depot
            prob += lpSum(assign_vars[i, 'O', k]
                          for i in vrp.LOCS) == 1

            # Specify all vehicles must leave the depot
            prob += lpSum(assign_vars['O', j, k]
                          for j in vrp.LOCS) == 1

        else:

            # # Moved this into the vrp.distcap set of constraints
            # # Specify that if a vehicle is used it must enter the depot
            # prob += lpSum(assign_vars[i, 'O', k]
            #               for i in vrp.LOCS) == use_vars[k]

            # Specify that if a vehicle is used it must leave the depot
            prob += lpSum(assign_vars['O', j, k]
                          for j in vrp.LOCS) == use_vars[k]

        # Condition for checking if the route taken by each vehicle does not exceed the allowed maximum
        # journey distance
        if vrp.distcap is not None:

            # For each vehicle k, ensure that the maximum distance travelled is less than the distance
            # capacity and 0 if that vehicle is not used.
            prob += lpSum(vrp.dist[i, j] * assign_vars[i, j, k]
                          for i in vrp.EXTLOCS
                          for j in vrp.EXTLOCS
                          if i != j) <= vrp.distcap * use_vars[k]

        else:

            # Strangely returns better solutions with this isolated here.
            # Specify that if a vehicle is used it must enter the depot
            if not vrp.allused:
                prob += lpSum(assign_vars[i, 'O', k]
                              for i in vrp.LOCS) == use_vars[k]

            # Cardinality of arcs for vehicles in use
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


# Solve the TSP
def solve(prob, options={}):

    # Set the options
    prob.options = options

    # When checking feasibility, use a callback
    # to check for subtours
    prob.is_solution_feasible = is_solution_feasible
    # When generating cuts, use a callback
    # to generate subtour elimination constraints
    prob.generate_cuts = generate_cuts
    prob.branch_method = my_branch                                            ###############

    if ("Root" in options) and (options["Root"] is not None):                 ###############
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

    dippyOpts = {}
    #               'CutCGL': 1, # <----- Cuts turned on
    # 'CutCGL': 0,  # <----- Cuts turned off
    #               'LogDumpModel': 5,
    #               'LogDebugLevel': 5,
    # }
    # Can use Cut Generator Library (CGL) cuts too
    if "Cuts" in options:
        if options["Cuts"] == "CGL":
            dippyOpts['CutCGL'] = 1
        elif options["Cuts"] == "Lazy Tight":
            dippyOpts['CutLazyTight'] = 1

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

    # Reads and displays the solution if one is found
    if xopt is not None:
        for var in prob.variables():
            if abs(xopt[var]) > options["Tol"]:
                print(var.name, "=", xopt[var])
    else:
        print("Dippy could not find and optimal solution")

    # Draw the final B-&-B tree if it is being displayed
    # if prob.display_mode != 'off':
    #     tree_nodes = prob.Tree.get_node_list()
    #     numNodes = len(tree_nodes)
    #     print("Number of nodes =", numNodes)
    #     print("Tree depth =", max(prob.Tree.get_node(n).attr['level'] for n in tree_nodes) - 1)
    #     if prob.display_mode in ['pygame', 'xdot', 'matplotlib']:
    #         prob.Tree.display(pause=True, wait_for_click=False)

    return xopt


# User callback for generating cuts
def generate_cuts(prob, sol):

    # No constraints added
    cons = []
    cons_added = 0

    # Get the assignment variables and values
    assign_vars = prob.assign_vars
    assign_vals = dict([((i, j, k), sol[assign_vars[i, j, k]]) for (i, j, k) in assign_vars.keys()])

    # Get the threshold for whether an arc should be considered
    # as part of the solution (almost = 1 by default)
    if "Tours" in prob.options:
        threshold = prob.options["Tours"]
    else:
        threshold = 1.0 - prob.tol  # Default is only consider integer arcs

    # Get the graphs for each vehicle
    nodes = prob.vrp.EXTLOCS[:]
    arcs = [(i, j, k) for (i, j, k) in assign_vars.keys() if sol[assign_vars[i, j, k]] > threshold]

    # Loop over the vehicles that are used
    for k in prob.vrp.VEHS:

        # Extract the arcs for each vehicle
        vehArcs = [x[:2] for x in arcs if x[2] == k]

        # Define the set of nodes that have not been put in a connected component
        not_connected = set(nodes[:])

        # While any nodes are not in a connected component
        while not_connected:

            # Get an unconnected node
            start = not_connected.pop()

            # Find a subtour from that starting node
            tNodes, tArcs = get_subtour(nodes, vehArcs, start)

            # If it is a subtour (and not a complete tour), add a subtour elimination constraint provided that
            # the depot is not included in the subtour,
            if len(tNodes) == len(tArcs) and len(tNodes) < len(nodes) and ('O' not in tNodes):
                cons_added += 1

                # If a subtour is found then that graph must be banned

                # Option 1
                cons.append(lpSum(assign_vars[i, j, k]
                                  for (i, j) in tArcs) <= len(tArcs) - 1)

                # Option 2
                # cons.append(lpSum(assign_vars[i, j, k]
                #                   for i in tNodes
                #                   for j in set(nodes).difference(tNodes)) +
                #             lpSum(assign_vars[j, i, k]
                #                   for i in tNodes
                #                   for j in set(nodes).difference(tNodes))
                #             >= 2)

                # print("Subtour elimination!", cons[-1])

                # Return one subtour elimination constraint at a time
                if cons_added == 1:
                    return cons

            # Remove the subtour nodes as they are now connected
            not_connected -= set(tNodes)

    if len(cons) > 0:
        return cons
    else:
        return None


# User callback for checking feasibility
def is_solution_feasible(prob, sol, tol):

    # Display feasibility checks if desired
    # assignments = get_assignments(prob, sol, tol)
    # prob.vrp.setSolution(assignments, tol)
    # prob.vrp.displaySolution(title="Feasibility Check", showProb=None)

    # Get the threshold for whether an arc should be considered
    # as part of the solution (almost = 1 by default)
    if "Tours" in prob.options:
        threshold = prob.options["Tours"]
    else:
        threshold = 1.0 - prob.tol  # Default is only consider integer arcs

    # Get the assignment variables
    assign_vars = prob.assign_vars
    assign_vals = dict([((i, j, k), sol[assign_vars[i, j, k]])
                       for (i, j, k) in assign_vars.keys()])

    # Get the graphs for each vehicle
    nodes = prob.vrp.EXTLOCS[:]
    arcs = [(i, j, k) for (i, j, k) in assign_vars.keys() if sol[assign_vars[i, j, k]] > threshold]

    # Loop over the vehicles that are used
    for k in prob.vrp.VEHS:

        # Extract the arcs for each vehicle to pass into get_subtour()
        vehArcs = [x[:2] for x in arcs if x[2] == k]

        # Define the set of nodes that have not been put in a connected component
        not_connected = set(nodes[:])

        # While any nodes are not in a connected component
        while not_connected:

            # Get an unconnected node
            start = not_connected.pop()

            # Find a subtour from that node
            tNodes, tArcs = get_subtour(nodes, vehArcs, start)
            # tNodes, tArcs = get_subtour(nodes, vehArcs, nodes[0])

            #   If a subtour is found then the solution is not feasible, so will declare it as such
            if (len(tNodes) == len(tArcs)) and (len(tNodes) < len(nodes)) and ('O' not in tNodes):
                # print("Solution has subtours!")
                return False

            # Remove the subtour nodes as they are now connected
            not_connected -= set(tNodes)

    # Otherwise it is feasible
    # print("Solution has no subtours!")
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

    if (bounds is None) and ("Priority" in options) and (options["Priority"] == "x"):
        # If symmetry found nothing / is off, run this
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

    else:
        if prob.node_heuristic:  # Are we using a node heuristic?
            if ("Node" in options) and (options["Node"] == "Rounding"):
                sol = improvement_tsp(prob, xhat)

    if sol is not None:
        return [sol]


def root_tsp(prob):

    # The maxdist root node heuristic is poor, don't use it
    if prob.vrp.distcap is not None:
        return None

    vrp         = prob.vrp
    use_vars    = prob.use_vars
    assign_vars = prob.assign_vars
    LOCS        = vrp.LOCS
    allused     = vrp.allused
    maxdist     = vrp.distcap
    x           = vrp.x
    y           = vrp.y

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


def improvement_tsp(prob, xhat):

    if prob.vrp.distcap is not None or prob.vrp.allused:
        return None

    vrp = prob.vrp
    use_vars = prob.use_vars
    assign_vars = prob.assign_vars

    sol = {}
    # Initialise sol
    for k in vrp.VEHS:
        sol[use_vars[k]] = 0
        for i in vrp.EXTLOCS:
            for j in vrp.EXTLOCS:
                if i != j:
                    sol[assign_vars[i, j, k]] = 0

    if "Tours" in prob.options:
        threshold = prob.options["Tours"]
    else:
        threshold = 1.0 - prob.tol  # Default is to only consider integer arcs

    # Split xhat into the vehicles and arcs
    vehicles = dict({(key, val) for (key, val) in xhat.items() if str(key)[0] == "x"})
    arcs     = dict({(key, val) for (key, val) in xhat.items() if str(key)[0] == "y"})

    # Force vehicles to be on / off
    for key, val in vehicles.items():
        if val >= 0.3:
            sol[key] = 1
        else:
            sol[key] = 0

    # If no arcs are on, turn on the first by default
    if len({val for (key, val) in sol.items() if val == 1}) == 0:
        sol[use_vars[1]] = 1

    # Force arcs to be on / off
    for key, val in arcs.items():

        if val >= 0.5:
            sol[key] = 1
        else:
            sol[key] = 0

    zz = dict({(key, val) for (key, val) in sol.items() if val > 0.95})

    if is_solution_feasible(prob, sol, prob.tol):
        return sol
    else:
        return None


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
