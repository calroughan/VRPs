# Get all the PuLP and Dippy functionality
from pulp import *
import coinor.dippy as dippy

# Get some useful math functions
from math import floor, ceil
# Get the drawing functionality
import matplotlib.pyplot as plt
# Get the TSP problem figure size and subtour function
from tsp_prob import FIGSIZE, get_subtour_tsp


# Formulate the TSP as the Assignment Problem,
# i.e., not subtour elimination constraints


def formulate_tsp(tsp, options={}):
    prob = dippy.DipProblem("TSP",
                            # display_mode='matplotlib',
                            display_mode='None',
                            display_interval=10)

    assign_vars = LpVariable.dicts("y",
                                   [(i, j) for i in tsp.EXTLOCS
                                    for j in tsp.EXTLOCS
                                    if i != j],
                                   cat=LpBinary)

    prob += lpSum(tsp.dist[i, j] * assign_vars[i, j]
                  for i in tsp.EXTLOCS
                  for j in tsp.EXTLOCS
                  if i != j), "min_dist"

    for j in tsp.EXTLOCS:
        # One arc in
        prob += lpSum(assign_vars[i, j] for i in tsp.EXTLOCS
                      if i != j) == 1
    for i in tsp.EXTLOCS:
        # One arc out
        prob += lpSum(assign_vars[i, j] for j in tsp.EXTLOCS
                      if j != i) == 1

        # Attach the problem data and variable dictionaries to the DipProblem
    prob.tsp = tsp
    prob.assign_vars = assign_vars

    if "Tol" in options:
        prob.tol = options["Tol"]
    else:
        prob.tol = pow(pow(2, -24), 2.0 / 3.0)

    return prob


# Solve the TSP
def solve_tsp(prob, options={}):
    # Set the options
    prob.options = options
    # When branching, use a callback (just for
    # showing a fractional TSP solution)
    prob.branch_method = user_branch_tsp                                                                             # #
    # When checking feasibility, use a callback
    # to check for subtours
    prob.is_solution_feasible = is_solution_feasible_tsp                                                             # #
    # When generating cuts, use a callback
    # to generate subtour elimination constraints
    prob.generate_cuts = generate_cuts_tsp                                                                           # #

    dippyOpts = {
        #               'CutCGL': 1, # <----- Cuts turned on
        'CutCGL': 0,  # <----- Cuts turned off
        #               'LogDumpModel': 5,
        #               'LogDebugLevel': 5,
    }
    # Can use Cut Generator Library (CGL) cuts too
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


def solve_and_display_tsp(prob, options={}):
    xopt = solve_tsp(prob, options)                                                                                  # #

    # Reads and displays the solution if one is found
    if xopt is not None:
        for var in prob.variables():
            if abs(xopt[var]) > options["Tol"]:
                print(var.name, "=", xopt[var])
    else:
        print("Dippy could not find and optimal solution")

        # Draw the final B-&-B tree if it is being displayed
    if prob.display_mode != 'off':
        tree_nodes = prob.Tree.get_node_list()
        numNodes = len(tree_nodes)
        print("Number of nodes =", numNodes)
        print("Tree depth =", max(prob.Tree.get_node(n).attr['level'] for n in tree_nodes) - 1)
        if prob.display_mode in ['pygame', 'xdot', 'matplotlib']:
            prob.Tree.display(pause=True, wait_for_click=False)

    return xopt


def user_branch_tsp(prob, sol):
    # User callback for branching doesn't find a bracnhc but...
    # Display the current node solution
    assignments = get_assignments_tsp(prob, sol, prob.tol)                                                           # #
    prob.tsp.setSolution(assignments, prob.tol)  # All non-zero arcs
    # prob.tsp.displaySolution(title="Branching", showProb=False)

    return None


def generate_cuts_tsp(prob, sol):
    # User callback for generating cuts

    # Extract and print the current solution value
    obj_val = 0.0
    for v in prob.objective:
        obj_val += prob.objective[v] * sol[v]
    print("In generate_cuts...")
    print("Solution value =", obj_val)
    sys.stdout.flush()

    # Extract and print the current solution
    for var in prob.variables():
        if abs(sol[var]) > prob.tol:
            print(var.name, "=", sol[var])

    # No constraints added
    cons = []
    cons_added = 0

    # Get the assignment variables
    assign_vars = prob.assign_vars
    # Get the threshold for whether an arc should be considered
    # as part of the solution (almost = 1 by default)
    if "Tours" in prob.options:
        threshold = prob.options["Tours"]
    else:
        threshold = 1.0 - prob.tol  # Default is only consider integer arcs

    # Get the solution nodes and arcs
    nodes = prob.tsp.EXTLOCS[:]
    arcs = [(i, j) for (i, j) in assign_vars.keys() if sol[assign_vars[i, j]] > threshold]

    # Define the set of nodes that have not been put in a connected component
    not_connected = set(nodes[:])
    # While any nodes are not in a connected component
    while not_connected:
        # Get an unconnected node
        start = not_connected.pop()
        # Find a subtour from that node
        tNodes, tArcs = get_subtour_tsp(nodes, arcs, start)                                                          # #
        # If it is a subtour (and not a complete tour)
        if len(tNodes) == len(tArcs) and \
                len(tNodes) < len(nodes):
            # Add a subtour elimination constraint
            cons_added += 1
            #        Option 1
            #        cons.append( lpSum(assign_vars[i, j]
            #                           for (i, j) in tArcs) \
            #                     <= len(tArcs) - 1 )
            #       Option 2
            cons.append(lpSum(assign_vars[i, j]
                              for i in tNodes
                              for j in set(nodes).difference(tNodes)) +
                        lpSum(assign_vars[j, i]
                              for i in tNodes
                              for j in set(nodes).difference(tNodes))
                        >= 2)
            print("Subtour elimination!", cons[-1])
            # Return one subtour elimination constraint at a time
            # (can comment out as multiple cuts can be added at once)
            if cons_added == 1: return cons
        # Remove the subtour nodes as they are now connected
        not_connected -= set(tNodes)

    if len(cons) > 0:
        return cons
    else:
        return None


def is_solution_feasible_tsp(prob, sol, tol):
    # User callback for checking feasibility

    # Extract and print the current solution value
    obj_val = 0.0
    for v in prob.objective:
        obj_val += prob.objective[v] * sol[v]
    print("In is_solution_feasible...")
    print("Solution value =", obj_val)
    sys.stdout.flush()

    # Extract and print the current solution
    for var in prob.variables():
        if abs(sol[var]) > tol:
            print(var.name, "=", sol[var])

    # Get the assignment variables
    assign_vars = prob.assign_vars

    # Display the current node solution
    assignments = get_assignments_tsp(prob, sol, tol)                                                                # #
    prob.tsp.setSolution(assignments, tol)
    # prob.tsp.displaySolution(title="Feasibility Check", showProb=None)

    # Get the threshold for whether an arc should be considered
    # as part of the solution (almost = 1 by default)
    if "Tours" in prob.options:
        threshold = prob.options["Tours"]
    else:
        threshold = 1.0 - tol  # Default is only consider integer arcs

    # Get the solution nodes and arcs
    nodes = prob.tsp.EXTLOCS[:]
    arcs = [(i, j) for (i, j) in assign_vars.keys() if sol[assign_vars[i, j]] > threshold]
    # Find a subtour from the first node
    tNodes, tArcs = get_subtour_tsp(nodes, arcs, nodes[0])                                                           # #
    #   print(tNodes, tArcs)
    # If it is a subtour then the problem is infeasible, so will generate cuts
    if (len(tNodes) == len(tArcs)) and (len(tNodes) < len(nodes)):
        print("Solution has subtours!")
        return False

    # No subtours!
    print("Solution has no subtours!")
    return True


# Get the arcs assigned to the tour by
# looking at the TSP formulation and the solution
def get_assignments_tsp(prob, sol, tol):
    assignments = {}
    # Get the arcs (= tup) and the assignment variables
    # for each arc
    for tup, var in prob.assign_vars.items():
        # Is there a current value for that arc?
        if sol[var] is not None:
            # Is it non-zero?
            if sol[var] > tol:
                # Add that assignment
                assignments[tup] = sol[var]

    # Return the current value of the arc assignments (might be fractional)
    return assignments
