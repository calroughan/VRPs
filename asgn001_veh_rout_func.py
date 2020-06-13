from PIL import Image

Image.MAX_IMAGE_PIXELS = None

from pulp import *
import coinor.dippy as dippy

from math import floor, ceil
import matplotlib.pyplot as plt
from veh_rout_prob import FIGSIZE, get_graphs, get_subtour

tol = pow(pow(2, -20), 2.0 / 3.0)
myopts = {
    "Tol": tol,
    # "Interval": 10,
    "Interval": 1000,
    # "Interval": None,
    "Antisymmetry": "on",
    "Aggregate": "on",
    "Tight": "on",
    "Priority": "x",
    # "Cuts": "CGL"
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

    prob += lpSum(vrp.dist[i, j] * assign_vars[i, j, k]
                  for i in vrp.EXTLOCS
                  for j in vrp.EXTLOCS
                  for k in vrp.VEHS if i != j), "min_dist"

    for j in vrp.LOCS:
        # One arc in
        prob += lpSum(assign_vars[i, j, k] for i in vrp.EXTLOCS
                      for k in vrp.VEHS
                      if i != j) == 1
    for i in vrp.LOCS:
        # One arc out
        prob += lpSum(assign_vars[i, j, k] for j in vrp.EXTLOCS
                      for k in vrp.VEHS
                      if j != i) == 1

    for k in vrp.VEHS:
        if vrp.allused:
            # If all vehicles must be used
            prob += lpSum(assign_vars['O', j, k] for j in vrp.LOCS) >= 1

        for j in vrp.LOCS:
            # Conserve vehicle flow
            prob += lpSum(assign_vars[ifrom, j, k] for ifrom in vrp.EXTLOCS
                          if ifrom != j) == \
                    lpSum(assign_vars[j, ito, k] for ito in vrp.EXTLOCS
                          if ito != j)

    #   # Number of arcs into the origin = number of vehicles being used
    #   prob += lpSum(assign_vars[i, 'O', k] for i in vrp.LOCS
    #                                        for k in vrp.VEHS) \
    #        == lpSum(use_vars[k] for k in vrp.VEHS)
    #   for k in vrp.VEHS:
    #     prob += lpSum(assign_vars[i, 'O', k] for i in vrp.LOCS) \
    #          == use_vars[k]
    #
    #   # Number of arcs out of the origin = number of vehicles being used
    #   prob += lpSum(assign_vars['O', j, k] for j in vrp.LOCS
    #                                      for k in vrp.VEHS) \
    #        == lpSum(use_vars[k] for k in vrp.VEHS)
    for k in vrp.VEHS:
        prob += lpSum(assign_vars[i, 'O', k] for i in vrp.LOCS) \
                == use_vars[k]
        prob += lpSum(assign_vars['O', j, k] for j in vrp.LOCS) \
                == use_vars[k]

    if ('Antisymmetry' in options) and (options['Antisymmetry'] == 'on'):
        # Order the use of the vehicles
        for k in vrp.VEHS:
            if k != vrp.VEHS[-1]:
                prob += use_vars[k] >= use_vars[k + 1]

    if vrp.distcap is not None:
        for k in vrp.VEHS:
            prob += lpSum(vrp.dist[i, j] * assign_vars[i, j, k]
                          for i in vrp.EXTLOCS
                          for j in vrp.EXTLOCS
                          if i != j) <= vrp.distcap * use_vars[k]
    else:
        for k in vrp.VEHS:
            prob += lpSum(assign_vars[i, j, k]
                          for i in vrp.EXTLOCS
                          for j in vrp.EXTLOCS
                          if i != j) <= len(vrp.EXTLOCS) * use_vars[k]

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
    prob.branch_method = my_branch
    prob.is_solution_feasible = is_solution_feasible
    prob.generate_cuts = generate_cuts

    dippyOpts = {
        #               'CutCGL': 1, # <----- Cuts turned on
        'CutCGL': 0,  # <----- Cuts turned off
        #               'LogDumpModel': 5,
        #               'LogDebugLevel': 5,
    }
    if ("Cuts" in options) and (options["Cuts"] == "CGL"):
        dippyOpts['CutCGL'] = 1
    if "Interval" in options:
        prob.display_interval = options["Interval"]

    plt.figure(figsize=FIGSIZE)
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
    '''
    obj_val = 0.0
    for v in prob.objective:
      obj_val += prob.objective[v] * sol[v]
    print("In generate_cuts...")
    print("Solution value =", obj_val)
    sys.stdout.flush()

    for var in prob.variables():
      if abs(sol[var]) > prob.tol:
        print(var.name, "=", sol[var])
    '''

    cons = []
    cons_added = 0

    assign_vars = prob.assign_vars
    assign_vals = dict([((i, j, k), sol[assign_vars[i, j, k]]) for (i, j, k) in assign_vars.keys()])
    if "Tours" in prob.options:
        threshold = prob.options["Tours"]
    else:
        threshold = 1.0 - prob.tol  # Default is only consider integer arcs
    vehNodes, vehArcs = get_graphs(prob.vrp, assign_vals, threshold)
    for k in prob.vrp.VEHS:
        if (sol[prob.use_vars[k]] > prob.tol):
            kNodes = vehNodes[k]
            kArcs = vehArcs[k]
            #        print(k, kNodes, kArcs)
            not_connected = set(kNodes)
            while not_connected:
                start = not_connected.pop()
                nodes, arcs = get_subtour(kNodes, kArcs, start)
                #          print(nodes, arcs)
                #          print(len(nodes), len(arcs), len(kNodes))
                if (len(nodes) == len(arcs)) and (len(nodes) < len(kNodes)):
                    cons_added += 1
                    #            prob.vrp.setSolution(assign_vals, prob.tol)
                    #            prob.vrp.displaySolution(title="Subtour example", showProb=False)
                    cons.append(lpSum(assign_vars[i, j, k]
                                      for (i, j) in kArcs) \
                                <= len(kArcs) - 1)
                #            print("Subtour elimination!", cons[-1])
                #            return cons
                not_connected -= set(nodes)

    if len(cons) > 0:
        return cons
    else:
        return None


def is_solution_feasible(prob, sol, tol):
    '''
    obj_val = 0.0
    for v in prob.objective:
      obj_val += prob.objective[v] * sol[v]
    print("In is_solution_feasible...")
    print("Solution value =", obj_val)
    sys.stdout.flush()
    
    for var in prob.variables():
      if abs(sol[var]) > prob.tol:
        print(var.name, "=", sol[var])
    '''

    assign_vars = prob.assign_vars
    '''
    # Display the current node solution
    prob.vrp.setSolution(assign_vars, sol, prob.tol)
    prob.vrp.displaySolution()
    '''
    assign_vals = dict([((i, j, k), sol[assign_vars[i, j, k]])
                        for (i, j, k) in assign_vars.keys()])
    if "Tours" in prob.options:
        threshold = prob.options["Tours"]
    else:
        threshold = 1.0 - prob.tol  # Default is only consider integer arcs
    vehNodes, vehArcs = get_graphs(prob.vrp, assign_vals, threshold)
    for k in prob.vrp.VEHS:
        if (sol[prob.use_vars[k]] > prob.tol):
            kNodes = vehNodes[k]
            kArcs = vehArcs[k]
            #        print(k, kNodes, kArcs)
            if len(kNodes) > 0:
                nodes, arcs = get_subtour(kNodes, kArcs, kNodes[0])
                #          print(nodes, arcs)
                if (len(nodes) == len(arcs)) and (len(nodes) < len(kNodes)):
                    #            print("Solution has subtours!")
                    return False

    #    print("Solution has no subtours!")
    return True


def get_assignments(prob, sol, tol):
    assignments = {}
    for tup, var in prob.assign_vars.items():
        if sol[var] is not None:
            if sol[var] > tol:
                assignments[tup] = sol[var]

    return assignments


def my_branch(prob, sol):
    #
    #   print "Custom branching, not doing anything yet!"
    #   bounds = None
    #
    #   return bounds
    options = prob.options

    #  print "Custom branching, aggregate branch and prioritising variables!"
    bounds = None

    #   if  ( ("Symmetry" in options) and (options["Symmetry"] == "on") ) \
    #   and ( ("Aggregate" in options) and (options["Aggregate"] == "on") ):
    if ("Aggregate" in options) and (options["Aggregate"] == "on"):
        bounds = symmetry(prob, sol)

    if "Priority" in options:
        if (options["Priority"] == "x") or (options["Priority"].startswith("xy")):
            if bounds is None:
                bounds = most_frac_use(prob, sol)

        if options["Priority"].startswith("y") or options["Priority"].startswith("xy"):
            if bounds is None:
                if "ym" in options['Priority']:
                    bounds = frac_assign(prob, sol, most=True)
                elif "yl" in options['Priority']:
                    bounds = frac_assign(prob, sol, most=False)
            if 'x' in options["Priority"]:
                if bounds is None:
                    bounds = most_frac_use(prob, sol)

    return bounds


def most_frac_use(prob, sol):
    # Get the attached data and variable dicts
    vrp = prob.vrp
    use_vars = prob.use_vars
    tol = prob.tol

    most = float('-inf')
    veh = None
    for k in vrp.VEHS:
        alpha = sol[use_vars[k]]
        up = ceil(alpha)  # Round up to next nearest integer
        down = floor(alpha)  # Round down
        frac = min(up - alpha, alpha - down)
        if frac > tol:  # Is fractional?
            if frac > most:
                most = frac
                veh = k

    down_lbs = {}
    down_ubs = {}
    up_lbs = {}
    up_ubs = {}
    if veh is not None:
        #    print veh, sol[use_vars[veh]]
        down_ubs[use_vars[veh]] = 0.0
        up_lbs[use_vars[veh]] = 1.0

        return down_lbs, down_ubs, up_lbs, up_ubs


def frac_assign(prob, sol, most=True):
    # Get the attached data and variable dicts
    vrp = prob.vrp
    assign_vars = prob.assign_vars
    tol = prob.tol

    if most:
        best = float('-inf')
    else:
        best = float('inf')
    assign = None
    for i in vrp.EXTLOCS:
        for j in vrp.EXTLOCS:
            if i != j:
                for k in vrp.VEHS:
                    up = ceil(sol[assign_vars[i, j, k]])  # Round up to next nearest integer
                    down = floor(sol[assign_vars[i, j, k]])  # Round down
                    frac = min(up - sol[assign_vars[i, j, k]], sol[assign_vars[i, j, k]] - down)
                    if frac > tol:  # Is fractional?
                        if most and (frac > best):
                            best = frac
                            assign = (i, j, k)
                        elif (not most) and (frac < best):
                            best = frac
                            assign = (i, j, k)

    down_lbs = {}
    down_ubs = {}
    up_lbs = {}
    up_ubs = {}
    if assign is not None:
        #    print assign, sol[assign_vars[assign]]
        down_ubs[assign_vars[assign]] = 0.0
        up_lbs[assign_vars[assign]] = 1.0

        return down_lbs, down_ubs, up_lbs, up_ubs


def symmetry(prob, sol):
    # Get the attached data and variable dicts
    vrp = prob.vrp
    use_vars = prob.use_vars
    tol = prob.tol

    alpha = sum(sol[use_vars[k]] for k in vrp.VEHS)
    #  print "# bins =", alpha
    up = int(ceil(alpha))  # Round up to next nearest integer
    down = int(floor(alpha))  # Round down
    frac = min(up - alpha, alpha - down)
    if frac > tol:  # Is fractional?
        #    print "Symmetry branch"

        down_lbs = {}
        down_ubs = {}
        up_lbs = {}
        up_ubs = {}
        for n in range(up - 1, len(vrp.VEHS)):
            down_ubs[use_vars[vrp.VEHS[n]]] = 0.0
        #    print down_ubs
        for n in range(up):
            up_lbs[use_vars[vrp.VEHS[n]]] = 1.0
        #    print up_lbs

        return down_lbs, down_ubs, up_lbs, up_ubs
