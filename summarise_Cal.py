from random import random, seed
import time

from veh_rout_prob import VRProb

# No need to import your veh_rout_func, that is done automatically based on your UPI. (Ignore the red function warnings,
# it "should" run

# Leave asgn001 here - It is Mike's benchmark code posted to Canvas
# Change the one below it to your UPI
upis = [
    "asgn001",
    "crou060"
]


def total_dist(vrp):
    total = 0
    for k in vrp.VEHS:
        for arc in vrp.assignment[k]:
            print(" ", arc[0], arc[1], vrp.dist[arc[0], arc[1]])
            total += vrp.dist[arc[0], arc[1]]

    return total


def cart1(theseed, disp=False):
    seed(theseed)

    # Python starts here
    numLocs = 5
    locs = list(range(1, numLocs + 1))
    itmBnds = [0, 4]

    x = dict([(i, round(random() * (itmBnds[1] - itmBnds[0]) + itmBnds[0], 2)) for i in locs])
    y = dict([(i, round(random() * (itmBnds[1] - itmBnds[0]) + itmBnds[0], 2)) for i in locs])
    x['O'] = 2
    y['O'] = 2

    ncurr = 5

    vrp = VRProb(LOCS=locs,
                 ncurr=ncurr,
                 x=x,
                 y=y)

    opts = myopts

    prob = formulate(vrp, options=opts)

    tstart = time.time()
    xopt = solve(prob, options=opts)
    dur = time.time() - tstart

    assignments = get_assignments(prob, xopt, opts["Tol"])

    vrp.setSolution(assignments, prob.tol)

    if disp:
        vrp.displaySolution(title="Solution, veh_rout_cart1, seed = %d" % theseed, showProb=False)
    tdist = total_dist(vrp)
    if disp: print("Overall distance", tdist)

    tree_nodes = prob.Tree.get_node_list()
    numNodes = len(tree_nodes)

    return tdist, numNodes, dur


def cart2(theseed, disp=True):
    seed(theseed)

    # Python starts here
    numLocs = 5
    locs = list(range(1, numLocs + 1))
    itmBnds = [0, 4]

    x = dict([(i, round(random() * (itmBnds[1] - itmBnds[0]) + itmBnds[0], 2)) for i in locs])
    y = dict([(i, round(random() * (itmBnds[1] - itmBnds[0]) + itmBnds[0], 2)) for i in locs])
    x['O'] = 2
    y['O'] = 2

    ncurr = 5

    vrp = VRProb(LOCS=locs,
                 ncurr=ncurr,
                 x=x,
                 y=y,
                 maxdist=6)

    opts = myopts

    prob = formulate(vrp, options=opts)

    tstart = time.time()
    xopt = solve(prob, options=opts)
    dur = time.time() - tstart

    assignments = get_assignments(prob, xopt, opts["Tol"])

    vrp.setSolution(assignments, prob.tol)

    if disp:
        vrp.displaySolution(title="Solution, veh_rout_cart1, seed = %d" % theseed, showProb=False)
    tdist = total_dist(vrp)
    if disp: print("Overall distance", tdist)

    tree_nodes = prob.Tree.get_node_list()
    numNodes = len(tree_nodes)

    return tdist, numNodes, dur


def cart3(theseed, disp=True):
    seed(theseed)

    # Python starts here
    numLocs = 5
    locs = list(range(1, numLocs + 1))
    itmBnds = [0, 4]

    x = dict([(i, round(random() * (itmBnds[1] - itmBnds[0]) + itmBnds[0], 2)) for i in locs])
    y = dict([(i, round(random() * (itmBnds[1] - itmBnds[0]) + itmBnds[0], 2)) for i in locs])
    x['O'] = 2
    y['O'] = 2

    ncurr = 3  # <--- NOTE less vehicles

    vrp = VRProb(LOCS=locs,
                 ncurr=ncurr,
                 x=x,
                 y=y,
                 useall=True)

    opts = myopts

    prob = formulate(vrp, options=opts)

    tstart = time.time()
    xopt = solve(prob, options=opts)
    dur = time.time() - tstart

    assignments = get_assignments(prob, xopt, opts["Tol"])

    vrp.setSolution(assignments, prob.tol)

    if disp:
        vrp.displaySolution(title="Solution, veh_rout_cart1, seed = %d" % theseed, showProb=False)
    tdist = total_dist(vrp)
    if disp: print("Overall distance", tdist)

    tree_nodes = prob.Tree.get_node_list()
    numNodes = len(tree_nodes)

    return tdist, numNodes, dur


results = {}
for upi in upis:
    try:
        exec("from %s_veh_rout_func import myopts, formulate, solve, get_assignments" % upi)
        myopts['Interval'] = None

        seeds = {
            # 1: [1, 5],
            # 2: [1, 3, 5, 7],
            # 3: [1, 5, 9, 13]
            1: [1, 5, 23, 42, 64, 100],
            2: [1, 5, 23, 42, 64, 100],
            3: [1, 5, 23, 42, 64, 100]
        }

        result = {}
        for inst in seeds:
            for s in seeds[inst]:
                if inst == 1:
                    result['cart%d_%d' % (inst, s)] = cart1(s, disp=False)
                elif inst == 2:
                    result['cart%d_%d' % (inst, s)] = cart2(s, disp=False)
                elif inst == 3:
                    result['cart%d_%d' % (inst, s)] = cart3(s, disp=False)
        results[upi] = result

    except:
        results[upi] = None

print("\n\n\n====================================\n")
titles = ['UPI:',
          'Cart_&_Seed:',
          'Total_distance:',
          'Nodes_in_tree:',
          'Time_taken_(s):']
first = []
second = []

for p in results[upis[0]]:

    try:

        first.append([upis[0],
                      p[4:],
                      round(results[upis[0]][p][0], 3),
                      results[upis[0]][p][1],
                      round(results[upis[0]][p][2], 3)
                      ])

    except:
        print(upis[0], " threw an error on cart and seed: ", p[4:], )
        first.append([upis[0], "Failed", "Failed", "Failed", "Failed"])

    try:
        second.append([upis[1],
                       p[4:],
                       round(results[upis[1]][p][0], 3),
                       results[upis[1]][p][1],
                       round(results[upis[1]][p][2], 3)
                       ])
    except:
        print(upis[1], " threw an error on cart and seed: ", p[4:], )
        second.append([upis[1], "Failed", "Failed", "Failed", "Failed"])


print("\n\nHere's how your code compared to the benchmark asgn001:\n")
for row in range(0, len(first)):
    for item in range(0, 5):
        print("".join(titles[item].ljust(20)),
              "".join(str(first[row][item]).ljust(12)),
              "".join(str(second[row][item]).ljust(12)))
    print("\n")

