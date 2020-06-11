# Get the TSP problem definition with some useful functions
from tsp_prob import TSPProb
# Get the formulation, solve, and assignment extraction functions
from tsp_func import formulate, solve_and_display, get_assignments

if __name__ == '__main__':
    # Python starts here

    # Define the coordinates for the locations to be visited in the TSP problem
    xcoords = [0, 1, 1, \
               4, 4, 4, 5, \
               5, 5]
    ycoords = [4, 2, 4, \
               1, 4, 5, 0, \
               2, 5]
    # Generate a list of locations for the coordinates
    LOCS = list(range(1, len(xcoords) + 1))
    # Generate a dcitionary of (x, y) pairs for each locations
    x = dict([(l, xcoords[i]) for i, l in enumerate(LOCS)])
    y = dict([(l, ycoords[i]) for i, l in enumerate(LOCS)])
    # Define the location of the origin (base) that is not part of LOCS
    x['O'] = 0
    y['O'] = 2

    # Create the TSP problem from the location list and teh x and y dictionaries
    tsp = TSPProb(LOCS=LOCS,
                  x=x,
                  y=y)

    # Draw the TSP
    tsp.drawProblem()

    # Set the zero tolerance
    tol = pow(pow(2, -20), 2.0 / 3.0)
    # Set the solve options including a threshold for considering
    # fractional assignments as part of a (sub)tour
    opts = {
        "Tol": tol,
        "Interval": 1000,
        #      "Tours": tol,
        #      "Tours": 0.5 - tol,
        "Tours": 1.0 - tol,
    }

    # Formulate the TSP
    prob = formulate(tsp, options=opts)

    # Solve the TSP and display the B-&-B tree
    xopt = solve_and_display(prob, options=opts)

    # Extract the arcs assigned to the tour
    assignments = get_assignments(prob, xopt, tol)

    # Set the tour as the current solution
    tsp.setSolution(assignments, tol)

    # Draw the problem and the solution
    tsp.displaySolution(title="Solution")
    # Just draw the solution (tour)
    tsp.displaySolution(title="Solution", showProb=False)
