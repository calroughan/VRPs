# Toggle for using networkx (Python's network module) to find cycles
# Might need "pip install networkx" in your virtual environment to use this
USE_NETWORKX = True

# Get a function for help calculating distances
from math import sqrt

# Get the network functionality and
# the drawing functionality
import networkx as nx
import matplotlib.pyplot as plt

# Some options for drawing pictures, be good if this
# gets improved one day
FIGSIZE = (3, 1.5)
FIGSTRETCH = 1.5
NODESIZE = 100  # Default = 300
FONTSIZE = 8  # Default = 12


# The definition of a TSP, including its initialization
# NOTE. This TSP is assumed to be fully connected, i.e.,
# there is a directed arcs between every pair of nodes


class TSPProb:

    def __init__(self, LOCS, x=None, y=None, dist=None):
        # Save the locations
        self.LOCS = LOCS
        # Create a new set of locations with the origin included
        self.EXTLOCS = LOCS[:]
        self.EXTLOCS.append('O')
        # Save the locations of the locations
        self.x = x
        self.y = y
        if (x is None) and (y is None) and (dist is None):
            # Must have x and y coordinates of a distance matrix
            raise Exception("No coordinates or distance matrix in VRPProb!")
        elif (dist is None):
            # If not distance atrix is supplied, then calculate using
            # Euclidean distance between Cartesian coordinates
            dist = {}
            for i in x.keys():
                for j in x.keys():
                    dist[i, j] = sqrt((x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2)
        # Save the distance matrix
        self.dist = dist

    # This function draws the problem, i.e., the locations, it
    # doesn't draw arcs as the network is assumed to be fully connected

    def drawProblem(self):
        # Needs coordinates, i.e., not just distances, to draw the problem
        if (self.x is None) and (self.y is None):
            print("No (x, y)-coordinates so can't draw VRPProb!")
        else:
            # Create a new directed graph
            G = nx.DiGraph()
            # Add all the locations including the origin
            G.add_nodes_from(self.EXTLOCS)
            # Get the positions from the (x, y) coordinates (dictionary)
            pos = dict([(i, (self.x[i], self.y[i])) for i in self.EXTLOCS])
            # Create the figure
            plt.figure(figsize=FIGSIZE)
            # Draw the nodes
            nx.draw(G, pos, with_labels=True, node_size=NODESIZE, font_size=FONTSIZE)
            # Rescale the figure
            xmin = min(self.x.values())
            xmax = max(self.x.values())
            ymin = min(self.y.values())
            ymax = max(self.y.values())
            xscale = FIGSTRETCH * (xmax - xmin)
            yscale = FIGSTRETCH * (ymax - ymin)
            xmid = (xmin + xmax) / 2.0
            ymid = (ymin + ymax) / 2.0
            plt.xlim(xmid - xscale / 2.0, xmid + xscale / 2.0)
            plt.ylim(ymid - yscale / 2.0, ymid + yscale / 2.0)
            # Show the figure
            plt.show()

    # This function sets the current tour for the TSP
    def setSolution(self, assignments, tol):
        # Add all arcs with assignment with value > tol to
        # the current tour for the TSP
        self.assignment = []
        for (i, j) in assignments:
            if assignments[i, j] > tol:
                self.assignment.append((i, j))

    # This function draws the problem (if showProb is True)
    # and the current tour for the TSP
    def displaySolution(self, title=None, showProb=True):
        # Set the arc colour
        color = 'b'
        # Create a new directed graph
        G = nx.DiGraph()
        # Add all the locations including the origin
        G.add_nodes_from(self.EXTLOCS)
        # Get the positions from the (x, y) coordinates (dictionary)
        pos = dict([(i, (self.x[i], self.y[i])) for i in self.EXTLOCS])
        # Get the current figure and resize it
        fig = plt.gcf()
        fig.set_size_inches(FIGSIZE[0], 2 * FIGSIZE[1])
        # Get some values for rescaling
        xmin = min(self.x.values())
        xmax = max(self.x.values())
        ymin = min(self.y.values())
        ymax = max(self.y.values())
        xscale = FIGSTRETCH * (xmax - xmin)
        yscale = FIGSTRETCH * (ymax - ymin)
        xmid = (xmin + xmax) / 2.0
        ymid = (ymin + ymax) / 2.0
        if showProb:
            # If the problem is being shown then...
            ax1 = plt.subplot(211)
            # Draw it with the node labels
            nx.draw_networkx_nodes(G, pos, node_size=NODESIZE)
            nx.draw_networkx_labels(G, pos, font_size=FONTSIZE)
            # And rescale
            ax1.set_xlim(xmid - xscale / 2.0, xmid + xscale / 2.0)
            ax1.set_ylim(ymid - yscale / 2.0, ymid + yscale / 2.0)
            ax1.axis('off')
            # Then move on to the tour plot
            ax2 = plt.subplot(212)
        else:
            # Otherwise just get the tour plot ready
            ax2 = plt.gca()
        # Draw all the arcs in the tour and calculate the distance
        total = 0
        tour = []
        for arc in self.assignment:
            # Print each arc as it is added
            print(" ", arc[0], arc[1], self.dist[arc[0], arc[1]])
            # Ad each arc
            G.add_edge(arc[0], arc[1])
            tour.append((arc[0], arc[1]))
            total += self.dist[arc[0], arc[1]]
        # Print the total distance
        print("  Total =", total)
        #    print(pos)
        #    print(tour)
        # Draw the tour
        nx.draw_networkx_edges(G, pos, edgelist=tour, edge_color=color)
        # Draw the nodes
        nx.draw_networkx_nodes(G, pos, node_size=NODESIZE)
        # Draw the node labels
        nx.draw_networkx_labels(G, pos, font_size=FONTSIZE)
        # Rescale the plot
        ax2.set_xlim(xmid - xscale / 2.0, xmid + xscale / 2.0)
        ax2.set_ylim(ymid - yscale / 2.0, ymid + yscale / 2.0)
        ax2.axis('off')
        if title:
            plt.title(title)
        # Show the plot(s)
        plt.show()


# This function finds a subtour in the graph defined by thenodes and thearcs
# starting from node
def get_subtour_tsp(thenodes, thearcs, node):
    # If using networkx then...
    if USE_NETWORKX:
        # Create a new directed graph
        G = nx.DiGraph()
        # Add the nodes
        G.add_nodes_from(thenodes)
        # Add the arcs
        G.add_edges_from(thearcs)
        # See if a cycle exists starting from node in the (di)graph
        try:
            arcs = nx.find_cycle(G, node)
            # Find the nodes to go with the subtour
            nodes = list(set([a[0] for a in arcs]).union(set([a[1] for a in arcs])))
        except nx.exception.NetworkXNoCycle:
            # If no subtour exists then just return the node
            nodes = [node]
            arcs = []
    else:
        # If not using networkx then...
        # Initialise the set of nodes with the starting node in it
        nodes = set([node])
        # Initialise the set of arcs as empty
        arcs = set([])
        # No nodes have been processed
        not_processed = set(thenodes)
        # The starting node needs to be processed
        to_process = set([node])

        # While any nodes are left to process...
        while to_process:
            # Get a node to process
            c = to_process.pop()
            # It has now been processed
            not_processed.remove(c)
            # Find all arcs from the node being processed
            # to the unprocessed nodes
            new_arcs = [(c, i) \
                        for i in not_processed \
                        if (c, i) in thearcs]
            # In both directions
            new_arcs.extend([(i, c) \
                             for i in not_processed \
                             if (i, c) in thearcs])
            # Find all nodes that are connected to the
            # node being processed
            new_nodes = [i for i in not_processed \
                         if (i, c) in new_arcs]
            new_nodes.extend([i for i in not_processed \
                              if (c, i) in new_arcs])
            # Add the new arcs to the set of arcs in the subtour (|= is union)
            arcs |= set(new_arcs)
            # Add the new nodes to the set of nodes in the subtour
            nodes |= set(new_nodes)
            # Process all the new nodes
            to_process |= set(new_nodes)

    # Convert the sets to lists
    nodes = list(nodes)
    arcs = list(arcs)

    return nodes, arcs
