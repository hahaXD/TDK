import networkx as nx

from parser import TMP_parseGraph, TMP_parseEdge


class Graph():
    """The causal graph.

    In causality, each event type is considered a node, and the direction of causation
    is denoted by a directed arrow between the two events. For example, we know that
    the sun rising causes a rooster to crow, and not the other way around::

        import pycausality
        g = pycausality.Graph("sun -> rooster")
        g.show()

    The Graph class is used to hold the causal graph and links to associated datasets.
    It is a subclass of the NetworkX DiGraph, with native support for constructs inherent
    in causal analysis.

    """

    def __init__(self, G, **attr):
        """The graph can be initialized in several ways, which correspond to
        different uses of PyCausality

        When used for exploratory data analysis, PyCausality includes a parser for a simplified
        graph input language built specifically for graphs::

            g = Graph("z->x; x->y; x--y")
            g.show()

        When loading a saved graph from file, you can initialize the Graph::

            g = Graph("./mygraph.txt")



        """
        if isinstance(G, str):
            self.nx = nx.DiGraph(**attr)
            TMP_parseGraph(self.nx, G)
        elif isinstance(G, nx.DiGraph):
            self.nx = G
        elif isinstance(G, Graph):
            self.nx = G.nx
        else:
            self.nx = nx.DiGraph(G, **attr)

    def _parse_edge(self, e):
        return TMP_parseEdge(e)

    @property
    def edges(self):
        """Returns the list of edges in the graph, including edges between latent variables"""
        return self.nx.edges

    @property
    def nodes(self):
        """Returns the list of nodes in the graph, including the latent variables"""
        return self.nx.nodes

    def copy(self):
        """Generates a deep copy of the graph and its properties.
        """
        return Graph(self.nx)

    def save(self, filename):
        """Saves the graph to a file, which can later be loaded with
        the load constructor::

            Graph("x->y").save("myfile.txt")

            g2 = Graph("myfile.txt")

        If you want to export the graph as an image or other file format,
        see the export function.
        """

    def export(self, filename, **attrs):
        """Exports the graph to the given filename. The export format is given
        by the file extension.

        Currently svg and png are supported for drawing the graph,
        with optional width and height::

            Graph("x->y").export("graph.png")

        If you want to save the graph so that it can later be read by PyCausality,
        please use the save function instead.
        """

    def add(self, stuff):
        """Allows you to add to the graph. This function is a catch-all for inserting
        into the graph. You can give it a list of nodes, nodes and edges, or parse strings::

            G.add(['x','y','alpha']) # Adds three nodes
            G.add([('x','y'),('y','x'), 'z']) # Adds 2 edges and a node
            G.add("a->b;c->d;c--d") # Adds nodes and edges in the parse string

        """
        TMP_parseGraph(self.nx, stuff)

    def remove(self, stuff):
        """Enables removing from the graph. This function is a catch-all for removing from
        the graph
        """

    def on_graph_add(self, callback):
        """Allows you to run a callback when the user is attempting to add nodes or edges
        to the graph. 

        Such callbacks can be used to check constraints on the graph. For example::

            def myCallback(G,nodes,edges):
                # TODO: An example
                return (nodes,edges)
            G = Graph()
            G.on_graph_add(myCallback)

        These callbacks enable live editing of the graph such that changes are immediately
        reflected in the editor UI.
        """
        pass

    def on_graph_remove(self, callback):
        """Allows you to run a callback when the user is attempting to remove nodes or edges
        from the graph.

        The callback can be used to check constraints on the graph, and either throw errors,
        or return a modified set of nodes and edges to remove.

        This enables you to handle graphs seen in causality. For example, when implementing
        a graph with transportability nodes, the callback can be used to remove
        the transportability node when the user removes its child.
        """
        pass

    def children(self, node):
        """Returns the children of the given node
        """
        return self.nx.successors(node)

    def parents(self, node):
        """Returns the parents to the given node
        """
        return self.nx.predecessors(node)

    def __str__(self):
        gstring = []
        for node in self.nx:
            if "latent" in self.nx.node[node]:
                gstring.append("--".join(list(self.nx.successors(node))))
            else:
                for succ in self.nx.successors(node):
                    gstring.append(node + "->"+ succ)
        return " ".join(gstring)
    def __repr__(self):
        return "Graph(" + str(self) + ")"
