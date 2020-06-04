from graph import Graph
from wright import generateCovarianceEquations, generateSimplifiedCovarianceEquations
import networkx as nx
import random as rnd

def edges_union_find(eql_edges):
    uf = nx.utils.UnionFind()
    for edge_pair in eql_edges:
        uf.union(edge_pair[0], edge_pair[1])
    return uf

def latentedge(G, edge):
    q = sorted(list(G.children(edge[0])))
    return tuple(['latent'] + q)


def generateEdges(G, dname="lambda_", lname="epsilon_"):
    e = {}
    for edge in G.edges():
        if "latent" in G.nodes[edge[0]]:
            le = latentedge(G, edge)
            e[le] = lname + le[1] + le[2]
        else:
            e[edge] = dname + edge[0] + edge[1]
    return e

def replaceDictVals(d, v, v2):
    d2 = {}
    for key, value in d.items():
        d2[key] = v2[v.index(value)]
    return d2

class Linear(Graph):
    """Linear Structural Equation Model.

    This is a parametric model, where it is assumed that
    the output of a node is a linear combination of its inputs. Each edge is associated with what is
    called a "Structural Parameter", which represents the weight that the tail node of the edge has
    on the head node, meaning that the structural parameter is the :math:`\\alpha` in :math:`Y=\\alpha X`.
    """

    def __init__(self, G):
        Graph.__init__(self, G)

    def solve_by_random_weights(self, idedge, eql_edge):

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.rational_field import QQ
        from sage.rings.ideal import Ideal
        from sage.symbolic.ring import var
        e  = generateEdges(self)
        e2 = generateEdges(self, "X_", "e_")
        variables = list(e2.values())
        ring = PolynomialRing(QQ, variables)
        gens = ring.gens()
        for e_t in e:
            e[e_t] = (rnd.random() - 0.5) * 10
        e[eql_edge[0]] = (rnd.random() - 0.5) * 10
        e[eql_edge[1]] = e[eql_edge[0]]
        e2 = replaceDictVals(e2, variables, gens)
        cov1 = generateCovarianceEquations(self, e)
        cov2 = generateCovarianceEquations(self, e2)
        eqns = [cov1[k] - cov2[k] for k in cov1]
        eqns.append(e2[eql_edge[0]] - e2[eql_edge[1]])
        dc_vars  = [e2[k] for k in e2 if k != idedge]
        core_basis = Ideal(eqns).elimination_ideal(dc_vars).groebner_basis()
        return core_basis, e2[idedge]

    def gidentify(self, idedge, eql_edges=None):
        import sage.all
        if eql_edges is not None:
            eql_edges_uf = edges_union_find (eql_edges)
        else:
            eql_edges_uf = []
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.rational_field import QQ
        from sage.rings.ideal import Ideal
        from sage.symbolic.ring import var
        e  = generateEdges(self)
        e2 = generateEdges(self, "X_", "e_")
        variables = list(e2.values()) + list(e.values())
        ring = PolynomialRing(QQ, variables)
        gens = ring.gens()
        e = replaceDictVals(e, variables, gens)
        e2 = replaceDictVals(e2, variables, gens)
        if eql_edges is not None:
            for eql_edge in eql_edges:
                e[eql_edge[1]] = e[eql_edges_uf[eql_edge[1]]]
                e[eql_edge[0]] = e[eql_edges_uf[eql_edge[0]]]
        cov1 = generateCovarianceEquations(self, e)
        cov2 = generateCovarianceEquations(self, e2)
        eqns = [cov1[k] - cov2[k] for k in cov1]
        dc_vars  = [e2[k] for k in e2 if k != idedge and k not in eql_edges_uf]
        core_basis = Ideal(eqns).elimination_ideal(dc_vars).groebner_basis()
        if eql_edges is not None:
            eql_vars = [e2[k] for k in e2 if k != idedge and k in eql_edges_uf]
            core_eqns = [k for k in core_basis]
            for eql_edge in eql_edges:
                core_eqns.append(e2[eql_edge[0]] - e2[eql_edge[1]])
            core_basis = Ideal(core_eqns).elimination_ideal(eql_vars).groebner_basis()
        return core_basis, e2[idedge], e
        #print (eqlvars, solvevars)
        #print (myideal.groebner_basis())
        #return (myideal.groebner_basis(), e2[idedge], e, myideal)
        # # Now replace the edge maps with the actual variables:
        # gens = ring.gens()
        # e = replaceDictVals(e, variables, gens)
        # e2 = replaceDictVals(e2, variables, gens)

        # # Finally, get the covariance equations from wright's rules
        # # for both sets of variables

        # cov1 = generateCovarianceEquations(G, e)
        # cov2 = generateCovarianceEquations(G, e2)

        # eqns = [cov1[k] - cov2[k] for k in cov1]
        # knownedges = [getedge(k, e2) for k in known]
        # solvevars = [e2[k] for k in e2 if k != idedge and e2[k] not in knownedges]
        # myideal = Ideal(eqns).elimination_ideal(solvevars)
        # return (myideal.groebner_basis(), e2[idedge], e, myideal)

        # e  = generateEdges(self)
        # e2 = generateEdges(self, "X_", "e_")
        # v_names = list(e.values()) + list(e2.values())
        # sym_v   = [sp.symbols(nm) for nm in v_names]
        # target_v = [sym_v[i] for i in range(0, len(v_names)) if "lambda_" in v_names[i] or "epsilon_" in v_names[i]]
        # e  = replaceDictVals(e, v_names, sym_v)
        # e2 = replaceDictVals(e2, v_names, sym_v)
        # cov1 = generateCovarianceEquations(self, e)
        # cov2 = generateCovarianceEquations(self, e2)
        # eqns = [cov1[k] - cov2[k] for k in cov1]
        # print (target_v)
        # print (sp.solve(eqns, target_v))



    def identifiable_edges(self, known=[]):
        """Returns a list of the edges that can be identified in the current graph
        using auxiliary variables :cite:`chen2017identification`. This method is not necessary for identification,
        so there might be some edges which are identifiable but not returned by this method.

        Returns:
            list: A list of identifiable edges in the graph
        """
        return getIdentifiable(self.nx, [self._parse_edge(i) for i in known])

    def paths(self, simplify=False, sage=False):
        """Returns the equations corresponding to wright's rules for the graph. The equations
        are in terms of sympy variables.


        Args:
            simplify (bool,optional): When set to True, returns pre-simplified equations which 
                are most amenable to computation. Rather than being the raw wright's rules, the
                equations substitute in the covariance between other variables to simplify computation
                where possible.
        Returns:
            paths (dict): The path equations
            covariances (dict): The corresponding covariance variables
            edges (dict): The dict of the variables corresponding to edges.
        """
        from .algorithms.wright import wright
        p, c, e = wright(self, simplify)
        if not sage:
            return (p, c, e)
        else:
            e_new = {}
            for key in e:
                e_new[key] = e[key]._sage_()
            c_new = {}
            for key in c:
                c_new[key] = c[key]._sage_()
            p_new = {}
            for key in p:
                try:
                    p_new[key] = p[key]._sage_()
                except:
                    p_new[key] = p[key]
            return (p_new, c_new, e_new)

    def bilinear_paths(self, sage=False, full=False):
        from .algorithms.wright import bilinearWright
        p, de, c, e, d = bilinearWright(self, full)
        if not sage:
            return (p, de, c, e, d)
        else:
            e_new = dconv(e)
            c_new = dconv(c)
            p_new = dconv(p)
            d_new = dconv(de)
            de_new = dconv(de)
            return (p_new, de_new, c_new, e_new, d_new)

    def bilinear_matrix(self, sage=False):
        from .algorithms.wright import bilinearMatrix
        v, c, d = bilinearMatrix(self)
        if not sage:
            return v, c, d

        # SAGE requires... more stuff
        from sage.all import SR, matrix
        l = v.shape[0]
        v_ = []
        c_ = []
        for i in range(l * l):
            v_.append(v[i]._sage_())
            c_.append(c[i]._sage_())

        return matrix(SR, l, l, v_), matrix(SR, l, l, c_), (dconv(d[0]),
                                                            dconv(d[1]))

    def identify(self, e, known=[], sage=True, full=False):
        a, b = self._parse_edge(e)
        val = ideqn(self.nx, a, b, [self._parse_edge(i) for i in known])
        if val is None:
            return None
        eqns, system = val
        if not full:
            if not sage:
                return eqns[(a, b)][2]
            return eqns[(a, b)][2]._sage_()

        if not sage:
            return eqns, system
        for i in eqns:
            curval = eqns[i]

            eqns[i] = (curval[0]._sage_(), curval[1]._sage_(),
                       curval[2]._sage_())
        return eqns, system

    def groebner(self, edge, format="covariance", known=[]):
        """Returns the Groebner basis of the solutions for the given edge.

        It performs variable elimination as detailed in :cite:`garcia2010identifying`. A set of equations
        in the goal variable is returned, whose solutions will be the possible values of the structural parameter consistent with
        an observed covariance matrix.

        If all of the bases returned are 0, it means that the structural parameter is not identifiable.

        .. note:: This code relies on SAGE_. If you have SAGE_ installed,
            you might be able to get it working in python, but there are issues importing SAGE_ modules
            on some systems (such as unset environmental variables).
            If you run into issues, it is recommended to run this code directly from the SAGE_ app or notebook.

        .. _SAGE: http://www.sagemath.org/

        Args:
            edge: The edge to generate a Groebner basis for.
            format (str,optional): The type of basis to return. There are two possible values,
                "structural" and "covariance". "covariance" returns the basis in terms of the
                covariance, which can be used to solve for the parameter from the observational
                distribution (this is the algorithm from the paper). The "structural" solves for
                the parameter in terms of the parameter, ie, given the structural parameters of the graph,
                it solves for the possible values of the graph from the covariance, but in terms of the parameters.
        Returns:
            The groebner basis for the solutions. Each equation is in one variable (X), and setting all of them equal to 0 gives
            the possible values of the edge consistent with the covariance.

        """
        from .algorithms.groebner import gidentify, gcidentify, gsidentify
        e = self._parse_edge(edge)
        known = [self._parse_edge(i) for i in known]
        if format == "covariance":
            return gcidentify(self, e, True, known=known)
        elif format == "covariance-full":
            return gcidentify(self, e, False, known=known)
        elif format == "simplified":
            return gsidentify(self, e, known=known)
        return gidentify(self, e, known=known)

    def __repr__(self):
        return "Linear(" + str(self) + ")"
