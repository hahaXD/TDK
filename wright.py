# This file implements an algorithm for getting the path equations for various covariances
# in a normalized Linear Acyclic SEM.

# While generating the raw equations is very simple, we want a more clever algorithm which
# directly generates simplified equations.
# In particular, we want the normalization equations to disappear, which can be done in
# acyclic models by simply never expanding past intersections of two paths,
# which is easily achieved by performing simultaneous expansion in the topological sort of
# both nodes.

import networkx as nx
from sympy import symbols, Matrix, eye, zeros


def latentedge(G, edge):
    return tuple(["latent"] + list(G.children(edge[0])))


def generateEdges(G, dname="lambda_", lname="epsilon_"):
    e = {}
    for edge in G.edges():
        if "latent" in G.nodes[edge[0]]:
            le = latentedge(G, edge)
            e[le] = lname + le[1] + le[2]
        else:
            e[edge] = dname + edge[0] + edge[1]
    return e


def ancestralEquation(G, e, n):
    A = {n: 1}
    for parent in G.parents(n):
        if not "latent" in G.nodes[parent]:
            p = e[(parent, n)]
            for key, value in ancestralEquation(G, e, parent).items():
                if key in A:
                    A[key] = A[key] + p * value
                else:
                    A[key] = p * value
    return A


def generateCovarianceValues(G, name="sigma", full=False):
    c = {}
    for node1 in G.nodes():
        if not "latent" in G.nodes[node1]:
            for node2 in G.nodes():
                if not "latent" in G.nodes[node2]:
                    if node1 != node2 and not (node2, node1) in c:
                        c[(node1, node2)] = name + "_" + node1 + node2
                        if full:
                            c[(node2, node1)] = name + "_" + node2 + node1
    return c


def replaceDictVals(d, v, v2):
    d2 = {}
    for key, value in d.items():
        d2[key] = v2[v.index(value)]
    return d2


def getDirectedPaths(G, e):
    directedPaths = {}

    topo_nodes = list(nx.topological_sort(G.nx))

    for node in topo_nodes:
        for ancestor in topo_nodes:
            if ancestor == node:
                break

            if not "latent" in G.nodes[node] and not "latent" in G.nodes[ancestor]:
                # To generate the directed paths, we take the directed paths to all of our parents,
                # and add them up together
                # First generate the directed paths up to this edge
                dp = 0
                for parent in G.parents(node):
                    if not "latent" in G.nodes[parent]:
                        edge = (
                            e[(node, parent)]
                            if (node, parent) in e
                            else e[(parent, node)]
                        )
                        if ancestor == parent:
                            dp += edge
                        elif (ancestor, parent) in directedPaths:
                            dp += edge * directedPaths[(ancestor, parent)]
                directedPaths[(ancestor, node)] = dp

    return directedPaths


def generateCovarianceEquations(G, e):
    directedPaths = {}
    covariances = {}

    topo_nodes = list(nx.topological_sort(G.nx))
    for node in topo_nodes:
        for ancestor in topo_nodes:
            if ancestor == node:
                break

            if not "latent" in G.nodes[node] and not "latent" in G.nodes[ancestor]:
                # To generate the directed paths, we take the directed paths to all of our parents,
                # and add them up together
                # First generate the directed paths up to this edge
                dp = 0
                for parent in G.parents(node):
                    if not "latent" in G.nodes[parent]:
                        edge = (
                            e[(node, parent)]
                            if (node, parent) in e
                            else e[(parent, node)]
                        )
                        if ancestor == parent:
                            dp += edge
                        elif (ancestor, parent) in directedPaths:
                            dp += edge * directedPaths[(ancestor, parent)]
                directedPaths[(ancestor, node)] = dp

                # Next, we want to generate the covariance equations
                # First we take the covariances coming from all parents
                cov = 0
                for parent in G.parents(node):
                    if not "latent" in G.nodes[parent]:
                        edge = (
                            e[(node, parent)]
                            if (node, parent) in e
                            else e[(parent, node)]
                        )
                        if ancestor == parent:
                            cov += edge
                        elif (ancestor, parent) in covariances:
                            cov += edge * covariances[(ancestor, parent)]
                        else:
                            cov += edge * covariances[(parent, ancestor)]
                    else:
                        # Latent variable: we multiply the directed path with the latent var
                        le = latentedge(G, (parent, node))
                        # Get the variable on the other side of the bidirected edge
                        bdnode = le[1]
                        if bdnode == node:  # Wrong parent!
                            bdnode = le[2]
                        if bdnode == ancestor:
                            cov += e[le]
                        else:
                            if (bdnode, ancestor) in directedPaths:
                                cov += e[le] * directedPaths[(bdnode, ancestor)]
                covariances[(ancestor, node)] = cov

    return covariances


def generateSimplifiedCovarianceEquations(G, e, c):
    directedPaths = {}
    covariances = {}

    topo_nodes = list(nx.topological_sort(G.nx))
    for node in topo_nodes:
        for ancestor in topo_nodes:
            if ancestor == node:
                break

            if not "latent" in G.nodes[node] and not "latent" in G.nodes[ancestor]:
                # To generate the directed paths, we take the directed paths to all of our parents,
                # and add them up together
                # First generate the directed paths up to this edge
                dp = 0
                for parent in G.parents(node):
                    if not "latent" in G.nodes[parent]:
                        edge = (
                            e[(node, parent)]
                            if (node, parent) in e
                            else e[(parent, node)]
                        )
                        if ancestor == parent:
                            dp += edge
                        elif (ancestor, parent) in directedPaths:
                            dp += edge * directedPaths[(ancestor, parent)]
                directedPaths[(ancestor, node)] = dp

                # Next, we want to generate the covariance equations
                # First we take the covariances coming from all parents
                cov = 0
                for parent in G.parents(node):
                    if not "latent" in G.nodes[parent]:
                        edge = (
                            e[(node, parent)]
                            if (node, parent) in e
                            else e[(parent, node)]
                        )
                        if ancestor == parent:
                            cov += edge
                        elif (ancestor, parent) in covariances and covariances[
                            (ancestor, parent)
                        ] == 0:
                            pass  # Don't add 0 covariances to the equations
                        elif (parent, ancestor) in covariances and covariances[
                            (parent, ancestor)
                        ] == 0:
                            pass
                        elif (ancestor, parent) in c:
                            cov += edge * c[(ancestor, parent)]
                        else:
                            cov += edge * c[(parent, ancestor)]
                    else:
                        # Latent variable: we multiply the directed path with the latent var
                        le = latentedge(G, (parent, node))
                        # Get the variable on the other side of the bidirected edge
                        bdnode = le[1]
                        if bdnode == node:  # Wrong parent!
                            bdnode = le[2]
                        if bdnode == ancestor:
                            cov += e[le]
                        else:
                            if (bdnode, ancestor) in directedPaths:
                                cov += e[le] * directedPaths[(bdnode, ancestor)]
                covariances[(ancestor, node)] = cov

    return covariances


def generateBilinearCovarianceEquations(G, e, c, d):
    """These covariance equations are split into two parts: The causal
    portion (ie, direct effects) and the observational portion,
    which corresponds to the covariance equations.
    """
    directedPaths = {}
    covariances = {}

    topo_nodes = list(nx.topological_sort(G.nx))
    for node in topo_nodes:
        for ancestor in topo_nodes:
            if ancestor == node:
                break

            if not "latent" in G.nodes[node] and not "latent" in G.nodes[ancestor]:
                # To generate the directed paths, we take the directed paths to all of our parents,
                # and add them up together
                # First generate the directed paths up to this edge
                dp = 0
                for parent in G.parents(node):
                    if not "latent" in G.nodes[parent]:
                        edge = (
                            e[(node, parent)]
                            if (node, parent) in e
                            else e[(parent, node)]
                        )
                        if ancestor == parent:
                            dp += edge
                        elif (ancestor, parent) in directedPaths and directedPaths[
                            (ancestor, parent)
                        ] != 0:
                            dp += edge * d[(ancestor, parent)]
                directedPaths[(ancestor, node)] = dp

                # Next, we want to generate the covariance equations
                # First we take the covariances coming from all parents
                cov = 0
                for parent in G.parents(node):
                    if not "latent" in G.nodes[parent]:
                        edge = (
                            e[(node, parent)]
                            if (node, parent) in e
                            else e[(parent, node)]
                        )
                        if ancestor == parent:
                            cov += edge
                        elif (ancestor, parent) in covariances and covariances[
                            (ancestor, parent)
                        ] == 0:
                            pass  # Don't add 0 covariances to the equations
                        elif (parent, ancestor) in covariances and covariances[
                            (parent, ancestor)
                        ] == 0:
                            pass
                        elif (ancestor, parent) in c:
                            cov += edge * c[(ancestor, parent)]
                        else:
                            cov += edge * c[(parent, ancestor)]
                    else:
                        # Latent variable: we multiply the directed path with the latent var
                        le = latentedge(G, (parent, node))
                        # Get the variable on the other side of the bidirected edge
                        bdnode = le[1]
                        if bdnode == node:  # Wrong parent!
                            bdnode = le[2]
                        if bdnode == ancestor:
                            cov += e[le]
                        else:
                            if (bdnode, ancestor) in directedPaths and directedPaths[
                                (bdnode, ancestor)
                            ] != 0:
                                cov += e[le] * d[(bdnode, ancestor)]
                covariances[(ancestor, node)] = cov

    return (covariances, directedPaths)


def generateFullBilinearCovarianceEquations(G, e, c, d):
    # This version generates all covariance equations that can be written in bilinear form.
    # Usually there are multiple ways to topologically sort the nodes. All this does is gives
    # the equations of all possible topological sorts
    covariances = {}
    # First use a topological sort to generate all directed paths and covariances.
    # We use the old bilinear code to do this
    bicov, directedPaths = generateBilinearCovarianceEquations(G, e, c, d)

    topo_nodes = list(nx.topological_sort(G.nx))
    for node in topo_nodes:
        if "latent" in G.nodes[node]:
            continue
        for ancestor in topo_nodes:
            if node == ancestor:
                continue  # We actually do all other equations
            if "latent" in G.nodes[ancestor]:
                continue
            cov = 0
            if (
                not (node, ancestor) in directedPaths
                or directedPaths[(node, ancestor)] == 0
            ):
                for parent in G.parents(node):
                    if "latent" in G.nodes[parent]:
                        # Latent variable: we multiply the directed path with the latent var
                        le = latentedge(G, (parent, node))
                        # Get the variable on the other side of the bidirected edge
                        bdnode = le[1]
                        if bdnode == node:  # Wrong parent!
                            bdnode = le[2]
                        if bdnode == ancestor:
                            cov += e[le]
                        else:
                            if (bdnode, ancestor) in directedPaths and directedPaths[
                                (bdnode, ancestor)
                            ] != 0:
                                cov += e[le] * d[(bdnode, ancestor)]
                    else:
                        edge = (
                            e[(node, parent)]
                            if (node, parent) in e
                            else e[(parent, node)]
                        )
                        if ancestor == parent:
                            cov += edge
                        elif (ancestor, parent) in bicov and bicov[
                            (ancestor, parent)
                        ] == 0:
                            pass  # Don't add 0 covariances to the equations
                        elif (parent, ancestor) in bicov and bicov[
                            (parent, ancestor)
                        ] == 0:
                            pass
                        elif (ancestor, parent) in c:
                            cov += edge * c[(ancestor, parent)]
                        else:
                            cov += edge * c[(parent, ancestor)]
            if cov != 0:
                covariances[(ancestor, node)] = cov
    return (covariances, directedPaths)


def wright(G, simplify=True):
    e = generateEdges(G)
    c = generateCovarianceValues(G)
    variables = list(e.values()) + list(c.values())
    gens = symbols(variables)

    e = replaceDictVals(e, variables, gens)
    c = replaceDictVals(c, variables, gens)

    cov1 = None
    if simplify:
        cov1 = generateSimplifiedCovarianceEquations(G, e, c)
    else:
        cov1 = generateCovarianceEquations(G, e)

    # Now set the c to have same keys as cov1
    for key in cov1:
        if not key in c:
            c[key] = c[(key[1], key[0])]
            del c[(key[1], key[0])]

    return (cov1, c, e)


def bilinearWright(G, full=False):
    e = generateEdges(G)
    c = generateCovarianceValues(G)
    d = generateCovarianceValues(G, "delta", True)
    variables = list(e.values()) + list(c.values()) + list(d.values())
    gens = symbols(variables)

    e = replaceDictVals(e, variables, gens)
    c = replaceDictVals(c, variables, gens)
    d = replaceDictVals(d, variables, gens)
    cov = None
    de = None
    if full:
        cov, de = generateFullBilinearCovarianceEquations(G, e, c, d)
    else:
        cov, de = generateBilinearCovarianceEquations(G, e, c, d)

    # Now set the c to have same keys as cov1
    for key in cov.keys():
        if not key in c:
            c[key] = c[(key[1], key[0])]
        if cov[key] == 0:
            del cov[key]

    # Remove the 0 directed equations
    for key in de.keys():
        if de[key] == 0:
            del de[key]

    return (cov, de, c, e, d)


def gv(c, n1, n2):
    if (n1, n2) in c:
        return c[(n1, n2)]
    return c[(n2, n1)]


def directedPaths(G, fromnode, tonode):
    e = generateEdges(G)
    variables = list(e.values())
    gens = symbols(variables)
    e = replaceDictVals(e, variables, gens)

    dp = getDirectedPaths(G, e)

    return dp[(fromnode, tonode)], e


def bilinearMatrix(G):
    # Not actually bilinear, it is just the covariance matrix

    # But to generate the matrix, we first need to figure out the number of nodes in the graph, and what they are called.
    n = []

    # While not necessary, it will be useful to see the structure of the DAG reflected in the covariance matrix entries
    topo_nodes = list(nx.topological_sort(G.nx))
    for node in topo_nodes:
        if "latent" in G.nodes[node]:
            continue
        n.append(node)

    # Now generate the variables for the latent variance of each node
    v = []
    for node in n:
        v.append("epsilon_" + node + node)

    # Generate the variables through which the changes will happen
    e = generateEdges(G)
    c = generateCovarianceValues(G)
    d = generateCovarianceValues(G, "delta", True)
    variables = v + list(e.values()) + list(c.values()) + list(d.values())
    gens = symbols(variables)

    e = replaceDictVals(e, variables, gens)
    c = replaceDictVals(c, variables, gens)
    d = replaceDictVals(d, variables, gens)

    for i in range(len(v)):
        v[i] = gens[i]

    # We now have a dict with all the relevant edges as symbols, as well as all the relevant directed edges.
    # We generate the matrix...

    # A lookup for indices
    topo_dict = {}
    for i in range(len(n)):
        topo_dict[n[i]] = i

    Cov = eye(len(n))
    BP = eye(len(n))
    """
    for i in range(len(n)):
        # Diagonal entries are the variance of that specific value, plus 
        # the variances of parents
        BP[i, i] = v[i]
    """

    # We want to generate the directed paths to add to the bilinear forms of the equations
    directedPaths = {}
    for node in topo_nodes:
        if "latent" in G.nodes[node]:
            continue
        for ancestor in topo_nodes:
            if "latent" in G.nodes[ancestor]:
                continue
            if node == ancestor:
                break

            # To generate the directed paths, we take the directed paths to all
            # of our parents, and add them up together:
            dp = 0
            for parent in G.parents(node):
                if "latent" in G.nodes[parent]:
                    continue
                edge = e[(parent, node)]
                if ancestor == parent:
                    dp += edge
                elif (ancestor, parent) in directedPaths:
                    dp += edge * directedPaths[(ancestor, parent)]
            if dp != 0:
                directedPaths[(ancestor, node)] = dp

    # Now go through the covariance matrix one node at a time
    for i in range(len(n)):
        for j in range(len(n)):
            if i == j:
                break
            Cov[i, j] = gv(c, n[i], n[j])

            # Now we go through the parents of j to generate the bilinear forms
            # for parent in G.parents(n[j]):
            node = n[i]
            ancestor = n[j]
            cov = 0
            for parent in G.parents(node):
                if "latent" in G.nodes[parent]:
                    # Latent variable: we multiply the directed path with the latent var
                    le = latentedge(G, (parent, node))
                    # Get the variable on the other side of the bidirected edge
                    bdnode = le[1]
                    if bdnode == node:  # Wrong parent!
                        bdnode = le[2]
                    if bdnode == ancestor:
                        cov += e[le]
                    else:
                        if (bdnode, ancestor) in directedPaths:
                            cov += e[le] * directedPaths[(bdnode, ancestor)]
                else:
                    edge = e[(parent, node)]
                    if ancestor == parent:
                        cov += edge
                    elif (
                        BP[topo_dict[parent], topo_dict[ancestor]]
                        != BP[topo_dict[ancestor], topo_dict[parent]]
                    ):  # check if the covariances are 0 at that element in the matrix
                        cov += edge * gv(c, ancestor, parent)
            BP[i, j] = cov

    # Now that we have one side of the matrix, go through the other side!
    # This time it is more difficult, because if the given node is actually a descendant,
    # then the equations are identical. Only non-descendants can have different equations.
    for j in range(len(n)):
        for i in range(len(n)):
            if i == j:
                break
            Cov[i, j] = gv(c, n[i], n[j])

            # First things first: if other element was 0, this one must be too
            if BP[j, i] == 0:
                continue

            # Next, we check if j is an ancestor of i
            if n[i] in nx.ancestors(G.nx, n[j]):
                BP[i, j] = BP[j, i]
                continue

            # Finally, this means that we have an asymmetric equation (corresponding to a different topological sort)
            # We generate the values here the same way as above
            node = n[i]
            ancestor = n[j]
            cov = 0
            for parent in G.parents(node):
                if "latent" in G.nodes[parent]:
                    # Latent variable: we multiply the directed path with the latent var
                    le = latentedge(G, (parent, node))
                    # Get the variable on the other side of the bidirected edge
                    bdnode = le[1]
                    if bdnode == node:  # Wrong parent!
                        bdnode = le[2]
                    if bdnode == ancestor:
                        cov += e[le]
                    else:
                        if (bdnode, ancestor) in directedPaths:
                            cov += e[le] * directedPaths[(bdnode, ancestor)]
                else:
                    edge = e[(parent, node)]
                    if ancestor == parent:
                        cov += edge
                    elif (
                        BP[topo_dict[parent], topo_dict[ancestor]] != 0
                        or BP[topo_dict[ancestor], topo_dict[parent]] != 0
                    ):  # check if the covariances are 0 at that element in the matrix
                        cov += edge * gv(c, ancestor, parent)
            BP[i, j] = cov

    return BP, Cov, (e, d)


if __name__ == "__main__":
    from sympy import pprint
    from ..linear import Linear

    # The expression can be extracted fro mthe sympy equations as shown here:
    # https://docs.sympy.org/latest/tutorial/manipulation.html
    pprint(directedPaths(Linear("x->y x--y z->x"), "z", "y"))

    # pprint(bilinearMatrix(Linear("1->2 1--2 1->3 1--3 1->4 1--4 1--5 4->5")))

