# HACK: This is here temporarily to have something working. This needs to be refactored ASAP
# Manually inputting the graph is a task I heavily dislike,
# so I am including a simple parser for a simplified version of
# the text graph input used in fusion.
# That is, we simply draw the arrows:
#
# z->x
# x->y
# x--y
#
# The above 3 lines represent the instrumental variable.
import pyparsing as pp
# Set up the variable names - var represents a node name
varname = pp.Combine(
    pp.Word(pp.alphanums + "_", exact=1) + pp.Optional(  #Allow number nodes
        pp.Word(pp.alphanums + "_")))
arrow = pp.Or(["--", "->"])
edge = pp.Group(varname + arrow + varname)
graphparser = pp.OneOrMore(edge)


def TMP_parseGraph(G, txt):
    parseresult = graphparser.parseString(txt)

    for edge in parseresult:
        if edge[1] == "->":
            G.add_edge(str(edge[0]), str(edge[2]))
        else:
            # Uh oh, latent alert!
            latentName = "U({},{})".format(edge[0], edge[2])
            G.add_edges_from([(latentName, str(edge[0])), (latentName,
                                                           str(edge[2]))])
            G._node[latentName]["latent"] = True
    return G


def TMP_parseEdge(e):
    if len(e) == 2 or len(e) == 3:
        return (str(e[0]), str(e[1]))
    else:
        e = edge.parseString(e)[0]
        if e[1] == "->":
            return (str(e[0]), str(e[2]))
        else:
            return e
