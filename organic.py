import linear as lr
import sys
import json
import pickle
with open (sys.argv[1], "r") as fp:
    graphs = json.load(fp)["graphs"]

n_graphs = len(graphs)
for g_id, g in enumerate(graphs):
    print ("Handeling {} / {} graphs".format(g_id, n_graphs))
    if len(g["arrow_edges"]) == 0 and len(g["latent_edges"]) == 0:
        continue
    l = lr.Linear(g["graph_str"])
    print (g["graph_str"])
    g_res = {}
    for edge in g["arrow_edges"]:
        basis, target, e = l.gidentify((edge[0], edge[1]))
        g_res[(edge[0], edge[1])] = (basis, target, e)
        if len(basis) > 0:
            assert (basis[-1].degree(target) != 0)
    for edge in g["latent_edges"]:
        l_edge = tuple(["latent"] + sorted(edge))
        basis, target, e = l.gidentify(l_edge)
        g_res[l_edge] = (basis, target, e)
        if len(basis) > 0:
            assert (basis[-1].degree(target) != 0)
    with ("organic/{}.pk".format(g_id), "w") as fp:
        pickle.dump(g_res, fp)
