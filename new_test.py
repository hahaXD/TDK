import itertools
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
    for pair in itertools.combinations(g["arrow_edges"], 2):
        print ("Identifying with Constraint ", pair)
        pair = [(pair[0][0],pair[0][1]),  (pair[1][0],pair[1][1])]
        g_res = {"eql_edges":pair}
        for edge in g["arrow_edges"]:
            edge = (edge[0], edge[1])
            l.g_relevent_basis(edge, pair)
        for edge in g["latent_edges"]:
            l_edge = tuple(["latent"] + sorted(edge))
            l.g_relevent_basis(l_edge, pair)
        with open("usda_basic/{}_{}.pickle".format(g_id, pair), "wb") as fp:
            pickle.dump(g_res, fp)
    for pair in itertools.combinations(g["latent_edges"], 2):
        l_pair = [tuple(["latent"] + sorted(pair[0])), tuple(["latent"] + sorted(pair[1]))]
        print ("Identifying with Constraint ", l_pair)
        g_res = {"eql_edges":l_pair}
        for edge in g["arrow_edges"]:
            edge = (edge[0], edge[1])
            l.g_relevent_basis(edge, l_pair)
        for edge in g["latent_edges"]:
            l_edge = tuple(["latent"] + sorted(edge))
            l.g_relevent_basis(l_edge, l_pair)
        with open("usda_basic/{}_{}.pickle".format(g_id, l_pair), "wb") as fp:
            pickle.dump(g_res, fp)


