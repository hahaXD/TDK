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
            g_res[edge] = l.solve_by_random_weights(edge, pair)
        for edge in g["latent_edges"]:
            l_edge = tuple(["latent"] + sorted(edge))
            g_res[l_edge] = l.solve_by_random_weights(l_edge, pair)
        with open("usda/{}_{}.pickle".format(g_id, pair), "wb") as fp:
            pickle.dump(g_res, fp)
    for pair in itertools.combinations(g["latent_edges"], 2):
        print ("Identifying with Constraint ", pair)
        l_pair = [tuple(["latent"] + sorted(pair[0])), tuple(["latent"] + sorted(pair[1]))]
        g_res = {"eql_edges":l_pair}
        for edge in g["arrow_edges"]:
            edge = (edge[0], edge[1])
            g_res[edge] = l.solve_by_random_weights(edge, l_pair)
        for edge in g["latent_edges"]:
            l_edge = tuple(["latent"] + sorted(edge))
            g_res[l_edge] = l.solve_by_random_weights(l_edge, l_pair)
        with open("usda/{}_{}.pickle".format(g_id, l_pair), "wb") as fp:
            pickle.dump(g_res, fp)


