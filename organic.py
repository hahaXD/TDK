import linear as lr
import sys
import json
with open (sys.argv[1], "r") as fp:
    graphs = json.load(fp)["graphs"]

n_graphs = len(graphs)
results = []
for g_id, g in enumerate(graphs):
    print ("Handeling {} / {} graphs".format(g_id, n_graphs))
    if len(g["arrow_edges"]) == 0 and len(g["latent_edges"]) == 0:
        results.append({"graph": g})
        continue
    l = lr.Linear(g["graph_str"])
    print (g["graph_str"])
    for edge in g["arrow_edges"]:
        basis, target, e = l.gidentify((edge[0], edge[1]))
        if len(basis) > 0:
            print (basis)
            import sys
            sys.exit(1)
