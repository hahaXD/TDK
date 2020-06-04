import pickle
import sys
from sage.all import *
import glob
summary = {}
for i in range(1, 4096):
    fname = "organic/{}.pk".format(i)
    with open(fname, "rb") as fp:
        result = pickle.load(fp)
        for key in result:
            dg = result[key][0][-1].degree(result[key][1])
            if dg not in summary:
                summary[dg] = 1
            else:
                summary[dg] += 1
print (summary)
