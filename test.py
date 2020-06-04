import pickle
import sys
from sage.all import *
import glob
summary = {}
for i in range(1, 4096):
    fnames = glob.glob("usda/{}_*.pk".format(i))
    for fname in fnames:
        with open(fname, "rb") as fp:
            result = pickle.load(fp)
            for key in result:
                result[key][0].degree(result[key][1])
                dg = result[key][0][0].degree(result[key][1])
                if dg not in summary:
                    if dg not in summary:
                        summary[dg] = 1
                    else:
                        summary[dg] += 1
print (summary)


