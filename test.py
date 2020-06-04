import pickle
import sys
from sage.all import *
import glob

fnames = glob.glob("usda/*.pk")
for fname in fnames:
    with open(fname, "rb") as fp:
        result = pickle.load(fp)
        for key in result:
            print (key, len(result[key][0]), result[key][0])
            if len(result[key][0]) > 1:
                break


