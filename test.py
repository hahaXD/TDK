import pickle
import sys
from sage.all import *
import glob

fnames = glob.glob("organic/*.pk")
for fname in fnames:
    with open(fname, "rb") as fp:
        result = pickle.load(fp)
        for key in result:
            print (key, result[key][0])


