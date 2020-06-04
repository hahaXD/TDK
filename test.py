import pickle
import sys

with open(sys.argv[1], "rb") as fp:
    result = pickle.load(fp)

for key in result:
    print (key, result[key][0])


