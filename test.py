# from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
# from sage.rings.rational_field import QQ
# from sage.rings.ideal import Ideal
# from sage.symbolic.ring import var
from sage.all import *
import sympy as sp

variables = ["x", "y"]
ring = PolynomialRing(QQ, variables)
gens = ring.gens()
B = Ideal([gens[0]**2 - gens[1], gens[1]*gens[0]]).groebner_basis()
k = "["
q ="["+",".join([ str(p) + "=0" for p in B]) + "]"
print (sage.symbolic.relation.string_to_list_of_solutions(q))


