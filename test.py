# from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
# from sage.rings.rational_field import QQ
# from sage.rings.ideal import Ideal
# from sage.symbolic.ring import var
from sage.all import *

variables = ["x", "y"]
ring = PolynomialRing(QQ, variables)
gens = ring.gens()
Ideal([gens[0]**2 - gens[1], gens[1]*gens[0]]).groebner_basis().solve()
