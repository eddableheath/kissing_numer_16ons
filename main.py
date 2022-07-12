# Main Code Playground
# Author: Edmund Dable-Heath
"""
    Based on kissing number computation from 16ons work with Andrew Mendelsohn. Looking to iteratively construct a
    subset A of the 16ons that maximises the energy and minimises the size of A+A.
"""

from classes import Sons


# testing/run
if __name__ == "__main__":
    x = Sons()
    print(x)
    print('.vec gives vector representation')
    print(x.vec)
    x_inv = x.inverse()
    print('creates new 16on with inverse')
    print(x_inv.vec)
    print('re multiplying should give identity')
    print(x.mult(x_inv).vec)
    x_conj = x.conj()
    print('creates new 16on as conjugate')
    print(x_conj.vec)
    print('finds the norm of the element multiplied with its own inverse')
    print(x.mult(x.inverse()).norm())
    y = Sons()
    z = x.smult(-5)
    print(z.vec)
