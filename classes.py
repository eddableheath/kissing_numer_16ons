# Classes for the project
# Author: Edmund Dable-Heath
"""
    Classes for project.
"""

import numpy as np
import algebras as alg

seed = 12345
rng = np.random.default_rng(seed)
machine_zero_threshold = 1e-14


class Sons:

    """
        Main 16-ons class, containing relevant algebraic rules.

        Usage: Give it an 16 dimensional input vector (ndarray) to create 16 on, or let it randomly generate one within
                given bounds, with the optional integer flag allowing for a restriction to the 16 on ring of integers.
    """

    def __init__(self, input_vector=None, bound=16, integer=True):
        if input_vector is not None:
            self.vec = input_vector
        else:
            if integer:
                self.vec = rng.integers(-bound, bound+1, 16)
            else:
                self.vec = rng.normal(0, bound, 16)
            while np.sum(self.vec[8:]) < machine_zero_threshold:
                if integer:
                    self.vec = rng.integers(-bound, bound+1, 16)
                else:
                    self.vec = rng.normal(0, bound, 16)

    def mult(self, x):
        """
            Multiply with other 16on x
        """
        a, b = self.vec[:8], self.vec[8:]         #Split 16ons into pairs of octonions
        c, d = x.vec[:8], x.vec[8:]
        z = np.zeros(16)
        z[:8] = alg.omult(a, c) - alg.omult(d, alg.ostar(b))
        z[8:] = alg.omult(c, b) + alg.omult(alg.omult(alg.ostar(a), alg.oinverse(b)), alg.omult(b, d))
        return Sons(input_vector=z)

    def smult(self, s):
        """
            Scalar multiply
        """
        return Sons(s * self.vec)

    def add(self, x):
        """
            Add another 16on x
        """
        return Sons(self.vec + x.vec)

    def conj(self):
        """
            16on conjugation
        """
        mask = -np.ones(16)
        mask[0] = 1
        return Sons(self.vec * mask)

    def norm(self):
        """
            Norm of element
        """
        return self.mult(self.conj()).vec[0]

    def inverse(self):
        """
            Inverse of element
        """
        return Sons(self.conj().vec / self.norm())


