# Algebraic rules
# Author: Edmund Dable-Heath
"""
     Quaternion and up algebraic rules:
"""

import numpy as np


# Quaternions
def qmult(x: np.array, y: np.array) -> np.array:
    """
        Quaternion multiplication
    """
    return np.array([
        x[0]*y[0] - x[1]*y[1] - x[2]*y[2] - x[3]*y[3],
        x[0]*y[1] + x[1]*y[0] + x[2]*y[3] - x[3]*y[2],
        x[0]*y[2] - x[1]*y[3] + x[2]*y[0] + x[3]*y[1],
        x[0]*y[3] + x[1]*y[2] - x[2]*y[1] + x[3]*y[0]
    ])


def qstar(x: np.array) -> np.array:
    """
        Quaternion conjugate
    """
    return x*np.array([1, -1, -1, -1])


def qnorm(x: np.array) -> np.array:
    """
        Quaternion norm
    """
    return qmult(x, qstar(x))


def qinverse(x: np.array) -> np.array:
    """
        Quaternion inverse
    """
    return qstar(x)/qnorm(x)[0]


# Octonions
def omult(x: np.array, y: np.array) -> np.array:
    """
        Octonion multiplication
    """
    a, b = x[:4], x[4:]     #S plit octonions into pairs of quaternions
    c, d = y[:4], y[4:]

    z = np.zeros(8)
    z[:4] = qmult(a, c) - qmult(d, qstar(b))
    z[4:] = qmult(qstar(a), d) + qmult(c, b)
    return z


def ostar(x: np.array) -> np.array:
    """
        Octonion conjugate
    """
    mask = -np.ones(8)
    mask[0] = 1
    return x*mask


def onorm(x: np.array) -> np.array:
    """
        Octonion norm
    """
    return omult(x, ostar(x))


def oinverse(x: np.array) -> np.array:
    """
        Octonion inverse
    """
    return ostar(x)/onorm(x)[0]