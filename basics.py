# common.py


import numpy as np

from scipy.integrate import quad
from scipy.special import j0, j1, struve
from scipy.special import ellipe as eE
from scipy.special import ellipk as eK

from functools import lru_cache


μ0 = 4 * np.pi * np.power(10, -7)


def L(coilA, limit=125):
    p = coilA.re / coilA.ri
    q = coilA.l / coilA.ri

    # U = lambda x: quad(lambda x: x * j1(x), x, p * x)[0] / np.power(x, 3)
    def U(x): return np.pi * (-j1(x) * struve(0, x) + p * j1(p * x) * struve(0, p * x) +
                              j0(x) * struve(1, x) - p * j0(p * x) * struve(1, p * x)) / (2 * np.power(x, 2))

    def integrationT(x): return np.power(U(x), 2) * \
        (q * x + np.power(np.e, (-1 * q * x)) - 1)
    T = quad(integrationT, 0, np.inf, limit=limit)[0]

    return 2 * np.pi * μ0 * np.power(coilA.nc, 2) * np.power(coilA.ri, 5) * T


def calcK(Ra, Rb, d):
    return np.sqrt((4 * Ra * Rb) / (np.power((Ra + Rb), 2) + np.power(d, 2)))


@lru_cache()
def M(Ra, Rb, d):
    k = calcK(Ra, Rb, d)

    return μ0 * np.sqrt(Ra * Rb) * ((2 / k - k) * eK(k) - (2 / k) * eE(k))


@lru_cache()
def dM(Ra, Rb, d):
    k = calcK(Ra, Rb, d)

    return (μ0 * k * d * (2 * (1 - np.power(k, 2)) * eK(k) - (2 - np.power(k, 2)) * eE(k))) / (4 * (1 - np.power(k, 2)) * np.sqrt(Ra * Rb))


class currentFilament():
    def __init__(self, ri, re, l, R, L, x0):
        self.ri = ri
        self.re = re
        self.l = l
        self.R = R
        self.L = L
        self.x = x0

        self.r = (self.ri + self.re) / 2

        Sc = l * (self.re - self.ri)
        self.nc = 1 / Sc    

    def updatePosition(self, delta):
        self.x += delta
