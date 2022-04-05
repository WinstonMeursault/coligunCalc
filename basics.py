# common.py


import numpy as np

from scipy.integrate import quad
from scipy.special import j0, j1, struve
from scipy.special import ellipe as eE
from scipy.special import ellipk as eK


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


def calcK(coilA, coilB, d):
    return np.sqrt((4 * coilA.r * coilB.r) / (np.power((coilA.r + coilB.r), 2) + np.power(d, 2)))


def M(coilA, coilB):
    d = abs(coilA.x - coilB.x)
    k = calcK(coilA, coilB, d)

    return μ0 * np.sqrt(coilA.r * coilB.r) * ((2 / k - k) * eK(k) - (2 / k) * eE(k))


def dM(coilA, coilB, d=None):
    d = abs(coilA.x - coilB.x)
    k = calcK(coilA, coilB, d)

    return (μ0 * k * d * (2 * (1 - np.power(k, 2)) * eK(k) - (2 - np.power(k, 2)) * eE(k))) / (4 * (1 - np.power(k, 2)) * np.sqrt(coilA.r * coilB.r))


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
