# common.py

from functools import lru_cache

import numpy as np
from scipy.integrate import quad
from scipy.special import ellipe as eE
from scipy.special import ellipk as eK
from scipy.special import j0, j1, struve

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


class drivingCoil():
    def __init__(self, rdi, rde, ld, n, resistivity, Swire, k):
        self.ri = rdi
        self.re = rde
        self.l = ld
        self.n = n                      # 线圈匝数
        self.x = 0.5 * self.l
        self.SR = resistivity           # 电阻率
        self.Swire = Swire              # 单根导线的截面积
        self.k = k                      # 驱动线圈填充率

        self.r = (self.ri + self.re) / 2
        Sd = (self.re - self.ri) * self.l

        self.nc = self.n / Sd

        self.R = self.R()
        self.L = L(self)

    def R(self):
        return (self.SR * self.k * np.pi * (np.power(self.re, 2) - np.power(self.ri, 2)) * self.l) / np.power(self.Swire, 2)


class armature():
    def __init__(self, rai, rae, la, resistivity, v0, ma, m, n, x0):
        self.ri = rai
        self.re = rae
        self.l = la
        self.SR = resistivity
        self.v0 = v0
        self.ma = ma
        self.m = m
        self.n = n

        self.x = x0

        self.currentFilaments = {}

        for i in range(1, self.m + 1):
            for j in range(1, self.n + 1):
                self.currentFilaments[i][j] = currentFilament(ri=self.currentFilamentR(j - 1), re=self.currentFilamentR(j),
                                                              l=self.l / self.m, R=None, L=None, x0=self.l * (i - 0.5) / self.m)

        self.R()
        self.L()

    def currentFilamentR(self, j):
        return self.ri + (self.re - self.ri) * j / self.n

    def updatePosition(self, delta):
        for i in range(1, self.m + 1):
            for j in range(1, self.n + 1):
                self.currentFilaments[i][j].updatePosition(delta)

    def R(self):
        deltaR = 2 * np.pi * self.SR * self.m / self.l
        R = {
            1: 2 * np.pi * self.SR * ((self.m / (2 * self.l)) + (self.m * self.n * self.ri / (self.l * (self.re - self.ri))))
        }

        for a in range(2, self.n + 1):
            R[a] = R[a - 1] + deltaR

        for i in range(1, self.m + 1):
            for j in range(1, self.n + 1):
                self.currentFilaments[i][j].R = R[j]

    def L(self):
        for i in range(1, self.m + 1):
            for j in range(1, self.n + 1):
                self.currentFilaments[i][j].L = L(self.currentFilaments[i][j])
