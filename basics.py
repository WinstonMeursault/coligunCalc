# common.py

from functools import lru_cache

import numpy as np
from scipy.integrate import quad
from scipy.special import ellipe as eE
from scipy.special import ellipk as eK
from scipy.special import j0, j1, struve

μ0 = 4 * np.pi * np.float_power(10, -7)


def calcL(ri, re, l, nc, limit=200):
    p = re / ri
    q = l / ri

    # U = lambda x: quad(lambda x: x * j1(x), x, p * x)[0] / np.power(x, 3)
    def U(x):
        return np.pi * (-j1(x) * struve(0, x) + p * j1(p * x) * struve(0, p * x) + j0(x) * struve(1, x) - p * j0(p * x) * struve(1, p * x)) / (2 * np.power(x, 2))

    def integrationT(x):
        return np.power(U(x), 2) * (q * x + np.power(np.e, (-1 * q * x)) - 1)

    T = quad(integrationT, 0, np.inf, limit=limit)[0]

    return 2 * np.pi * μ0 * np.power(nc, 2) * np.power(ri, 5) * T


def calcK(Ra, Rb, d):
    return np.sqrt((4 * Ra * Rb) / (np.power((Ra + Rb), 2) + np.power(d, 2)))


@lru_cache()
def calcM(Ra, Rb, d):
    k = calcK(Ra, Rb, d)

    return μ0 * np.sqrt(Ra * Rb) * ((2 / k - k) * eK(k) - (2 / k) * eE(k))


@lru_cache()
def calcdM(Ra, Rb, d):
    k = calcK(Ra, Rb, d)

    return (μ0 * k * d * (2 * (1 - np.power(k, 2)) * eK(k) - (2 - np.power(k, 2)) * eE(k))) / (4 * (1 - np.power(k, 2)) * np.sqrt(Ra * Rb))


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

        self.nc = self.n / ((self.re - self.ri) * self.l)

        self.R = (self.SR * self.k * np.pi * (np.power(self.re, 2) -
                  np.power(self.ri, 2)) * self.l) / np.power(self.Swire, 2)
        self.L = calcL(self.ri, self.re, self.l, self.nc)


class armature():
    def __init__(self, rai, rae, la, resistivity, v0, ma, m, n, x0, limit = 200):
        self.ri = rai
        self.re = rae
        self.l = la
        self.SR = resistivity
        self.v0 = v0
        self.ma = ma
        self.m = m
        self.n = n
        self.x = x0

        self.__currentFilamentL = self.l / self.m
        self.__currentFilamentNc = 1 / (self.__currentFilamentL * (self.__currentFilamentRe(1) - self.__currentFilamentRi(1)))

        self.R = self.__R()
        self.L = self.__L(limit)

    def __currentFilamentRi(self, j):
        return self.ri + (self.re - self.ri) * (j - 1) / self.n

    def __currentFilamentRe(self, j):
        return self.ri + (self.re - self.ri) * j / self.n

    def currentFilamentAR(self, j):
        '''Average radius'''
        return self.__currentFilamentRi(j) + 0.5 * self.__currentFilamentL

    def currentFilamentX(self, i):
        '''Calculate position of currentFilament'''
        return self.x - 0.5*self.l + (i - 0.5) * self.__currentFilamentL

    def updatePosition(self, delta):
        '''update position of armature'''
        self.x += delta

    def __R(self):
        deltaR = 2 * np.pi * self.SR * self.m / self.l
        R = [2 * np.pi * self.SR * ((self.m / (2 * self.l)) + (self.m * self.n * self.ri / (self.l * (self.re - self.ri))))]

        for k in range(0, self.n - 1):
            R.append(R[k] + deltaR)

        return R * self.m

    def __L(self, limit):
        L = []

        for l in range(1, self.n + 1):
            L.append(calcL(self.__currentFilamentRi(l), self.__currentFilamentRe(l),
                     self.__currentFilamentL, self.__currentFilamentNc, limit))

        return L * self.m
