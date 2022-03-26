# common.py

from math import e, pi, pow, sqrt

import numpy as np

from scipy.integrate import quad
from scipy.special import j0, j1, struve
from scipy.special import ellipe as eE
from scipy.special import ellipk as eK


μ0 = 4 * pi * pow(10, -7)


def L(coilA, limit=125):
    p = coilA.re / coilA.ri
    q = coilA.l / coilA.ri

    # U = lambda x: quad(lambda x: x * j1(x), x, p * x)[0] / pow(x, 3)
    def U(x): return pi * (-j1(x) * struve(0, x) + p * j1(p * x) * struve(0, p * x) +
                           j0(x) * struve(1, x) - p * j0(p * x) * struve(1, p * x)) / (2 * pow(x, 2))

    def integrationT(x): return pow(U(x), 2) * \
        (q * x + pow(e, (-1 * q * x)) - 1)
    T = quad(integrationT, 0, np.inf, limit=limit)[0]

    return 2 * pi * μ0 * np.power(coilA.nc, 2) * np.power(coilA.ri, 5) * T


def calcK(coilA, coilB, d):
    return sqrt((4 * coilA.r * coilB.r) / (pow((coilA.r + coilB.r), 2) + pow(d, 2)))


def M(coilA, coilB, d=None):
    if d == None:
        d = abs(coilA.x - coilB.x)
    k = calcK(coilA, coilB, d)

    return μ0 * sqrt(coilA.r * coilB.r) * ((2 / k - k) * eK(k) - (2 / k) * eE(k))


def dM(coilA, coilB, d=None):
    if d == None:
        d = abs(coilA.x - coilB.x)
    k = calcK(coilA, coilB, d)

    return (μ0 * k * d * (2 * (1 - pow(k, 2)) * eK(k) - (2 - pow(k, 2)) * eE(k))) / (4 * (1 - pow(k, 2)) * sqrt(coilA.r * coilB.r))


class drivingCoil():
    def __init__(self, rdi, rde, ld, n, x0, resistivity, s, k):
        self.ri = rdi
        self.re = rde
        self.l = ld
        self.n = n                      # 线圈匝数
        self.x = x0
        self.SR = resistivity           # 电阻率
        self.s = s                      # 单根导线的截面积
        self.k = k                      # 驱动线圈填充率

        self.r = (self.ri + self.re) / 2
        Sd = (self.re - self.ri) * self.l

        self.nc = self.n / Sd

        self.R = self.R()
        self.L = L(self)

    def R(self):
        return (self.SR * self.k * pi * (pow(self.re, 2) - pow(self.ri, 2)) * self.l) / pow(self.s, 2)


class currentFilament():
    def __init__(self, ri, re, l, R, L):
        self.ri = ri
        self.re = re
        self.l = l
        self.R = R
        self.L = L

        self.r = (self.ri + self.re) / 2
        
        Sc = l * (self.re - self.ri)
        self.nc = 1 / Sc


class armature():
    def __init__(self, rai, rae, la, x0, resistivity, m, n):
        self.ri = rai
        self.re = rae
        self.l = la
        self.x = x0
        self.SR = resistivity
        self.m = m
        self.n = n

        self.currentFilaments = {}

        for i in range(1, self.m + 1):
            for j in range(1, self.n + 1):
                self.currentFilaments[i][j] = currentFilament(ri=self.currentFilamentR(j - 1), 
                                                              re=self.currentFilamentR(j), 
                                                              l=self.l / self.m, 
                                                              R=None, 
                                                              L=None)

        self.R()
        self.L()

    def currentFilamentR(self, j):
        return self.ri + (self.re - self.ri) * j / self.n

    def updatePosition(self, delta):
        self.x += delta

    def R(self):
        deltaR = 2 * pi * self.SR * self.m / self.l
        R = {
            1: 2 * pi * self.SR * ((self.m / (2 * self.l)) + (self.m * self.n * self.ri / (self.l * (self.re - self.ri))))
        }

        for a in range(2, self.n + 1):
            R[a] = R[a - 1] + deltaR

        for k in self.currentFilaments:
            for l in range(1, self.n + 1):
                k[l].R = R[l]

    def L(self):
        for i in range(1, self.m + 1):
            for j in range(1, self.n + 1):
                self.currentFilaments[i][j].L = L(self.currentFilaments[i][j])
