# singleStage.py

import numpy as np

from basics import *


class drivingCoil():
    def __init__(self, rdi, rde, ld, n, resistivity, Swire, k, s):
        self.ri = rdi
        self.re = rde
        self.l = ld
        self.n = n                      # 线圈匝数
        self.x = s
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
    def __init__(self, rai, rae, la, resistivity, v0, ma, m, n):
        self.ri = rai
        self.re = rae
        self.l = la
        self.SR = resistivity
        self.v0 = v0
        self.ma = ma
        self.m = m
        self.n = n

        self.x = 0

        self.currentFilaments = {}

        for i in range(1, self.m + 1):
            for j in range(1, self.n + 1):
                self.currentFilaments[i][j] = currentFilament(ri=self.currentFilamentR(j - 1),
                                                              re=self.currentFilamentR(
                                                                  j),
                                                              l=self.l / self.m,
                                                              R=None,
                                                              L=None,
                                                              x0=self.l * (i - 0.5) / self.m - 0.5 * self.l)

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

        for k in self.currentFilaments:
            for l in range(1, self.n + 1):
                k[l].R = R[l]

    def L(self):
        for i in range(1, self.m + 1):
            for j in range(1, self.n + 1):
                self.currentFilaments[i][j].L = L(self.currentFilaments[i][j])


class singleStageCoilgun():
    def __init__(self, drivingCoil, armature, U, V):
        self.drivingCoil = drivingCoil
        self.armature = armature
        self.U = U
        self.V = V


class multiStageCoilgun():
    def __init__(self):
        pass
