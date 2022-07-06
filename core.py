# singleStage.py

import numpy.matlib 
import numpy as np

from basics import *


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


class singleStageCoilgun():
    def __init__(self, drivingCoil, armature, U, C, deltaT):
        self.drivingCoil = drivingCoil
        self.armature = armature
        self.C = C
        self.deltaT = deltaT

        self.beforeRegister = {"Uc": None, "Id": None}
        self.nowRegister = {"Uc": None, "Id": None}

        # 计算常数矩阵[R]       电阻
        self.armatureR = []
        for i in range(1, self.armature.m + 1):
            for j in range(1, self.n + 1):
                self.armatureR.append(self.armature.currentFilaments[i][j].R)

        self.R = np.diag([self.drivingCoil.R] + self.armatureR)
        del self.armatureR

        # 计算常数矩阵[L]       电感
        self.armatureL = []
        for i in range(1, self.armature.m + 1):
            for j in range(1, self.n + 1):
                self.armatureL.append(self.armature.currentFilaments[i][j].L)

        self.L = np.diag([self.drivingCoil.L] + self.armatureL)
        del self.armatureL

        # 计算常数矩阵[M]       各电阻丝之间的互感
        # TODO
        self.M = np.zeros((self.armature.m * self.armature.n + 1, self.armature.m * self.armature.n + 1))

        # 初始化时变矩阵[M1]    驱动线圈与各电阻丝之间的互感
        self.M1 = np.zeros((self.armature.m * self.armature.n + 1, self.armature.m * self.armature.n + 1))
        for i in range(1, self.armature.m + 1):
            for j in range(1, self.armature.n + 1):
                self.M1[0][(i - 1) * self.armature.n + j] = M(self.drivingCoil.r, self.armature.currentFilaments[i][j].r, 
                                                                abs(self.armature.currentFilaments[i][j].r - self.drivingCoil.x))

        # 初始化时变矩阵[dM1/dx]    驱动线圈与各电阻丝之间的互感梯度
        self.dM = np.zeros((self.armature.m * self.armature.n + 1, self.armature.m * self.armature.n + 1))
        for i in range(1, self.armature.m + 1):
            for j in range(1, self.armature.n + 1):
                self.M1[0][(i - 1) * self.armature.n + j] = dM(self.drivingCoil.r, self.armature.currentFilaments[i][j].r, 
                                                                abs(self.armature.currentFilaments[i][j].r - self.drivingCoil.x))

        # 初始化时变矩阵[U]     电容电压
        self.U = np.matrix([U] + [0] * self.armature.m * self.armature.n).T 

        # 初始化待求矩阵[I] / [dI]     驱动回路电流及其导
        self.I  = np.matlib.zeros((self.armature.m * self.armature.n + 1, 1))
        self.Id = self.U / (self.L - self.M1)
    
    def __delattr__(self, __name):
        M.cache_clear()
        dM.cache_clear()

        super().__delattr__(__name)

    def __update(self):
        # Update [M1]

        # Update [dM1/dx]

        # Update [U]
        Uc = self.beforeRegister["Uc"] - \
            self.deltaT * self.beforeRegister["Id"]

        self.nowRegister["Uc"] = Uc
        self.U = np.matrix([Uc] + [0] * self.armature.m * self.armature.n).T

        self.beforeRegister = self.nowRegister

    def __runT(self):
        pass

    def run(self,x):    
        pass


class multiStageCoilgun():
    def __init__(self):
        pass