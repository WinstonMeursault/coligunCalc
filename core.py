# singleStage.py

import numpy as np
import numpy.matlib

from basics import *


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
        self.M = np.zeros((self.armature.m * self.armature.n + 1,self.armature.m * self.armature.n + 1))

        # 初始化时变矩阵[M1]    驱动线圈与各电阻丝之间的互感
        self.M1 = np.zeros((self.armature.m * self.armature.n + 1, self.armature.m * self.armature.n + 1))
        for i in range(1, self.armature.m + 1):
            for j in range(1, self.armature.n + 1):
                self.M1[0][(i - 1) * self.armature.n + j] = M(self.drivingCoil.r, self.armature.currentFilaments[i][j].r,
                                                              abs(self.armature.currentFilaments[i][j].r - self.drivingCoil.x))

        self.M1 = self.M1 + self.M1.T

        # 初始化时变矩阵[dM1/dx]    驱动线圈与各电阻丝之间的互感梯度
        self.dM = np.zeros((self.armature.m * self.armature.n + 1, self.armature.m * self.armature.n + 1))
        for i in range(1, self.armature.m + 1):
            for j in range(1, self.armature.n + 1):
                self.M1[0][(i - 1) * self.armature.n + j] = dM(self.drivingCoil.r, self.armature.currentFilaments[i][j].r,
                                                               abs(self.armature.currentFilaments[i][j].r - self.drivingCoil.x))

        self.dM = self.dM + self.dM.T

        # 初始化时变矩阵[U]     电容电压
        self.U = np.matrix([U] + [0] * self.armature.m * self.armature.n).T

        # 初始化待求矩阵[I] / [dI]     驱动回路电流及其导
        self.I = np.matlib.zeros((self.armature.m * self.armature.n + 1, 1))
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

    def run(self, x):
        pass


class multiStageCoilgun():
    def __init__(self):
        pass
