# singleStage.py

import pstats
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

        self.armatureR = []
        for i in range(1, self.armature.m + 1):
            for j in range(1, self.n + 1):  
                self.armatureR.append(self.armature.currentFilaments[i][j].R)

        self.R = np.diag([self.drivingCoil.R] + self.armatureR)
        del self.armatureR

        self.armatureL = []
        for i in range(1, self.armature.m + 1):
            for j in range(1, self.n + 1):
                self.armatureL.append(self.armature.currentFilaments[i][j].L)

        self.L = np.diag([self.drivingCoil.L] + self.armatureL)
        del self.armatureL

        self.M = np.zeros((self.armature.m * self.armature.n + 1,self.armature.m * self.armature.n + 1))
        for i in range(1, self.armature.m + 1):
            for j in range(1, self.armature.n + 1):
                for k in range(1, self.armature.m + 1):
                    for l in range(1 ,self.armature.n + 1):
                        self.M[(i - 1) * self.armature.n + j][(k - 1) * self.armature.n + l] = M(self.drivingCoil.r, 
                                        self.armature.currentFilaments[k][L], 
                                        abs(self.armature.currentFilaments[k][l].x - self.drivingCoil.x))
        for x in range(1, self.armature.m + 1):
            for y in range(1, self.armature.n + 1):
                if x == y:
                    self.M[x][y] = 0

        self.__updatedM()
        
        self.__updatedM1()

        self.U = np.matrix([U] + [0] * self.armature.m * self.armature.n).T

        self.I = np.matlib.zeros((self.armature.m * self.armature.n + 1, 1))
        self.Id = self.U / (self.L - self.M1)

    def __delattr__(self, __name):
        M.cache_clear()
        dM.cache_clear()

        super().__delattr__(__name)

    def __updatedM(self):
        self.M1 = np.zeros((self.armature.m * self.armature.n + 1, self.armature.m * self.armature.n + 1))
        for i in range(1, self.armature.m + 1):
            for j in range(1, self.armature.n + 1):
                self.M1[0][(i - 1) * self.armature.n + j] = M(self.drivingCoil.r, self.armature.currentFilaments[i][j].r,
                                                              abs(self.armature.currentFilaments[i][j].x - self.drivingCoil.x))
        self.M1 = self.M1 + self.M1.T

    def __updatedM1(self):
        self.dM = np.zeros((self.armature.m * self.armature.n + 1, self.armature.m * self.armature.n + 1))
        for i in range(1, self.armature.m + 1):
            for j in range(1, self.armature.n + 1):
                self.M1[0][(i - 1) * self.armature.n + j] = dM(self.drivingCoil.r, self.armature.currentFilaments[i][j].r,
                                                               abs(self.armature.currentFilaments[i][j].x - self.drivingCoil.x))
        self.dM = self.dM + self.dM.T

    def __update(self):
        self.__updatedM()

        self.__updatedM()

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
