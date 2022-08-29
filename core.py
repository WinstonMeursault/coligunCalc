from copy import deepcopy

import numpy as np
import numpy.matlib

from basics import *


class singleStageCoilgun():
    def __init__(self, drivingCoil, armature, U, C, deltaT):
        self.drivingCoil = drivingCoil
        self.armature = armature
        self.C = C
        self.dt = deltaT

        self.t = 0

        self.U = np.matrix([U] + [0] * self.armature.m * self.armature.n).T
        self.I = np.matrix(np.zeros((self.armature.m * self.armature.n + 1, 1)))
        self.Id = np.matrix(np.zeros((self.armature.m * self.armature.n + 1, 1)))

        self.cache = None
        self.__cache()

        self.R = np.matrix(np.diag([self.drivingCoil.R] + self.armature.R))

        self.L = np.matrix(np.diag([self.drivingCoil.L] + self.armature.L))

        self.M = np.matrix(np.zeros((self.armature.m * self.armature.n + 1, self.armature.m * self.armature.n + 1)))
        for i in range(1, self.armature.m + 1):
            for j in range(1, self.armature.n + 1):
                for k in range(1, self.armature.m + 1):
                    for l in range(1, self.armature.n + 1):
                        self.M[(i - 1) * self.armature.n + j, (k - 1) * self.armature.n + l] = \
                            calcM(self.drivingCoil.r, self.armature.currentFilamentAR(l), abs(self.armature.currentFilamentX(k) - self.drivingCoil.x))
        for x in range(1, self.M.shape[0]):
            for y in range(1, self.M.shape[1]):
                if x == y:
                    self.M[x, y] = 0
                    
        self.__updateM1()
        self.__updatedM1()

        self.Id = self.U / (self.L - self.M1)

        self.F = 0
        self.a = 0
        self.Va = 0

    def __delattr__(self, __name):
        calcM.cache_clear()
        calcdM.cache_clear()

        super().__delattr__(__name)

    def __updateM1(self):
        self.M1 = np.matrix(np.zeros((self.armature.m * self.armature.n + 1, self.armature.m * self.armature.n + 1)))
        for i in range(1, self.armature.m + 1):
            for j in range(1, self.armature.n + 1):
                self.M1[0, (i - 1) * self.armature.n + j] = calcM(self.drivingCoil.r, self.armature.currentFilamentAR(j),
                                                                  abs(self.armature.currentFilamentX(i) - self.drivingCoil.x))
        self.M1 = self.M1 + self.M1.T

    def __updatedM1(self):
        self.dM1 = np.matrix(np.zeros((self.armature.m * self.armature.n + 1, self.armature.m * self.armature.n + 1)))
        for i in range(1, self.armature.m + 1):
            for j in range(1, self.armature.n + 1):
                self.dM1[0, (i - 1) * self.armature.n + j] = \
                    calcdM(self.drivingCoil.r, self.armature.currentFilamentAR(j), abs(self.armature.currentFilamentX(i) - self.drivingCoil.x))
        self.dM1 = self.dM1 + self.dM1.T

    def __cache(self):
        self.cache = None               # 防止self.cache缓存上一次的self.chache
        self.cache = deepcopy(self)

    def __update(self):
        self.U = np.matrix([self.cache.U[0, 0] - self.dt * self.cache.Id[0, 0]] + [0] * self.armature.m * self.armature.n).T

        self.Id = np.linalg.inv(self.L - self.M1) * (self.U + self.Va * self.dM1 * self.I - self.R * self.I - self.M * self.I)
        self.I = self.cache.I + self.dt * self.cache.Id 

        self.F = 0
        for i in range(1, self.armature.m + 1):
            for j in range(1, self.armature.n + 1):
                self.F += self.dM1[0, self.armature.n * (i - 1) + j] * self.I[self.armature.n * (i - 1) + j, 0] 
        self.F = -1 * self.I[0, 0] * self.F         # BUG * self.F ???
        
        self.a = self.F / self.armature.ma
        self.Va = self.cache.Va + self.cache.a * self.dt
        self.armature.x = self.cache.armature.x + self.cache.Va * self.dt
        
        # print("[DEBUG] F = " + str(self.F) + "N, Va = " + str(self.Va) + "m/s, x = " + str(self.armature.x) + "m.")

        self.__updateM1()
        self.__updatedM1()

    def run(self, xn):
        while self.armature.x < xn:
            self.__cache()
            self.__update() 

        η = (self.armature.ma * (np.power(self.Va, 2) - np.power(self.armature.v0, 2))) / (self.C * np.power(self.U[0, 0], 2))

        return (η, self.Va)


class multiStageCoilgun():
    def __init__(self):
        pass

# ----------------------------------------------------------------TESTS ----------------------------------------------------------------


def TestA():
    import time
    tic = time.perf_counter()

    print("[TEST]TEST A: READY----------------------------------------------------------------")

    dc = drivingCoil(0.015125, 0.029125, 0.063, 63, 0.0000000175, 0.000028, 0.7)
    print("[OK]dc initialized")

    # a = armature(0.012, 0.015, 0.07, 0.0000000175, 0, 6.384, 70000, 3000, 0.0505)
    a = armature(0.012, 0.015, 0.07, 0.0000000175, 0, 6.384, 70, 3, 0.0505)
    print("[OK]a initialized")

    sscg = singleStageCoilgun(dc, a, 8000, 0.001, 0.001)
    print("[OK]sscg initialized")

    (η, v) = sscg.run(1)
    print("[RESULT]η = " + str(η * 100) + "%, RIGHT = 8.71%")
    print("[RESULT]v = " + str(v) + "m/s, RIGHT = 66.07m/s")
    print("[INFO]OVER")

    toc = time.perf_counter()
    t = toc - tic
    print("[DEBUG]TIME: " + str(t) + "s, " + str(t / 60) + "min.")

def TestB():
    import time
    tic = time.perf_counter()

    print("[TEST]TEST B: READY----------------------------------------------------------------")

    dc = drivingCoil(0.032, 0.05, 0.04, 40, 0.0000000175, 0.000028, 0.7)
    print("[OK]dc initialized")

    a = armature(0.02, 0.03, 10,  0.0000000175, 0, 0.6, 40, 10, 0.02)
    print("[OK]a initialized")

    sscg = singleStageCoilgun(dc, a, 4000, 0.001, 0.001)
    print("[OK]sscg initialized")

    (η, v) = sscg.run(1)
    print("[RESULT]η = " + str(η * 100) + "%, RIGHT = 5.38%")
    print("[RESULT]v = " + str(v) + "m/s, RIGHT = 38m/s")
    print("[INFO]OVER")

    toc = time.perf_counter()
    t = toc - tic
    print("[DEBUG]TIME: " + str(t) + "s, " + str(t / 60) + "min.")
    
def TestC():
    import time
    tic = time.perf_counter()

    print("[TEST]TEST C: READY----------------------------------------------------------------")

    dc = drivingCoil(0.037, 0.05, 0.052, 66, 0.0000000175, 0.000028, 0.7)
    print("[OK]dc initialized")

    a = armature(0.023, 0.036, 0.05, 0.0000000175, 0, 0.6, 50, 13, 0.02)
    print("[OK]a initialized")

    sscg = singleStageCoilgun(dc, a, 4000, 0.001, 0.001)
    print("[OK]sscg initialized")

    (η, v) = sscg.run(1)
    print("[RESULT]η = " + str(η * 100) + "%, RIGHT = 13.6%")
    print("[RESULT]v = " + str(v) + "m/s, RIGHT = 61m/s")
    print("[INFO]OVER")

    toc = time.perf_counter()
    t = toc - tic
    print("[DEBUG]TIME: " + str(t) + "s, " + str(t / 60) + "min.")


if __name__ == "__main__":
    TestA()
    TestB()
    TestC()
