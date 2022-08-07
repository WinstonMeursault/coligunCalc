import numpy as np
import numpy.matlib

from basics import *


class singleStageCoilgun():
    def __init__(self, drivingCoil, armature, U, C, deltaT):
        self.drivingCoil = drivingCoil
        self.armature = armature
        self.C = C
        self.dt = deltaT

        self.cache = {"Uc(n-1)" : None, "Id(n-1)" : None}

        self.R = np.diag([self.drivingCoil.R] + self.armature.R)

        self.L = np.diag([self.drivingCoil.L] + self.armature.L)

        self.M = np.zeros((self.armature.m * self.armature.n + 1, self.armature.m * self.armature.n + 1))
        for i in range(1, self.armature.m + 1):
            for j in range(1, self.armature.n + 1):
                for k in range(1, self.armature.m + 1):
                    for l in range(1, self.armature.n + 1):
                        self.M[(i - 1) * self.armature.n + j][(k - 1) * self.armature.n + l] = calcM(self.drivingCoil.r, self.armature.currentFilamentAR(l),
                                                                                                     abs(self.armature.currentFilamentX(k) - self.drivingCoil.x))
        for x in range(0, self.M.shape[0]):
            for y in range(0, self.M.shape[1]):
                if x == y:
                    self.M[x][y] = 0

        self.__updatedM()
        self.__updateM1()

        self.U = np.matrix([U] + [0] * self.armature.m * self.armature.n).T

        self.I = np.zeros((self.armature.m * self.armature.n + 1, 1))
        self.Id = self.U / (self.L - self.M1)

    def __delattr__(self, __name):
        calcM.cache_clear()
        calcdM.cache_clear()

        super().__delattr__(__name)

    def __updatedM(self):
        self.M1 = np.zeros((self.armature.m * self.armature.n + 1, self.armature.m * self.armature.n + 1))
        for i in range(1, self.armature.m + 1):
            for j in range(1, self.armature.n + 1):
                self.M1[0][(i - 1) * self.armature.n + j] = calcM(self.drivingCoil.r, self.armature.currentFilamentAR(j),
                                                                  abs(self.armature.currentFilamentX(i) - self.drivingCoil.x))
        self.M1 = self.M1 + self.M1.T

    def __updateM1(self):
        self.dM = np.zeros((self.armature.m * self.armature.n + 1, self.armature.m * self.armature.n + 1))
        for i in range(1, self.armature.m + 1):
            for j in range(1, self.armature.n + 1):
                self.M1[0][(i - 1) * self.armature.n + j] = calcdM(self.drivingCoil.r, self.armature.currentFilamentAR(j),
                                                                   abs(self.armature.currentFilamentX(i) - self.drivingCoil.x))
        self.dM = self.dM + self.dM.T

    def cache(self):
        self.cache["Uc(n-1)"] = self.U[0][0]
        self.cache["Id(n-1)"] = self.Id[0][0]

    def __update(self):
        self.__updatedM()
        self.__updateM1()

        # Update [U]
        Ucn = self.Register["Uc(n-1)"] - self.dt * self.cache["Id(n-1)"]
        self.U = np.matrix([Ucn] + [0] * self.armature.m * self.armature.n).T

    def __runT(self):
        pass

    def run(self, x = None):
        pass


class multiStageCoilgun():
    def __init__(self):
        pass

# ----------------------------------------------------------------TESTS ----------------------------------------------------------------


if __name__ == "__main__":
    import time
    tic = time.perf_counter()

    print("[INFO]READY...")

    dc = drivingCoil(0.015125, 0.029125, 0.063, 63, 0.0000000175, 0.000028, 0.7)
    print("[OK]dc initialized")

    # a = armature(0.012, 0.015, 0.07, 0.0000000175, 0, 6.384, 70000, 3000, 0.075)
    a = armature(0.012, 0.015, 0.07, 0.0000000175, 0, 6.384, 50, 50, 0.0505)
    print("[OK]a initialized")

    sscg = singleStageCoilgun(dc, a, 8200, 0.001, 0.001)
    print("[OK]sscg initialized")

    # (Ita, v) = sscg.run()
    # print("[RESULT]Ita = " + str(ita * 100) + "%, RIGHT = 8.71%")
    # print("[RESULT]v = " + str(v) + "m/s, RIGHT = 66.07m/s")
    print("[INFO]OVER")

    toc = time.perf_counter()
    t = toc - tic
    print("[DEBUG]TIME: " + str(t) + "s, " + str(t / 60) + "min.")
