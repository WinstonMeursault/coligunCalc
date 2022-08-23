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

        self.cache = {"U(n-1)": None, "I(n-1)": None, "Id(n-1)": None}
        self.__cache()

        self.R = np.matrix(np.diag([self.drivingCoil.R] + self.armature.R))

        self.L = np.matrix(np.diag([self.drivingCoil.L] + self.armature.L))

        self.M = np.matrix(np.zeros((self.armature.m * self.armature.n + 1, self.armature.m * self.armature.n + 1)))
        for i in range(1, self.armature.m + 1):
            for j in range(1, self.armature.n + 1):
                for k in range(1, self.armature.m + 1):
                    for l in range(1, self.armature.n + 1):
                        self.M[(i - 1) * self.armature.n + j, (k - 1) * self.armature.n + l] = calcM(self.drivingCoil.r, self.armature.currentFilamentAR(l),
                                                                                                     abs(self.armature.currentFilamentX(k) - self.drivingCoil.x))
        for x in range(1, self.M.shape[0]):
            for y in range(1, self.M.shape[1]):
                if x == y:
                    self.M[x, y] = 0

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
                self.dM1[0, (i - 1) * self.armature.n + j] = calcdM(self.drivingCoil.r, self.armature.currentFilamentAR(j),
                                                                    abs(self.armature.currentFilamentX(i) - self.drivingCoil.x))
        self.dM1 = self.dM1 + self.dM1.T

    def __cache(self):
        self.cache["U(n-1)"] = self.U
        self.cache["I(n-1)"] = self.I
        self.cache["Id(n-1)"] = self.Id

    def __update(self):
        self.__updateM1()
        self.__updatedM1()

        self.U = np.matrix([self.cache["U(n-1)"][0, 0] - self.dt * self.cache["Id(n-1)"][0, 0]] + [0] * self.armature.m * self.armature.n).T

        self.Id = (self.U + self.Va * self.dM1 * self.I - self.R * self.I - self.M * self.I) / (self.L - self.M1)
        self.I = self.cache["I(n-1)"] + self.dt * self.cache["Id(n-1)"]

        self.F = 0
        for i in range(1, self.armature.m * self.armature.n + 1):
            self.F += self.dM1[0, i] * self.I[0, i]
        self.F = -1 * self.I[0, 0] * self.F

    def run(self, xn):
        while self.armature.x < xn:
            self.__cache()
            self.__update()

        η = (self.armature.ma * (np.power(self.Va, 2) - np.power(self.armature.v0, 2))) / (self.C * np.power(self.U, 2))

        return (η, self.Va)


class multiStageCoilgun():
    def __init__(self):
        pass

# ----------------------------------------------------------------TESTS ----------------------------------------------------------------


def TestA():
    import time
    tic = time.perf_counter()

    print("[INFO]READY...")

    dc = drivingCoil(0.015125, 0.029125, 0.063, 63,
                     0.0000000175, 0.000028, 0.7)
    print("[OK]dc initialized")

    # a = armature(0.012, 0.015, 0.07, 0.0000000175, 0, 6.384, 70000, 3000, 0.075)
    a = armature(0.012, 0.015, 0.07, 0.0000000175, 0, 6.384, 10, 10, 0.0505)
    print("[OK]a initialized")

    sscg = singleStageCoilgun(dc, a, 8200, 0.001, 0.001)
    print("[OK]sscg initialized")

    (η, v) = sscg.run(1.2)
    print("[RESULT]η = " + str(η * 100) + "%, RIGHT = 8.71%")
    print("[RESULT]v = " + str(v) + "m/s, RIGHT = 66.07m/s")
    print("[INFO]OVER")

    toc = time.perf_counter()
    t = toc - tic
    print("[DEBUG]TIME: " + str(t) + "s, " + str(t / 60) + "min.")


if __name__ == "__main__":
    TestA()
