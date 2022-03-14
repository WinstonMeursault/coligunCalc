# common.py

from math import e, pi, pow, sqrt

import numpy as np

from scipy.integrate import quad
from scipy.special import j0, j1, struve
from scipy.special import ellipe as eE
from scipy.special import ellipk as eK


μ0 = 4 * pi * pow(10, -7)


class coil():
    def __init__(self, x0, nc, rci, rco, l):
        self.x = x0
        self.nc = nc    #BUG
        self.rci = rci
        self.rco = rco
        self.l = l
        
        self.Ra = (rci + rco) / 2
        self.p = self.rco / self.rci
        self.q = self.l / self.rci
        
        # U = lambda x: quad(lambda x: x * j1(x), x, p * x)[0] / pow(x, 3)
        U = lambda x: pi * (-j1(x) * struve(0, x)  +  self.p * j1(self.p * x) * struve(0, self.p * x)  + j0(x) * struve(1, x)  -  self.p * j0(self.p * x) * struve(1, self.p * x)) / (2 * pow(x,2))
        self.integrationT = lambda x: pow(U(x), 2) * (self.q * x  +  pow(e,(-1 * self.q * x))  -  1)
        
        self.R = self.R()
        self.L = self.L()

    def T(self):
        return quad(self.integrationT, 0, np.inf, limit=100)

    def calcK(self, others, d):
        return sqrt((4 * self.Ra * others.Ra) / (pow((self.Ra + others.Ra), 2) + pow(d, 2)))

    def updatePosition(self, delta):
        self.x += delta

    def R(self):
        pass

    def L(self):
        return 2 * pi * μ0 * np.power(self.nc, 2) * np.power(self.rci, 5) * self.T()

    def M(self, others):
        d = abs(self.x - others.x)
        k = self.calcK(others, d)

        return μ0 * sqrt(self.Ra * others.Ra) * ((2 / k - k) * eK(k) - (2 / k) * eE(k))

    def dM(self, others):
        d = abs(self.x - others.x)
        k = self.calcK(others, d)

        return (μ0 * k * d * (2 * (1 - pow(k, 2)) * eK(k) - (2 - pow(k, 2)) * eE(k))) / (4 * (1 - pow(k, 2)) * sqrt(self.Ra * others.Ra))


class drivingCoil(coil):
    def __init__(self, x0, nc, rci, rco, l):
        coil.__init__(self, x0, nc, rci, rco, l)
    
    def R(self):
        # TODO

        pass


class armature(coil):
    def __init__(self, x0, nc, rci, rco, l):
        coil.__init__(self, x0, nc, rci, rco, l)
    
    def R(self):
        # TODO

        pass
