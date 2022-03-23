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
        self.nc = nc
        self.rci = rci
        self.rco = rco
        self.l = l
        
        self.Ra = (rci + rco) / 2
        
        self.R = self.R()
        self.L = self.L()

    def calcK(self, others, d):
        return sqrt((4 * self.Ra * others.Ra) / (pow((self.Ra + others.Ra), 2) + pow(d, 2)))

    def updatePosition(self, delta):
        self.x += delta

    def R(self):
        pass

    def L(self):
        p = self.rco / self.rci
        q = self.l / self.rci
        
        # U = lambda x: quad(lambda x: x * j1(x), x, p * x)[0] / pow(x, 3)
        U = lambda x: pi * (-j1(x) * struve(0, x)  +  p * j1(p * x) * struve(0, p * x)  + j0(x) * struve(1, x)  -  p * j0(p * x) * struve(1, p * x)) / (2 * pow(x,2))
        
        integrationT = lambda x: pow(U(x), 2) * (q * x  +  pow(e,(-1 * q * x))  -  1)
        T = quad(integrationT, 0, np.inf, limit=100)
        
        return 2 * pi * μ0 * np.power(self.nc, 2) * np.power(self.rci, 5) * T()

    def M(self, others):
        d = abs(self.x - others.x)
        k = self.calcK(others, d)

        return μ0 * sqrt(self.Ra * others.Ra) * ((2 / k - k) * eK(k) - (2 / k) * eE(k))

    def dM(self, others):
        d = abs(self.x - others.x)
        k = self.calcK(others, d)

        return (μ0 * k * d * (2 * (1 - pow(k, 2)) * eK(k) - (2 - pow(k, 2)) * eE(k))) / (4 * (1 - pow(k, 2)) * sqrt(self.Ra * others.Ra))


class drivingCoil(coil):
    def __init__(self,rdi, rde, ld, n):
        self.rdi = rdi
        self.rde = rde
        self.ld = ld
        self.n = n
        
        Sd = 0.5 * (self.rde - self.rdi) * self.ld
        nc = self.n / Sd
        
        #TODO 将Sd/nc移入Class coil
        
        coil.__init__(self, 0, nc, rdi, rde, ld)
    
    def R(self):
        # TODO

        pass

class armature(coil):
    def __init__(self, x0, nc, rci, rco, l):
        coil.__init__(self, x0, nc, rci, rco, l)
    
    def R(self):
        # TODO

        pass
