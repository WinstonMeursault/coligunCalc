# common.py

from math import e, pi, pow, sqrt

import numpy as np

from scipy.integrate import quad
from scipy.special import j0, j1, struve
from scipy.special import ellipe as eE
from scipy.special import ellipk as eK


μ0 = 4 * pi * pow(10, -7)


def L(coilA, limit = 125):
    p = coilA.re / coilA.ri
    q = coilA.l / coilA.ri
        
    # U = lambda x: quad(lambda x: x * j1(x), x, p * x)[0] / pow(x, 3)
    U = lambda x: pi * (-j1(x) * struve(0, x)  +  p * j1(p * x) * struve(0, p * x)  + j0(x) * struve(1, x)  -  p * j0(p * x) * struve(1, p * x)) / (2 * pow(x,2))
        
    integrationT = lambda x: pow(U(x), 2) * (q * x  +  pow(e,(-1 * q * x))  -  1)
    T = quad(integrationT, 0, np.inf, limit = limit)[0]
        
    return 2 * pi * μ0 * np.power(coilA.nc, 2) * np.power(coilA.ri, 5) * T

def calcK(coilA, coilB, d):
    return sqrt((4 * coilA.Ra * coilB.Ra) / (pow((coilA.Ra + coilB.Ra), 2) + pow(d, 2)))

def M(coilA, coilB):
    d = abs(coilA.x - coilB.x)
    k = calcK(coilA, coilB, d)

    return μ0 * sqrt(coilA.Ra * coilB.Ra) * ((2 / k - k) * eK(k) - (2 / k) * eE(k))

def dM(coilA, coilB):
    d = abs(coilA.x - coilB.x)
    k = calcK(coilA, coilB, d)

    return (μ0 * k * d * (2 * (1 - pow(k, 2)) * eK(k) - (2 - pow(k, 2)) * eE(k))) / (4 * (1 - pow(k, 2)) * sqrt(coilA.Ra * coilB.Ra))


class drivingCoil():
    def __init__(self,rdi, rde, ld, n, x0):
        self.ri = rdi
        self.re = rde
        self.l = ld
        self.n = n
        self.x0 = x0
        
        self.Ra = (self.ri + self.re) / 2
        
        Sd = (self.re - self.ri) * self.l
        self.nc = self.n / Sd
        
        self.L = L(self)

    def R(self):
        # TODO

        pass

class armature():
    def __init__(self, x0, nc, rci, rco, l):
        pass
    
    def R(self):
        # TODO

        pass
