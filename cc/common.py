# common.py

from math import pi, sqrt, pow

from scipy.special import ellipk as eK
from scipy.special import ellipe as eE

μ0 = 4 * pi * pow(10, -7)


class coil():
    def __init__(self, x0, Ra, Rd):
        self.Ra = Ra
        self.Rd = Rd
        self.x = x0

    def calcK(self, d):
        return sqrt((4 * self.Ra * self.Rd) / (pow((self.Ra + self.Rd), 2) + pow(d, 2)))

    def updatePosition(self, delta):
        self.x += delta
    
    def R(self):
        pass

    def L(self):
        # TODO
        
        pass

    def M(self, others):
        d = abs(self.x - others.x)
        k = self.calcK(d)

        return μ0 * sqrt(self.Ra * self.Rd) * ((2 / k - k) * eK(k) - (2 / k) * eE(k))

    def dM(self, others):
        d = abs(self.x - others.x)
        k = self.calcK(d)

        return (μ0 * k * d * (2 * (1 - pow(k, 2)) * eK(k) - (2 - pow(k, 2)) * eE(k))) / (4 * (1 - pow(k, 2)) * sqrt(self.Ra * self.Rd))

class drivingCoil(coil):
    def R(self):
        # TODO
        
        pass

class armature(coil):
    def R(self):
        # TODO
        
        pass