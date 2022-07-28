from functools import lru_cache

from numpy import inf
from scipy.integrate import quad
from scipy.special import ellipe as eE
from scipy.special import ellipk as eK
from scipy.special import j0, j1, struve

cdef long double pi = 3.141592653589793
cdef long double e = 2.71828182845678
cdef long double Mu0 = 4 * pi * pow(10, -7)


# U = lambda x: quad(lambda x: x * j1(x), x, p * x)[0] / pow(x, 3)
cdef U(long double x, double p):
    return pi * (-j1(x) * struve(0, x) + p * j1(p * x) * struve(0, p * x) + j0(x) * struve(1, x) - p * j0(p * x) * struve(1, p * x)) / (2 * pow(x, 2))

cdef integrationT(long double x, double p, double q):
    return pow(U(x, p), 2) * (q * x + pow(e, (-1 * q * x)) - 1)

cdef calcL(double ri, double re, double l, double nc,int limit=200):
    cdef double p = re / ri
    cdef double q = l / ri

    return 2 * pi * Mu0 * pow(nc, 2) * pow(ri, 5) * quad(integrationT, 0, inf, limit=limit, args = (p, q))[0]


cdef calcK(double Ra, double Rb, double d):
    return pow((4 * Ra * Rb) / (pow((Ra + Rb), 2) + pow(d, 2)), 0.5)


#@lru_cache()
cdef calcM(double Ra, double Rb, double d):
    cdef double k = calcK(Ra, Rb, d)

    return Mu0 * pow(Ra * Rb, 0.5) * ((2 / k - k) * eK(k) - (2 / k) * eE(k), 0.5)


#@lru_cache()
cdef calcdM(double Ra, double Rb, double d):
    cdef double k = calcK(Ra, Rb, d)

    return (Mu0 * k * d * (2 * (1 - pow(k, 2)) * eK(k) - (2 - pow(k, 2)) * eE(k))) / (4 * (1 - pow(k, 2)) * pow(Ra * Rb, 0.5))


class drivingCoil():
    def __init__(self, double rdi, double rde, double ld, int n, double resistivity, double Swire, float k):
        self.ri = rdi
        self.re = rde
        self.l = ld
        self.n = n                         # 线圈匝数
        self.x = 0.5 * self.l
        self.SR = resistivity           # 电阻率
        self.Swire = Swire              # 单根导线的截面积
        self.k = k                      # 驱动线圈填充率

        self.r = (self.ri + self.re) / 2

        self.nc = self.n / ((self.re - self.ri) * self.l)

        self.R = (self.SR * self.k * pi * (pow(self.re, 2) - pow(self.ri, 2)) * self.l) / pow(self.Swire, 2)
        self.L = calcL(self.ri, self.re, self.l, self.nc)


class armature():
    def __init__(self, double rai, double rae, double la, double resistivity, double v0, double ma, int m, int n, double x0, unsigned short int	 limit = 200):
        self.ri = rai
        self.re = rae
        self.l = la
        self.SR = resistivity
        self.v0 = v0
        self.ma = ma
        self.m = m
        self.n = n
        self.x = x0

        self.__currentFilamentL = self.l / self.m
        self.__currentFilamentNc = 1 / (self.__currentFilamentL * (self.__currentFilamentRe(1) - self.__currentFilamentRi(1)))

        self.R = self.__R()
        self.L = self.__L(limit)

    def __currentFilamentRi(self, int j):
        return self.ri + (self.re - self.ri) * (j - 1) / self.n

    def __currentFilamentRe(self, int j):
        return self.ri + (self.re - self.ri) * j / self.n

    def currentFilamentAR(self, int j):
        '''Average radius'''
        return self.__currentFilamentRi(j) + 0.5 * self.__currentFilamentL

    def currentFilamentX(self, int i):
        '''Calculate position of currentFilament'''
        return self.x - 0.5*self.l + (i - 0.5) * self.__currentFilamentL

    def updatePosition(self, double delta):
        '''update position of armature'''
        self.x += delta

    def __R(self):
        cdef long double deltaR = 2 * pi * self.SR * self.m / self.l
        cdef list R = [2 * pi * self.SR * ((self.m / (2 * self.l)) + (self.m * self.n * self.ri / (self.l * (self.re - self.ri))))]

        for k in range(0, self.n - 1):
            R.append(R[k] + deltaR)

        return R * self.m

    def __L(self, unsigned short int limit):
        cdef list L = []

        for l in range(1, self.n + 1):
            L.append(calcL(self.__currentFilamentRi(l), self.__currentFilamentRe(l), self.__currentFilamentL, self.__currentFilamentNc, limit))

        return L * self.m
