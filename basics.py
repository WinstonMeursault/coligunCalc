from functools import lru_cache
from typing import List

import numpy as np
from scipy.integrate import quad, nquad
from scipy.special import ellipe as eE
from scipy.special import ellipk as eK
from scipy.special import j0, j1, struve

μ0 = 0.0000012566370614359173


def calcL(ri: float, re: float, l: float, nc: int, limit: int = 200) -> float:
    """计算线圈自感

    Args:
        ri (float): 线圈内径
        re (float): 线圈外径
        l (float): 线圈长度
        nc (int): 线圈匝数
        limit (int, optional): Numpy quad参数, 若警告区间划分数过少可适当提高. Defaults to 200.

    Returns:
        float: 线圈自感
    """    
    
    p = re / ri
    q = l / ri

    # U = lambda x: quad(lambda x: x * j1(x), x, p * x)[0] / np.power(x, 3)
    def U(x: float) -> float:
        return np.pi * (-j1(x) * struve(0, x) + p * j1(p * x) * struve(0, p * x) + j0(x) * struve(1, x) - p * j0(p * x) * struve(1, p * x)) / (2 * np.power(x, 2))  

    def integrationT(x: float) -> float:
        return np.power(U(x), 2) * (q * x + np.power(np.e, (-1 * q * x)) - 1)

    T = quad(integrationT, 0, np.inf, limit=limit)[0]

    return 2 * np.pi * μ0 * np.power(nc, 2) * np.power(ri, 5) * T


def calcK(Ra: float, Rb: float, d: float) -> float:
    return np.sqrt((4 * Ra * Rb) / (np.power(Ra + Rb, 2) + np.power(d, 2)))


@lru_cache()
def calcM(Rai: float, Rae: float, La: float, Na: int, Rbi: float, Rbe: float, Lb: float, Nb: int, d: float) -> float:
    """计算两线圈互感(设为A B线圈)

    Args:
        Rai (float): A线圈内径
        Rae (float): A线圈外径
        La (float): A线圈长度
        Na (int): A线圈匝数
        Rbi (float): B线圈内径
        Rbe (float): B线圈外径
        Lb (float): B线圈长度
        Nb (int): B线圈匝数
        d (float): A B线圈中心面的距离

    Returns:
        float: A B两线圈之间的互感(返回在距离为d时的数值)
    """      

    def integrationM(r1: float, z1: float, r2: float, z2: float) -> float:
        ra = Rai + r1
        rb = Rbi + r2
        dz = abs(abs(d - 0.5 * La + 0.5 * Lb) + z1 - z2)

        k = calcK(ra, rb, dz)
        return μ0 * np.sqrt(ra * rb) * ((2 / k - k) * eK(k) - (2 / k) * eE(k))

    return (Na * Nb * nquad(integrationM, [[0, Lb], [Rbi, Rbe], [0, La], [Rai, Rbe]])[0]) / ((Rae - Rai) * La * (Rbe - Rbi) * Lb)


@lru_cache()
def calcdM(Rai: float, Rae: float, La: float, Na: int, Rbi: float, Rbe: float, Lb: float, Nb: int, d: float) -> float:  
    """计算两线圈互感梯度(设为A B线圈)

    Args:
        Rai (float): A线圈内径
        Rae (float): A线圈外径
        La (float): A线圈长度
        Na (int): A线圈匝数
        Rbi (float): B线圈内径
        Rbe (float): B线圈外径
        Lb (float): B线圈长度
        Nb (int): B线圈匝数
        d (float): A B线圈中心面的距离

    Returns:
        float: A B两线圈之间的互感梯度(返回在距离为d时的数值)
    """      

    def integrationdM(r1: float, z1: float, r2: float, z2: float) -> float:
        ra = Rai + r1
        rb = Rbi + r2
        dz = abs(abs(d - 0.5 * La + 0.5 * Lb) + z1 - z2)

        k = calcK(ra, rb, dz)
        return (μ0 * k * d * (2 * (1 - np.power(k, 2)) * eK(k) - (2 - np.power(k, 2)) * eE(k))) / (4 * (1 - np.power(k, 2)) * np.sqrt(ra * rb))

    return (Na * Nb * nquad(integrationdM, [[0, Lb], [Rbi, Rbe], [0, La], [Rai, Rbe]])[0]) / ((Rae - Rai) * La * (Rbe - Rbi) * Lb)


class drivingCoil():
    def __init__(self, rdi: float, rde: float, ld: float, n: int, resistivity: float, Swire: float, k: float) -> None:
        self.ri = rdi
        self.re = rde
        self.l = ld
        self.n = n                      # 线圈匝数
        self.x = 0.5 * self.l
        self.SR = resistivity           # 电阻率
        self.Swire = Swire              # 单根导线的截面积
        self.k = k                      # 驱动线圈填充率

        self.r = (self.ri + self.re) / 2

        self.nc = self.n / ((self.re - self.ri) * self.l)

        self.R = (self.SR * self.k * np.pi * (np.power(self.re, 2) - np.power(self.ri, 2)) * self.l) / np.power(self.Swire, 2)
        self.L = calcL(self.ri, self.re, self.l, self.nc)


class armature():
    def __init__(self, rai: float, rae: float, la: float, resistivity: float, v0: float, ma: float, m: int, n: int, x0: float, limit: int = 200) -> None:
        self.ri = rai
        self.re = rae
        self.l = la
        self.SR = resistivity
        self.v0 = v0
        self.ma = ma
        self.m = m
        self.n = n
        self.x = x0

        self.__deltaRN = (self.re - self.ri) / self.n
        self.currentFilamentL = self.l / self.m
        self.__currentFilamentNc = 1 / (self.currentFilamentL * self.__deltaRN)

        self.R = self.__R()
        self.L = self.__L(limit)

    def currentFilamentRi(self, j: int) -> float: 
        """计算电流丝内径

        Args:
            j (int): 电流丝径向编号

        Returns:
            float: 电流丝内径
        """             
        return self.ri + (self.re - self.ri) * (j - 1) / self.n

    def currentFilamentRe(self, j: int) -> float:
        """计算电流丝外径

        Args:
            j (int): 电流丝径向编号

        Returns:
            float: 电流丝外径
        """        
        return self.ri + (self.re - self.ri) * j / self.n

    def currentFilamentAR(self, j: int) -> float:
        """计算电流丝的平均半径

        Args:
            j (int): 电流丝径向编号

        Returns:
            float: 电流丝平均半径
        """        
        return self.ri + self.__deltaRN * (j - 0.5)

    def currentFilamentX(self, i: int) -> float:
        """计算电流丝绝对位置

        Args:
            i (int): 电流丝轴向编号

        Returns:
            float: 电流丝的绝对位置
        """        
        return self.x - 0.5 * self.l + (i - 0.5) * self.currentFilamentL

    def updatePosition(self, delta: float) -> None:
        """更新电枢位置(将位置更新为原位置与delta之和)

        Args:
            delta (float): 电枢位移
        """        
        self.x += delta

    def __R(self) -> List[float]:
        deltaR = 2 * np.pi * self.SR * self.m / self.l
        R = [2 * np.pi * self.SR * ((self.m / (2 * self.l)) + (self.m * self.n * self.ri / (self.l * (self.re - self.ri))))]

        for k in range(1, self.n):
            R.append(R[k - 1] + deltaR)

        return R * self.m

    def __L(self, limit: int) -> List[float]:
        L = []

        for l in range(1, self.n + 1):
            L.append(calcL(self.currentFilamentRi(l), self.currentFilamentRe(l), self.currentFilamentL, self.__currentFilamentNc, limit))

        return L * self.m
