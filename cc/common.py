# common.py

import numpy as np
from scipy.integrate import quad
from scipy.special import ellipk, ellipe, j0, j1, struve
from functools import lru_cache

pi = 3.141592653589793
e  = 2.718281828459045
Mu0 = 4 * pi * np.float_power(10, -7)

parameterSets = {
    "r1": 1, "r2": 1, "r3": 1, "r4": 1,
    "ld": 1, "la": 1,
    "RhoD": 0, "RhoA": 1,
    "k": 1, "Sw": 0, "nd": 1,
    "m": 1, "n": 1,
}

p = parameterSets


def Rd():
    V = pi * (p["r2"] * p["r2"] - p["r1"] * p["r1"]) * p["ld"]
    l = p["k"] * V / p["Sw"]
    Rd = p["RhoD"] * l / p["SW"]

    return Rd


@lru_cache(maxsize=1)
def Ra(j):
    if j == 1:
        rj = p["r3"] + j * (p["r4"] - p["r3"]) / p["n"]
        r = rj - (p["r4"] - p["r3"]) / (2 * p["n"])
        l = 2 * pi * r
        Sa = p["la"] * (p["r4"] - p["r3"]) / (p["m"] * p["n"])
        Ra = p["RhoA"] * l / Sa

        return Ra
    else:
        return Ra(j - 1) + p["m"] / p["la"]


def T(q, p):
    def integrand(x):
x = np.linspace(-10, 10, 1000)
    return quad(integrand, 0, np.inf)


def L(nc, rci, rco, lc):
    return 2 * pi * Mu0 * np.power(nc, 2) * np.power(rci, 5) * T(rco / rci, lc / rci)

def Ld():
    return L(p["nc"], p["r1"], p["r2"], p["ld"])


da = (p["r4"] - p["r3"]) / p["n"]
lc = p["la"] / p["m"]


def La(j):
    rci = p["r3"] + j * da

    return L(1, rci, rci + da, lc)


# DEBUG
print("[OK]INITED")
print(T(1.4, 1.86))