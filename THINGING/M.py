from math import pi, sqrt, pow

import matplotlib.pyplot as plt
from numpy.core.function_base import linspace

from scipy.special import ellipk as K
from scipy.special import ellipe as E

μ0 = 4 * pi * pow(10, -7)


def cK(Ra, Rd, d):
    return sqrt((4 * Ra * Rd) / (pow((Ra + Rd), 2) + pow(d, 2)))


def cM(Ra, Rd, d):
    k = cK(Ra, Rd, d)

    return μ0 * sqrt(Ra * Rd) * ((2 / k - k) * E(k) - (2 / k) * E(k))


def cdM(Ra, Rd, d):
    k = cK(Ra, Rd, d)

    return (μ0 * k * d * (2 * (1 - pow(k, 2)) * K(k) - (2 - pow(k, 2)) * E(k))) / (4 * (1 - pow(k, 2)) * sqrt(Ra * Rd))


if __name__ == "__main__":
    Ra = 0.1
    Rd = 0.12
    d = linspace(-1 , 1, 11000)   

    plt.xlabel("M")
    plt.ylabel("d")
    plt.plot(d, [ cM(Ra, Rd, t) for t in d])
    plt.show()

    plt.xlabel("M")
    plt.ylabel("dM")
    plt.plot(d, [cdM(Ra ,Rd, t) for t in d])
    plt.show()