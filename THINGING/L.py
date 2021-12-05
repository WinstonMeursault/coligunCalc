from math import pi, e, pow
import numpy as np
from scipy.integrate import quad
from scipy.special import j1

μ0 = 4 * pi * pow(10, -7)


def T(p, q):
    integrationU = lambda x: x * j1(x)
    U = lambda x: quad(integrationU, x, p * x)[0] / pow(x, 3)
    
    integrationT = lambda x: pow(U(x), 2) * (q * x  +  pow(e,(-1 * q * x))  -  1)
    T = quad(integrationT, 0, np.inf)
    
    return T

def L(nc, rci, rco, lc):
    return 2 * pi * μ0 * np.power(nc, 2) * np.power(rci, 5) * T(rco / rci, lc / rci)

def testT():
    result, integrationErr = T(1.4, 1.86)
    calculationErr = abs(result - 0.1149)
    err = integrationErr + calculationErr
    
    print("----------------------------------------------------------------")
    print("[RESULT]testT: " + str(result))
    
    print("----------------------------------------------------------------")
    print("[TEST]integrationErr = " + str(integrationErr))
    print("[TEST]calculationErr = " + str(calculationErr))
    print("[TEST]err = " + str(err))
    
    print("----------------------------------------------------------------")
    
    if err <= 0.00005:
        print("[PASS]Test T")
    else:
        print("[FAIL]Test T")
    
if __name__ == "__main__":
    testT()