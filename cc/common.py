# common.py

import numpy as np
from scipy.integrate import quad
from scipy.special import ellipk, ellipe, j0, j1, struve
from functools import lru_cache

pi = np.pi
e  = np.e
Mu0 = 4 * pi * np.float_power(10, -7)

