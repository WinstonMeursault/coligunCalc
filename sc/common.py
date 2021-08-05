# common.py

import scipy
from functools import lru_cache

π = 3.141592653589793
μ0 = 4 * π * pow(10, -7)


def Rd(ρd, k, r1, r2, lg, S):
    '''
    驱动线圈电阻                     
    ρd      驱动线圈电阻率             
    k       驱动线圈填充系数(驱动线圈截面绕制导线的总截面积与驱动线圈截面积之比)    
    r1      驱动线圈内径 
    r2      驱动线圈外径                
    lg      驱动线圈轴线长度             
    S       一根绕制驱动线圈的导线的截面积 
    RETURN:计算结果
    '''

    return (ρd * k * π * (r2 * r2 - r1 * r1) * lg) / S


def Ra(m, n, i, j, ρa, r, la):
    '''
    电枢电流丝电阻                   
    m,n     电流丝分割参数            
    i,j     需计算电阻的电流丝的位置参数
    ρa      电枢电阻率                
    r       半径参数组                   
    la      电枢轴线长度               
    '''

    return 2 * π * ρa * ((m * j) / la - 0.5 * m * la + (m * n * r["r3"]) / la * (r["r4"] - r["r3"]))


def x(F, mp, v0=0, x0):
    a = F/mp
    v = v0 + scipy.integrate.quad(a, 0, t)[0]
    x = x0 + scipy.integrate.quad(v, 0, t)[0]
    return x
