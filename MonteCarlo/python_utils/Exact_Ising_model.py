import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import scipy.integrate as integrate
import mpmath as mp

"""
Exact Calculation of the 2D Ising model.
"""
def C(T,J=1):
    b = 1/T #beta
    k = 2*np.sinh(2*b*J)/(np.cosh(2*b*J)**2)   #kappa
    K1, _ = integrate.quad(lambda x: 1/np.sqrt(1-k**2*np.sin(x)**2),0,np.pi/2)
    E1, _ = integrate.quad(lambda x: np.sqrt(1-k**2*np.sin(x)**2),0,np.pi/2)
    tanh = np.tanh(2*b*J)
    return 4/np.pi*((b*J/tanh)**2)*(K1-E1-(1-tanh**2)*(np.pi/2+(2*tanh**2-1)*K1)) # boltzman k = 1

def C2(T,J=1):
    b = 1/T
    Tc = 2/np.log(1+np.sqrt(2))
    return -2/np.pi*(2*J*b)**2*np.log(np.abs(1-T/Tc))

def C3(T,J=1):
    b = 1/T #beta
    k = 2*np.sinh(2*b*J)/(np.cosh(2*b*J)**2)   #kappa
    mp.mp.dps = 50
    K1 = mp.quad(lambda x: 1/mp.sqrt(1-k**2*mp.sin(x)**2),[0,mp.pi/2])
    E1 = mp.quad(lambda x: mp.sqrt(1-k**2*mp.sin(x)**2),[0,mp.pi/2])
    tanh = np.tanh(2*b*J)
    return 4/np.pi*((b*J/tanh)**2)*(K1-E1-(1-tanh**2)*(np.pi/2+(2*tanh**2-1)*K1)) # boltzman k = 1

"""
Exact critical temperature of 2D Ising model.
It's about 2.269185314213022
"""
def Exact_critical():
    return 2/np.log(1+np.sqrt(2))

"""
It returns Exact Calculation of the 2D Ising model.
Parameters
Ebin: Binning of exact magnetization calculation
Ebin2: Binning of exact specific heat calculation
"""
def Exact_calc(Ebin = 200, Ebin2 = 200):
    ExactT = np.zeros(Ebin)
    ExactM = np.zeros_like(ExactT)

    Critical_Temp = Exact_critical()

    for i in range(len(ExactM)):
        ExactT[i] = (Critical_Temp-0.01)/Ebin*(i+1)
        ExactM[i] = (1-np.sinh(2*(1/ExactT[i]))**-4)**(1/8)

    ExactT = np.append(ExactT,Critical_Temp)
    ExactM = np.append(ExactM,0)

    Ebin2 = 200
    ExactT2 = np.zeros(Ebin2)
    ExactC = np.zeros_like(ExactT2)

    for i in range(len(ExactC)):
        ExactT2[i] = 5/Ebin2*(i+1)
        ExactC[i] = C(ExactT2[i])
    return [ExactT, ExactM], [ExactT2, ExactC]