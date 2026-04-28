import numpy as np
import matplotlib.pyplot as plt
from plotting import *

def dsigma_dcosth(theta):  # Gives dsigmadcostheta in mubarn!!! LO!!!
    Emu = 100
    mp = 0.938272088
    mmu = 0.105658375
    pmuvec = np.sqrt(Emu**2 - mmu**2)

    a = Emu / mp
    m = mmu**2 / (4 * mp**2)

    xi = (2 - 2*np.cos(theta)*np.sqrt(1 - m*np.sin(theta)**2)) / (1 + 4*m*np.cos(theta)**2)
    x = 2 * xi / (2 + a * xi)

    Q2 = x * pmuvec**2
    alpha = 1/137.035999084

    tau = Q2 / (4 * mp**2)
    s = 2*Emu*mp + mp**2 + mmu**2

    R = (pmuvec**2 - tau*(s - 2*mp**2*(1+tau))) / (pmuvec**2*(1+tau))
    eps = (Emu**2 - tau*(s - mmu**2)) / (pmuvec**2 - tau*(s - 2*mp**2*(1+tau)))

    mup = 2.79284734
    GE = (1 + Q2/0.71)**(-2)
    GM = mup * GE


    #print("GE: ",GE)
    #print("GM: ",GM)
    #print("Q2: ",Q2)

    #ddQ2 = (2+a*xi)**2/(4*pmuvec**2) * (1-2*m*xi*np.cos(theta))/(2*(1-m*xi**2)) * ddcostheta
    #ddcostheta = 1 / ((2+a*xi)**2/(4*pmuvec**2) * (1-2*m*xi*np.cos(theta))/(2*(1-m*xi**2))) ddQ2
    dsigmadQ2 = (4*np.pi*alpha**2 / Q2**2) * R * (eps*GE**2 + tau*GM**2)

    jac = (2+a*xi)**2/(4*pmuvec**2) * (1-2*m*xi*np.cos(theta))/(2*(1-m*xi**2))

    return dsigmadQ2 / jac * (0.389379365 * 1e3)