import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from pymule import *
import numpy as np
from plotting import *
from theo_calc import *

plt.rcParams['text.usetex'] = True

homedir = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/"

setup(folder=homedir + "combined_xi_out/out")
sigF = alpha * sigma('mp2mpF')
sigR = alpha * sigma('mp2mpR')
fig, res = mergefkswithplot(sigF, sigR)
fig.savefig("Figures/xi_independence.pdf", bbox_inches='tight')