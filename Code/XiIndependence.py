from pymule import *
import numpy as np
from plotting import *
from theo_calc import *

homedir = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/"

lo_outs = [
    "17_06_200MeV_Q2big_xi001_final",
    "17_06_200MeV_Q2big_xi01_final"
    ]#,
    #"17_06_200MeV_Q2big_xi0001"]

# Lade F und R aus allen drei Runs
sigF_all = None
sigR_all = None

for run in lo_outs:
    setup(folder=homedir + run + "/out")
    sF = alpha * sigma('mp2mpF')
    sR = alpha * sigma('mp2mpR')
    if sigF_all is None:
        sigF_all = sF
        sigR_all = sR
    else:
        sigF_all = sigF_all + sF
        sigR_all = sigR_all + sR

# Jetzt mergefkswithplot mit den summierten XiRecords
fig, res = mergefkswithplot(
    sigF_all,
    sigR_all
)