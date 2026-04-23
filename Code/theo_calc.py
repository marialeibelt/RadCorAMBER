import numpy as np
import matplotlib.pyplot as plt
from plotting import *

theta = np.linspace(1.35e-3, 1.65e-3, 500)  # sinnvoller Bereich
costheta = np.cos(theta)
Emu = 100       #GeV
mp = 0.93827    #GeV
mmu = 0.10566   #GeV
pmu = np.sqrt(Emu**2-mmu**2)
print("pmu: ",pmu)

a = Emu / mp
m = mmu**2 / (4*mp**2)

xi = (2-2*np.cos(theta)*np.sqrt(1-m*np.sin(theta)**2))/(1+4*m*np.cos(theta)**2)
x = 2 * xi / (2 + a * xi)
print("xi: ",xi)
print("x: ",x)

pmuvec = np.sqrt(Emu**2-mmu**2)
Q2 = x * pmuvec**2
print(Q2)

Ekin = Q2/(2*mp)
Emuprime = Emu-Ekin
alpha=1/137.036
tau = Q2/(4*mp**2)
s = 2*Emu*mp+mp**2+mmu**2
R = (pmuvec**2-tau*(s-2*mp**2*(1+tau))) / (pmuvec**2*(1+tau))
eps = (Emu**2-tau*(s-mmu**2))/(pmuvec**2-tau*(s-2*mp**2*(1+tau)))

mup=2.793
GE = (1 + Q2/0.71)**(-2)
GM = mup*GE

dsigmadQ2 = (4*np.pi*alpha**2/Q2**2) * R * (eps*GE**2 + tau*GM**2)
#ddQ2 = (2+a*xi)**2/(4*pmu**2) * (1-2*m*xi*np.cos(theta))/(2*(1-m*xi**2)) * ddcostheta
#ddcostheta = 1 / ((2+a*xi)**2/(4*pmu**2) * (1-2*m*xi*np.cos(theta))/(2*(1-m*xi**2))) ddQ2
dsigmadcostheta = dsigmadQ2 * 1 / ((2+a*xi)**2/(4*pmu**2) * (1-2*m*xi*np.cos(theta))/(2*(1-m*xi**2)))
#print("dsigmadcostheta: ", dsigmadcostheta)
conv = 0.389379 * 1e3 
dsigmadcostheta_mub = dsigmadcostheta * conv

fig, axes = create_figure(figsize=(6,5))
ax = axes[0,0]

# sortieren (wichtig für saubere Linie)
#idx = np.argsort(costheta)
#ax.plot(costheta[idx], dsigmadcostheta_mub[idx])
idx = np.argsort(theta)
ax.plot(theta[idx], dsigmadcostheta_mub[idx])
ax.set_ylim(0.)
style_axis(
    ax,
    xlabel=r"$\cos\theta$",
    ylabel=r"$\frac{d\sigma}{d\cos\theta}\ (\mu b\mathrm{barn})$",
    title=r"$Muon Scattering Angle (lab)",
)

save_figure(fig, "dsigma_dcostheta", outdir="plots")