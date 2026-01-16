from pymule import *
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt

lo_outs = ["mcmule-release",
           "mp2mp_NLO_22_12",
           "mp2mp_NLO_12_01",
           "mp2mp_NLO_13_01"]

nlo_outs = ["mp2mp_testNLO",
            "mp2mp_NLO_22_12",
           "mp2mp_NLO_12_01",
           "mp2mp_NLO_13_01"]

savenames = ["combined"]

#Specify paths
lo_i = 3
nlo_i = 3
savename = savenames[0]+"_"+nlo_outs[nlo_i]




# Setup
setup(folder="/home/marialei/"+lo_outs[lo_i]+"/out") #LO
lo = mergefks(sigma('mp2mp0'))*alpha**2*conv 
setup(folder="/home/marialei/"+nlo_outs[nlo_i]+"/out") #NLO 
nlo = mergefks(sigma('mp2mpR'),sigma('mp2mpF'),anyxi=sigma('mp2mpA'))*alpha**3*conv

full = lo + nlo

#print(lo.value)
#print(lo["Emu"])
#print(nlo.value)
#print(nlo["Emu"])





# Plotting
fig, axes = plt.subplots(1, 2, figsize=(16, 6))  # 2 Subplots untereinander

# Observable 1: th3
obs1 = 'th3'
obs1_name = r'Scattering Angle'
plt.sca(axes[0]) 
th3_lo_func = errorband(lo[obs1])
th3_lo = lo[obs1]
th3_nlo_func = errorband(nlo[obs1])
th3_nlo = nlo[obs1]
print("th3_nlo: ",th3_nlo)

th3_added = addplots(th3_lo, th3_nlo)
th3_added_func = errorband(th3_added)
#print("\n th3_added: ",th3_added)

#th3_devided = dividenumbers(th3_nlo,th3_lo) #gives back error like divideplots
#print("\n th3_devided: ",th3_devided)


thK = mergebins(divideplots(th3_nlo, th3_added), 5) #K factor nlo/(nlo+lo)  #common K definition
#thK_func = errorband(thK)
print("\n thK: ",thK)

#thK2 = mergebins(divideplots(th3_nlo,th3_lo, offset=+0), 5) #nlo/lo +offset = K-1 +offset #different K definition!!
#thK2_func = errorband(thK2)
#print("\n thK2: ",thK2)


axes[0].grid(True, alpha=0.7, ls='dotted')
axes[0].set_xlabel(r'$\theta_3 (rad)$', fontsize=12)
axes[0].set_ylabel('Counts',fontsize=12)
axes[0].set_title(obs1_name)
axes[0].set_yscale('log')
#axes[0].legend(["ThetaK", "ThetaK offset"], framealpha=0, fontsize=12)
axes[0].legend(['LO','NLO',"LO with NLO correction"], framealpha=0, fontsize=12)
axes[0].set_xlim(1.3e-3,1.7e-3)
axes[0].xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
axes[0].ticklabel_format(style='sci', axis='x', scilimits=(-3,3))


# Observable 2: Emu
obs2 = 'Emu'
obs2_name = 'Energy of the Scattered Muon'
plt.sca(axes[1])

Emu_lo = lo[obs2]
Emu_lo_func = errorband(Emu_lo)
Emu_nlo = nlo[obs2]
Emu_nlo_func = errorband(Emu_nlo)

Emu_full = full[obs2]
#Emu_added = Emu_lo + Emu_nlo #!that does not work like that ->use addplots
Emu_full_func = errorband(Emu_full)

Emu_devided = dividenumbers(Emu_nlo,Emu_lo)
print("\n Emu_devided: ",Emu_devided)

axes[1].grid(True, alpha=0.7, ls='dotted')
axes[1].set_xlabel(r'$E_{\mu} (MeV)$', fontsize=12)
axes[1].set_ylabel('Counts',fontsize=12)
axes[1].set_title(obs2_name)
axes[1].set_yscale('log')
axes[1].legend(['LO','NLO correction', "LO corrected with NLO"], framealpha=0, fontsize=12)
#axes[1].set_xlim(0.98e5,1.03e5)
axes[1].xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
axes[1].ticklabel_format(style='sci', axis='x', scilimits=(-3,3))


mulify(fig, delx=0., dely=0.)

fig.tight_layout()
fig.savefig('/home/marialei/AMBER_RadCor/Figures/'+savename+'.pdf')
plt.draw()
