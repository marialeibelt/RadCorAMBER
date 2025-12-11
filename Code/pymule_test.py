from pymule import *

def OnePlot(xlabel,ylabel,xscale,yscale,nxaxis,xlabel_secondary,xscale_secondary,Savename):
    ax1.grid(True,alpha=0.7,ls='dotted')
    ax1.legend(fontsize = 9)
    ax1.set_xlabel(xlabel, fontsize=9)
    ax1.set_ylabel(ylabel,fontsize=9)
    if xscale == "log":
        ax1.set_xscale("log")
    if yscale == "log":
        ax1.set_yscale("log")
    if nxaxis == 2:
        ax1_secondary = ax1.twiny()
        ax1_secondary.set_xlabel(xlabel_secondary, fontsize=9)
        xlim1_secondary = secondary_xaxis_transform(np.array(ax1.get_xlim()))
        ax1_secondary.set_xlim(xlim1_secondary)
    if xscale_secondary == "log":
        ax1_secondary.set_xscale("log")
    plt.rcParams.update({'font.size': 7, 'font.family' : 'serif'})
    plt.tight_layout()
    fig.savefig("/home/marialei/AMBER_RadCor/Figures/" + Savename + ".pdf", bbox_inches = "tight")
    fig.savefig("/home/marialei/AMBER_RadCor/Figures/" + Savename + ".png", dpi=700, bbox_inches = "tight")
    plt.show()


plt.rcParams['text.usetex'] = False


setup(folder="/home/marialei/mcmule-release/out") 
#correct out folder mp2mp

lo = mergefks(sigma('mp2mp0'))*alpha**2*conv 
#nlo = mergefks(sigma('mp2mpR'),sigma('mp2mpF'),sigma('mp2mpA'))*alpha**3*conv

print(lo.value)
#print(lo["Emu"],lo["th3"])



#plotting
fig, axs = plt.subplots(1,2,1,col=blue)
plt.sca(axs)
obs = 'th3'

errorband(lo[obs])

axs.legend(['something'], loc='lower center',framealpha=0)
axs.set_title('title')
plt.xlim(0., 0.)
plt.ylim(0., 0.)
#mulify(fig, delx=0., dely=0.)

plt.draw()

fig.tight_layout()

fig.savefig('example.pdf')