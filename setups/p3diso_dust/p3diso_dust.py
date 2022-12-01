import sys, os
try:
    FARGO_DIR = os.environ["FARGO_PATH"]
    sys.path.append(os.path.join(FARGO_DIR, "python"))
except KeyError:
    pass

from fargo_lib import FARGO_Sim
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


try:
    Sim = FARGO_Sim(sys.argv[1])
except IndexError:
    Sim = FARGO_Sim('outputs/p3diso_dust')
try:
    snap_num = int(sys.argv[2])
except IndexError:
    snap_num = 0

    
r  = Sim.rmed
th = Sim.thetamed

rho_g = Sim.load_field('gasdens', snap_num)
vth_g  = Sim.load_field('gasvz',  snap_num)
vr_g   = Sim.load_field('gasvy',  snap_num)

plt.ioff()

plt.subplot(311)
plt.pcolormesh(r, np.cos(th), rho_g.T, norm=LogNorm())
plt.colorbar(label=r'$\rho_g$')
plt.xlabel(r"$r$")
plt.ylabel(r"$z/R$")

plt.subplot(312)
plt.pcolormesh(r, np.cos(th), vth_g.T)
plt.colorbar(label=r'$v_\theta$')
plt.xlabel(r"$r$")
plt.ylabel(r"$z/R$")

plt.subplot(313)
plt.pcolormesh(r, np.cos(th), vr_g.T)
plt.colorbar(label=r'$v_r$')
plt.xlabel(r"$r$")
plt.ylabel(r"$z/R$")


plt.tight_layout()

for i in range(int(Sim.get_parameter('NDUST'))):
    rho = Sim.load_field('dustdens{}_'.format(i), snap_num)
    vth = Sim.load_field('dustvz{}_'.format(i),   snap_num)     
    vr  = Sim.load_field('dustvy{}_'.format(i),   snap_num)     

    plt.figure()
    
    plt.subplot(311)
    plt.pcolormesh(r, np.cos(th), rho.T, norm=LogNorm(), vmin=1e-13)
    plt.colorbar(label=r'$\rho_{{d_{}}}$'.format(i))
    plt.xlabel(r"$r$")
    plt.ylabel(r"$z/R$")

    plt.subplot(312)
    plt.pcolormesh(r, np.cos(th), vth.T)
    plt.colorbar(label=r'$v_\theta$')
    plt.xlabel(r"$r$")
    plt.ylabel(r"$z/R$")

    plt.subplot(313)
    plt.pcolormesh(r, np.cos(th), vr.T)
    plt.colorbar(label=r'$v_r$')
    plt.xlabel(r"$r$")
    plt.ylabel(r"$z/R$")

    plt.tight_layout()

plt.show()

