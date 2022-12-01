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
    Sim = FARGO_Sim('outputs/vertical_equilibrium')
try:
    snap_num = int(sys.argv[2])
except IndexError:
    snap_num = 0

    
r  = Sim.rmed
th = Sim.thetamed

rho_g = Sim.load_field('gasdens', snap_num)
vth_g = Sim.load_field('gasvz',   snap_num)

plt.ioff()


plt.subplot(211)
plt.pcolormesh(r, th/np.pi-0.5, rho_g.T, norm=LogNorm())
plt.colorbar(label=r'$\rho_g$')
plt.xlabel(r"$r$")
plt.ylabel(r"$z/R$")

plt.subplot(212)
plt.pcolormesh(r, th/np.pi-0.5, vth_g.T)
plt.colorbar(label=r'$v_\theta$')
plt.xlabel(r"$r$")
plt.ylabel(r"$z/R$")

plt.tight_layout()

plt.show()

