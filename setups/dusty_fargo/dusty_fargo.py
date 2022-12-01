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

    
r   = Sim.Rmed
phi = Sim.phimed 

rho_g  = Sim.load_field('gasdens', snap_num)
vr_g   = Sim.load_field('gasvy',  snap_num)
vphi_g = Sim.load_field('gasvx',  snap_num)

vphi_g =  (vphi_g + Sim.get_parameter("OMEGAFRAME") * r) - r**-0.5

plt.ioff()

plt.subplot(311)
plt.pcolormesh(r, phi/np.pi, rho_g, norm=LogNorm())
plt.colorbar(label=r'$\rho_g$')
plt.xlabel(r"$r$")
plt.ylabel(r"$\phi$")

plt.subplot(312)
plt.pcolormesh(r, phi/np.pi, vr_g)
plt.colorbar(label=r'$v_r$')
plt.xlabel(r"$r$")
plt.ylabel(r"$\phi$")

plt.subplot(313)
plt.pcolormesh(r, phi/np.pi, vphi_g)
plt.colorbar(label=r'$\delta v_\phi$')
plt.xlabel(r"$r$")
plt.ylabel(r"$\phi$")

plt.tight_layout()

for i in range(int(Sim.get_parameter('NDUST'))):
    rho  = Sim.load_field('dustdens{}_'.format(i), snap_num)
    vr   = Sim.load_field('dustvy{}_'.format(i),   snap_num)     
    vphi = Sim.load_field('dustvx{}_'.format(i),   snap_num)

    vphi =  (vphi + Sim.get_parameter("OMEGAFRAME") * r) - r**-0.5

    plt.figure()
    
    plt.subplot(311)
    plt.pcolormesh(r, phi/np.pi, rho, norm=LogNorm())
    plt.colorbar(label=r'$\rho_{{d_{}}}$'.format(i))
    plt.xlabel(r"$r$")
    plt.ylabel(r"$\phi$")

    plt.subplot(312)
    plt.pcolormesh(r, phi/np.pi, vr)
    plt.colorbar(label=r'$v_r$')
    plt.xlabel(r"$r$")
    plt.ylabel(r"$\phi$")

    plt.subplot(313)
    plt.pcolormesh(r, phi/np.pi, vphi)
    plt.colorbar(label=r'$\delta v_\phi$')
    plt.xlabel(r"$r$")
    plt.ylabel(r"$\phi$")

    plt.tight_layout()

plt.show()

