import sys, os
try:
    FARGO_DIR = os.environ["FARGO_PATH"]
    sys.path.append(os.path.join(FARGO_DIR, "python"))
except KeyError:
    pass

from fargo_lib import FARGO_Sim
import numpy as np
import matplotlib.pyplot as plt


try:
    Sim = FARGO_Sim(sys.argv[1])
except IndexError:
    Sim = FARGO_Sim('outputs/dusty_shock')
try:
    snap_num = int(sys.argv[2])
except IndexError:
    snap_num = 5

    
x = Sim.zmed

rho_g = Sim.load_field('gasdens', snap_num)
vx_g  = Sim.load_field('gasvz',   snap_num)

plt.ioff()
plt.subplot(211)
plt.plot(x, rho_g, c='k', label='gas')
plt.plot(x, rho_g, '+', c='k')
for i in range(Sim.get_parameter("NDUST")):
    rho_d = Sim.load_field('dustdens{}_'.format(i), snap_num)
    rho_d /= Sim.get_parameter("DUST_TO_GAS")
    plt.plot(x, rho_d, label='dust_{}'.format(i))
plt.xlabel(r"$x$")
plt.ylabel(r"$\rho$")
plt.legend()

plt.subplot(212)
plt.plot(x, vx_g, c='k', label='gas')
for i in range(Sim.get_parameter("NDUST")):
    vx_d = Sim.load_field('dustvz{}_'.format(i), snap_num)
    plt.plot(x, vx_d, label='dust_{}'.format(i))
plt.xlabel(r"$x$")
plt.ylabel(r"$v_x$")

plt.show()
