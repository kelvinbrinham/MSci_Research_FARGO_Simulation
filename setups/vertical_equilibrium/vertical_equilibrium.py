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
    Sim = FARGO_Sim('outputs/vertical_equilibrium')
try:
    snap_num = int(sys.argv[2])
except IndexError:
    snap_num = 0


def equilibrium_profile(St, rho0, R=1):
    h     = R*float(Sim.get_parameter("ASPECTRATIO"))
    alpha = float(Sim.get_parameter("NU")) / (h*h*R**-1.5)

    alpha /= float(Sim.get_parameter("SCHMIDT"))

    z = Sim.zmed
    x = 0.5*z*z/(h*h)
    return rho0 * np.exp(-x - (St/alpha)*(np.exp(x) - 1))

    
x = Sim.zmed

rho_g = Sim.load_field('gasdens', snap_num)
vx_g  = Sim.load_field('gasvz',   snap_num)

plt.ioff()
plt.subplot(311)
plt.plot(x, rho_g, c='k', label='gas')
for i in range(Sim.get_parameter("NDUST")):
    rho_d = Sim.load_field('dustdens{}_'.format(i), snap_num)
    rho_d /= Sim.get_parameter("DUST_TO_GAS")
    l, = plt.plot(x, rho_d, label='dust_{}'.format(i))
    plt.plot(x, equilibrium_profile(10**(-2-i), rho_d[0]),
             c=l.get_color(), ls='--')
plt.xlabel(r"$x$")
plt.ylabel(r"$\rho$")
plt.legend()

plt.subplot(312)
plt.semilogy(x, rho_g, c='k', label='gas')
for i in range(Sim.get_parameter("NDUST")):
    rho_d = Sim.load_field('dustdens{}_'.format(i), snap_num)
    rho_d /= Sim.get_parameter("DUST_TO_GAS")
    l, = plt.plot(x, rho_d, label='dust_{}'.format(i))
    plt.plot(x, equilibrium_profile(10**(-2-i), rho_d[0]),
             c=l.get_color(), ls='--')
plt.xlabel(r"$x$")
plt.ylabel(r"$\rho$")
plt.ylim(ymin=1e-12)


plt.subplot(313)
plt.plot(x, vx_g, c='k', label='gas')
for i in range(Sim.get_parameter("NDUST")):
    vx_d = Sim.load_field('dustvz{}_'.format(i), snap_num)
    plt.plot(x, vx_d, label='dust_{}'.format(i))
plt.xlabel(r"$x$")
plt.ylabel(r"$v_x$")

plt.show()
