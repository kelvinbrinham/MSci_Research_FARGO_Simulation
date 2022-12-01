import numpy as np
import matplotlib.pyplot as plt
import sys
import os

##################
# Read in the grid info
#################
with open('dimensions.dat') as f:
    HEAD = f.readline().strip().split()
    for i, s in enumerate(HEAD):
        if 'NGHY' in s:
            Nghost = int(f.readline().strip().split()[i])
    
r   = np.genfromtxt('domain_y.dat')[Nghost:-Nghost]
rMin = r[:-1]
rMed   = 0.5*(r[1:] + r[:-1])
phi =  np.genfromtxt('domain_x.dat')
phiMin = phi[:-1]
phiMed = 0.5*(phi[1:] + phi[:-1])


Nr   = len(rMin)
Nphi = len(phiMin)

Omega = rMin**-1.5

##################
# Get the Stoke's Number
##################
St = 1. / np.genfromtxt('drag_coefficients.dat')
try:
    NDust = len(St)
except TypeError:
    NDust = 1
    St = [St,]


#################
# Find out which snapshot we want
#################
SnapNum = int(sys.argv[1])

################
# Read in the radial velocity data
################
vr_g = np.fromfile('gasvy'+str(SnapNum)+'.dat').reshape(Nr, Nphi).mean(1)
vp_g = np.fromfile('gasvx'+str(SnapNum)+'.dat').reshape(Nr, Nphi).mean(1)

v_d = [] 
vp_d = []
for i in range(NDust):
    vdi = np.fromfile('dustvy'+str(i)+'_'+str(SnapNum)+'.dat')
    v_d.append(vdi.reshape(Nr, Nphi).mean(1))

    vdi = np.fromfile('dustvx'+str(i)+'_'+str(SnapNum)+'.dat')
    vdi = vdi.reshape(Nr, Nphi).mean(1)
    vp_d.append(vdi)


###############
# Construct the analytical solution
###############
cs    = np.fromfile('gasenergy'+str(SnapNum)+'.dat').reshape(Nr, Nphi)[:,0]
Sigma = np.fromfile('gasdens'  +str(SnapNum)+'.dat').reshape(Nr, Nphi)[:,0]

H = cs / Omega
# Get alpha
Alpha = 0
for l in open('variables.par'):
    if 'NU' in l:
        Nu = float(l.strip().split()[1])
    elif 'ALPHA' in l:
        Alpha = float(l.strip().split()[1])
    elif 'FLARINGINDEX' in l:
        kFlare = float(l.strip().split()[1])
    elif 'SIGMASLOPE'in l:
        kSigma = float(l.strip().split()[1])
    elif 'ASPECTRATIO' in l:
        H_R = float(l.strip().split()[1])

H = H_R * rMin**(1 + kFlare)
cs = H * rMin**-1.5

q = 2*kFlare - 1
p = kSigma - (3 + q) / 2.

if Alpha:
    Nu = Alpha * cs * H

kDrift = (1 - 2*kFlare + kSigma)

vr_visc =  -3 * (1 - kSigma + 2*kFlare) * Nu / rMed
vr_drift = []
vp_drift = []
for St_i in St:
    vd = kDrift * H_R**2 * rMin[1:]**(2*kFlare-0.5)
    vnu = 0.5*(vr_visc[1:] + vr_visc[:-1])
    vr_drift.append((vnu/St_i - vd)/(St_i + 1./ St_i))
    
    vp_drift.append((0.5*vd*St_i - 0.5*vnu)/(St_i + 1./ St_i))
    

##############
# Plot the velocities
#############
plt.subplot(131)
l, = plt.semilogx(rMin, vr_g)
plt.semilogx(rMed, vr_visc,l.get_color() +'--')

plt.subplot(132)
for n in range(NDust):
    ax, = plt.semilogx(rMin, v_d[n], label=str(St[n]))
    plt.semilogx(rMin[1:], vr_drift[n], ax.get_color()+'--')

plt.subplot(133)
for n in range(NDust):
    ax, = plt.semilogx(rMin, vp_d[n]-vp_g, label=str(St[n]))
    plt.semilogx(rMin[1:], vp_drift[n], ax.get_color()+'--')
    plt.semilogx(rMin, -v_d[n]*St[n]/2, ax.get_color()+':')

plt.legend(loc='lower right')
plt.show()


