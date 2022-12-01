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
from scipy.fftpack import fftn, ifftn
from scipy.special import ellipk, ellipe
from scipy.integrate import quad, dblquad
from scipy.interpolate import RegularGridInterpolator


def gradient(y, x):
    grad = np.empty_like(y)

    h = np.diff(x)
    hd = h[:-1]
    hs = h[1:]

    grad[1:-1] = \
        (hs**2*y[2:] + (hd**2 - hs**2)*y[1:-1] - hd**2*y[:-2])/(hs*hd*(hs+hs))
    
    grad[0]  = (y[ 1] - y[ 0]) / (x[ 1] - x[ 0])
    grad[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])

    return grad


def kern(u, phi, eps, side):
    '''FFT kernel for disc self-gravity'''
    
    if side == 'forward':
        soft = eps*eps * np.exp(-u)
    elif side == 'backward':
        soft = eps*eps * np.exp( u) 
    elif side == 'symmetric':
        soft = 0.5*eps*eps*np.cosh(u)
    else:
        raise ValueError("Soft. parameter side={} not recognised".format(side))
    
    kern = 1. / np.sqrt(2*np.cosh(u) - 2*np.cos(phi) + soft)

    kern[kern > 1e5] = 0
    return kern

def compute_SG_potential(r, phi, sigma,eps=0.01, side='backward'):

    # Pad the arrays for fft    
    r0 = r[0]
    
    Nr = len(r)
    u = np.log(r/r0)

    norm = (u[1] - u[0])*(phi[1]-phi[0]) 

    umax = u[-1] + u[1] 
    u = np.concatenate([u, u - umax])
    u, phi = np.meshgrid(u, phi, indexing='ij')  

    Sigma = np.concatenate([sigma, np.zeros_like(sigma)])

    fft_kern = fftn(kern(u, phi, eps, side))
    fft_sig  = fftn(Sigma*np.exp(1.5*u))

    return - np.exp(-u/2)*ifftn(fft_sig*fft_kern).real * norm * r0


def SG_pot_direct(r0, r1, ri, Sigma, eps=0.01):

    def k(u, norm):
        csh = np.cosh(u)
        soft = eps*eps*np.exp(u)
        return lambda phi: norm / np.sqrt(2*(csh - np.cos(phi)) + soft)

    ui = np.log(ri)
    def fr(u):
        f0 = Sigma(np.exp(u)) * np.exp(1.5*u)
        return quad(k(ui-u, f0), 0, 2*np.pi)[0]
    
    return -quad(fr, np.log(r0), np.log(r1))[0] / np.sqrt(ri)


def interp_to_grid(R, theta, data, Npts=1025):

    Rmax = R[-1]
    Nr = len(R)
    Nt = len(theta)

    x = np.linspace(-Rmax, Rmax, Npts)
    y = np.linspace(-Rmax, Rmax, Npts)

    xs, ys = np.meshgrid(x, y)

    Rgrid = (xs**2 + ys**2)**0.5
    th_grid = np.arctan2(ys, xs)

    # Extend the density grid to handle periodicity in theta:
    values = np.empty([Nr, Nt + 2], dtype='f8')
    values[:, 1:-1] = data.T
    values[:,0] = data[-1]
    values[:,-1] = data[0]
    th = np.concatenate( [[theta[-1]-2*np.pi], theta, [theta[0] + 2*np.pi]])
    
    # Interpolate the points
    interp = RegularGridInterpolator([R, th], values,
                                     bounds_error=False, fill_value=0)

    grid_pts = np.array([Rgrid.ravel(), th_grid.ravel()]).T
    return x, y, interp(grid_pts).reshape(Npts, Npts)


if __name__ == "__main__":


    argv = sys.argv[1:]

    try:
        i = argv.index('-d')
        DIR = argv[i+1]
        argv = argv[:i] + argv[i+2:]
    except ValueError:
        DIR = './'


    try:
        snap_num = int(argv[0])
    except IndexError:
        snap_num = 0

    Sim = FARGO_Sim(DIR)


    r = Sim.Rmed
    phi = Sim.phimed 

    rho_g  = Sim.load_field('gasdens', snap_num)
    e      = Sim.load_field('gasenergy', snap_num)
    vr     = Sim.load_field('gasvy', snap_num)
    vphi   = Sim.load_field('gasvx', snap_num) + r
    pot    = Sim.load_field('potential', snap_num)

    gamma = Sim.get_parameter('GAMMA')

    Q = (e*(gamma-1)*gamma/rho_g)**0.5 * r**-1.5 / (np.pi * rho_g)

    print('Mean Toomre Q:\n\t', Q.mean())
    dV = np.outer(np.diff(Sim.phi), 0.5*np.diff(Sim.R**2))
    print('Total Mass:\n\t',  (dV*rho_g).sum())
    

    Sigma0 = Sim.get_parameter('SIGMA0')
    R_crit = Sim.get_parameter('SIGMACUTOFF')
    p_Sig  = Sim.get_parameter('SIGMASLOPE')    
    def Sigma(R):
        return Sigma0 *  R**-p_Sig * np.exp(-R/R_crit) 
    
    eps = Sim.get_parameter('SELFGRAVITYSOFTENING')
    
    plt.figure(1)
    plt.subplot(411)
    plt.loglog(r, rho_g.mean(0), label='gas')
    plt.loglog(r, Sigma(r))
    plt.xlabel('r') ; plt.ylabel('Sigma')

    plt.subplot(412)
    plt.loglog(r, Q.mean(0))
    plt.xlabel('r') ; plt.ylabel('Q')


    pot_fft = compute_SG_potential(r, phi, rho_g.T, eps)[:len(r)]

    rd = np.logspace(np.log10(r[0]), np.log10(r[-1]), 25)
    r0, r1 = Sim.R[0], Sim.R[-1]
    pot_direct = \
        np.array([SG_pot_direct(r0, r1, ri, Sigma, eps) for ri in rd])

    plt.subplot(413)
    gr  = - np.diff(pot+r**-1, axis=1).mean(0)/np.diff(r)
    gr_fft =  - np.diff(pot_fft, axis=0).mean(1)/np.diff(r)
    plt.loglog(Sim.R[1:-1],  gr, 'b--')
    plt.loglog(Sim.R[1:-1], -gr, 'b-')
    plt.loglog(Sim.R[1:-1],  gr_fft, 'k--')
    plt.loglog(Sim.R[1:-1], -gr_fft, 'k-')
    plt.xlabel('r') ; plt.ylabel('g_r')
    
    plt.subplot(414)
    plt.plot(r, vphi.mean(0) - r**-0.5)
    plt.xlabel('r') ; plt.ylabel('$\delta v_\phi$')
    
    plt.figure()
    plt.subplot(411)
    plt.pcolormesh(r, phi/np.pi, rho_g, 
                   vmin=max(1e-3*rho_g.max(), rho_g.min()), norm=LogNorm()) 
    plt.colorbar(label='$\Sigma$')
    plt.xscale('log') ; plt.xlim(r0, r1)
    plt.subplot(412)
    plt.pcolormesh(r, phi/np.pi, Q, vmax=min(Q.max(),30), norm=LogNorm())
    plt.colorbar(label='$Q$')
    plt.xscale('log') ; plt.xlim(r0, r1)
    plt.subplot(413)
    plt.pcolormesh(r, phi/np.pi, vr)
    plt.colorbar(label='$v_r$')
    plt.xscale('log') ; plt.xlim(r0, r1)
    plt.subplot(414)
    plt.pcolormesh(r, phi/np.pi, vphi - r**-0.5)
    plt.colorbar(label='$\delta v_\phi$')
    plt.xscale('log') ; plt.xlim(r0, r1)

    for i in range(int(Sim.get_parameter('NDUST'))):
        rho  = Sim.load_field('dustdens{}_'.format(i), snap_num)
        vr   = Sim.load_field('dustvy{}_'.format(i), snap_num)
        vphi = Sim.load_field('dustvx{}_'.format(i), snap_num) + r 

        plt.figure(1)
        plt.subplot(311)
        plt.loglog(r, rho.mean(0), label='dust_{}'.format(i))

        plt.figure()
        plt.subplot(311)
        plt.pcolormesh(r, phi/np.pi, rho, 
                       vmin=max(1e-3*rho_g.max(), rho_g.min()), norm=LogNorm()) 
        plt.colorbar(label='$\Sigma$')
        plt.xscale('log') ; plt.xlim(r0, r1)
        plt.subplot(312)
        plt.pcolormesh(r, phi/np.pi, vr)
        plt.colorbar(label='$v_r$')
        plt.xscale('log') ; plt.xlim(r0, r1)
        plt.subplot(313)
        plt.pcolormesh(r, phi/np.pi, vphi / r**-0.5 -1)
        plt.colorbar(label='$\delta v_\phi$')
        plt.xscale('log') ; plt.xlim(r0, r1)

    plt.figure(1)
    plt.subplot(411)
    plt.legend()

    plt.figure()
    plt.subplot(111, aspect='equal')
    x, y, Sig = interp_to_grid(r, phi, rho_g)
    plt.pcolormesh(x, y, Sig, norm=LogNorm())
    plt.xlabel('x [au]') ; plt.ylabel('y [au]')
    plt.xlim(x.min(), x.max()) ; plt.ylim(y.min(), y.max())
    plt.colorbar(label='Sigma')
    plt.xlim(x.min(), x.max())


    plt.show()
