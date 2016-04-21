#! /usr/env ipython

from pylab import *
import numpy as np
from scipy.interpolate import interp1d

def ddy_cgrid_centered(y, Y):
    dy = np.zeros(y.shape, np.float) #we know it will be this size
    h = Y[1] - Y[0] #this assumes the points are evenely spaced!
    dy[:,2:-2] = (y[:,0:-4] - 8 * y[:,1:-3] + 8 * y[:,3:-1] - y[:,4:]) / (12.0 * h)

    # simple differences at the end-points
    dy[:,0] = (y[:,1] - y[:,0])/(Y[1] - Y[0])
    dy[:,1] = (y[:,2] - y[:,1])/(Y[2] - Y[1])
    dy[:,-2] = (y[:,-2] - y[:,-3]) / (Y[-2] - Y[-3])
    dy[:,-1] = (y[:,-1] - y[:,-2]) / (Y[-1] - Y[-2])
    return dy

def ddz_cgrid_centered(y, Z):
    dy = np.zeros(y.shape, np.float) 
    dy[0] = (y[0] - y[1])/(Z[0] - Z[1])
    for i in range(1,len(y)-1):
        dy[i] = (y[i+1] - y[i-1])/(Z[i+1]-Z[i-1])
    dy[-1] = (y[-1] - y[-2])/(Z[-1] - Z[-2])
    return dy

def ddz_cgrid_centered_any(y,ZF):
    dy = np.zeros(y.shape, np.float) 
    dy[0] = (y[0] - y[1])/(ZF[0] - ZF[1])
    for i in range(1,len(y)-1):
        dy[i] = (y[i+1] - y[i-1])/(ZF[i+1]-ZF[i-1])
    dy[-1] = (y[-1] - y[-2])/(ZF[-1] - ZF[-2])
    return dy

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def get_psi_iso_z(PSISO, Th, Z, depth=2985):
    """Put the output from psi_iso into Z coordinates."""

    # figure out the depth of each layer
    h = Th
    # psi_iso is defined at the *bottom* of each layer,
    # therefore we want the depth at the bottom of the layer
    z = np.cumsum(h, axis=0) - depth
    # interpolate to center z points
    nz = len(Z)
    ny = len(Th[1,:])
    psi_iso_z = zeros((nz,ny))
    for j in arange(ny-1):
        psi_iso_z[:,j] = interp(Z[:],z[:,j], PSISO[:,j])
    return psi_iso_z
