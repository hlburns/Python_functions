#! /usr/env ipython                                                                                                  

from pylab import *
import numpy as np
from scipy.interpolate import interp1d
from Flat_channel import *
from useful import *
dx = 5000

def get_psi_isoz_full(lvrho,Th,Depth, Z):
    """Put the output from psi_iso into Z coordinates.
       Inclduing topography, still zonal average
    """
    VT = np.sum(lvrho*dx,axis=2)
    Depth = np.mean(Depth,axis=1)
    VTfdz = np.cumsum(VT[::-1, :], axis=0)  # sum up the water column
    psi = numba_regridy(VTfdz[::-1,:])/10**6
    # figure out the depth of each layer                                                                            
    h = Th.mean(axis=2)
    # psi_iso is defined at the *bottom* of each layer,                                                             
    # therefore we want the depth at the bottom of the layer                                                        
    z = np.cumsum(h, axis=0) 
    # interpolate to center z points                                                                                
    nz = len(Z)
    ny = len(Th[1,:,1])
    psi_iso_z = zeros((nz,ny))
    for j in arange(ny-1):
            layer_depth = z[:,j] - Depth[j]
            psi_iso_z[:,j] = interp(Z[:],layer_depth[:], psi[:,j])
    return psi_iso_z

def get_psi_isoz_full_nonzone(lvrho,Th,msk,Depth):
    """Put the output from psi_iso into Z coordinates.
       Non zonal average!!! Messier but, useful start to along stream lines!
    """
    VTfdz = np.cumsum(lvrho[::-1, :, :], axis=0)  # sum up the water column
    psi = numba_regridy(VTfdz[::-1,:,:])
    # figure out the depth of each layer                                                                            
    h = Th
    # psi_iso is defined at the *bottom* of each layer,                                                             
    # therefore we want the depth at the bottom of the layer                                                        
    z = np.cumsum(h, axis=0) 
    # interpolate to center z points                                                                                
    nz = len(Z)
    ny = len(Th[1,:,1])
    nx = len(Th[1,1,:])
    psi_iso_z = zeros((nz,ny,nx))
    for j in arange(ny-1):
        for i in arange(nx):
            layer_depth = z[:,j,i] - Depth[j,i]
            psi_iso_z[:,j,i] = interp(Z[:],layer_depth[:], psi[:,j,i])
    psi_iso_z = np.sum(msk*psi_iso_z*dx,axis=2)/10**6
    return psi_iso_z
