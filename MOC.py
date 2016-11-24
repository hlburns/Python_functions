#! /usr/env ipython
import numpy as np
import numpy.ma as ma

# Vtav = timeaveraged velocity
# dz cell height 
# dx cell width - you will need to convert from lat lon

def MOC(Vtav,dz,dx):
    '''MOC from time average Velocity, cell height and width
       MOC(Vtave, dz, dx)

    '''
    # zonally integrate V * dx 
    Vtav = np.nansum(Vtav*dx,axis = 2) 
    # multiply by dz 
    # either multiply the 1D dz array like I've done here
    # or tile the dz array to match your velocity array
    # V*np.tile(dz,(1,ny,nx)         
    psi_zone = np.apply_along_axis(np.multiply, 0, Vtav, dz)
    # integrate up the water column
    psi_vet_int = np.cumsum(-psi_zone[::-1,:], axis=0)
    # In MITgcm Psi is at my W points and I'm missing my bottom value
    # Pad with 0 as Psi is 0 at the bottom
    npad = ((0,1), (0,0))
    psi_pad = np.pad(psi_vert_int, pad_width=npad, mode='constant', 
                     constant_values=0)
    # Convert to Sv
    Psi = psi_pad/10**6
    return Psi

