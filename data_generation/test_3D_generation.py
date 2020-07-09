import time

import numpy as np
from PyMRStrain import *

from utils.generation import (generate_3D_cspamm, generate_3D_dense,
                              generate_3D_ref)
from utils.im_3D_parameters import (FOV, M0, T1, center, keC, keD,
                                    noise_levels, offset, resC, resD,
                                    slice_thickness)

if __name__=="__main__":

    # Parameters
    p = Parameters(R_en=0.025, R_ep=0.035, tau=0.01, time_steps=20, 
                  phi_en=-10*np.pi/180, phi_ep=-4*np.pi/180,
                  S_ar=1.1, S_en=0.7, psi=0, xi=0.5, sigma=2,
                  R_inner=0.015, R_outer=0.045, tA=0.15, tB=0.35, tC=0.5)
    save_pyobject(p, 'p.pkl')
    p=load_pyobject('p.pkl')

    # Flip angles
    # alpha_n = 21*np.pi/180  # flip angle at the last time step
    # TR = 1.0/p.time_steps   # repetition time
    # C  = np.exp(-TR/T1[2])  # constant
    # flip_angle = np.zeros([p.time_steps],dtype=np.float)
    # flip_angle[:] = alpha_n
    # for i in range (p.time_steps-2,-1,-1):
    #   flip_angle[i] = np.arctan(np.sin(alpha_n)*C)
    #   alpha_n = flip_angle[i]
    flip_angle = 15*np.pi/180

    # Field inhomogeneity
    phi = lambda X, Y: 0*(X+Y)/0.1*0.2

    # Imaging parameters
    iparam = {'FOV': FOV, 'center': center, 'offset': offset, 'T1': T1,
              'M0': M0, 'slice_thickness': slice_thickness, 'phi': phi,
              'flip_angle': flip_angle}

    # Spins
    spins = Spins(Nb_samples=1000000, parameters=p)

    # Create phantom object
    phantom = Phantom(spins, p, patient=False, write_vtk=False,
                      phi_en_apex=20*np.pi/180, z_motion=True)

    # Generate images
    std_noise = generate_3D_dense(p, iparam, phantom, resD, keD, rel_std=noise_levels[1])
    generate_3D_ref(p, iparam, phantom, resD, name='RDI')
    std_noise = generate_3D_cspamm(p, iparam, phantom, resC, keC, std_noise=std_noise)
    generate_3D_ref(p, iparam, phantom, resC, name='RCI')
