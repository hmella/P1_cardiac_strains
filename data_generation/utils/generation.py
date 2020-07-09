import os

import matplotlib.pyplot as plt
from PyMRStrain import *

from utils.im_3D_parameters import (echo_train_length, filter3, filter_lift3,
                                    filter_width3, off_resonance, receiver_bw,
                                    spatial_shift)
from utils.im_parameters import filter, filter_lift, filter_width


# Parameters generation
def generate_phantoms(nb_samples,ini=0,fin=0):

    # Create folders
    if not os.path.isdir('inputs/parameters'):
        os.makedirs('inputs/parameters',exist_ok=True)
    if not os.path.isdir('inputs/spins'):
        os.makedirs('inputs/spins',exist_ok=True)

    for d in range(ini,fin):

        # Create parameters
        p = Parameters(time_steps=20)
        p.h = 0.008
        p.sigma = np.random.uniform(0.5, 1.7)
        p.xi = np.random.uniform(0.5,1.5)

        # Create spins
        s = Spins(Nb_samples=nb_samples,parameters=p)

        # Write parameters and spins
        save_pyobject(p,'inputs/parameters/p_{:03d}.pkl'.format(d))
        save_pyobject(s,'inputs/spins/spins_{:03d}.pkl'.format(d), sep_proc=True)


# CSPAMM images generation
def generate_cspamm(resolutions, frequencies, patients, ini=0, fin=0):

    # Create folder
    if not os.path.isdir('inputs/kspaces'):
        os.makedirs('inputs/kspaces',exist_ok=True)
    if not os.path.isdir('inputs/masks'):
        os.makedirs('inputs/masks',exist_ok=True)

    # Resolutions loop
    for (rn, r) in enumerate(resolutions):

        # Create image
        I = CSPAMMImage(FOV=np.array([0.1, 0.1, 0.008]),
                  center=np.array([0.0,0.0,0.0]),
                  resolution=r,
                  encoding_frequency=np.array([0,0,0]),
                  T1=[1e-10,1e-10,0.85],
                  M0=[0,0,1],
                  flip_angle=15*np.pi/180,
                  encoding_angle=90*np.pi/180,
                  kspace_factor=2,
                  slice_thickness=0.008,
                  oversampling_factor=2,
                  phase_profiles=r[1])

        # Filter specifications
        I.filter       = filter
        I.filter_width = filter_width
        I.filter_lift  = filter_lift

        # Frequencies loop
        for (fn, f) in enumerate(frequencies):

            # Set encoding frequency
            I.encoding_frequency = f
         
            # Data loop
            for d in range(ini,fin):

                # Debug
                MPI_print('[CSPAMM] Generating freq. {:d}, data {:d}'.format(fn,d))

                # Load spins and parameter objects
                spins = load_pyobject('inputs/spins/spins_{:03d}.pkl'.format(d),sep_proc=True)
                param = load_pyobject('inputs/parameters/p_{:03d}.pkl'.format(d))

                # Create phantom
                phantom = Phantom(spins, param, patient=patients[d], z_motion=False)

                # Artifact
                artifact = None

                # Generate kspaces
                NSA_1, NSA_2, mask = I.generate(artifact, phantom, param, debug=False)

                plt.imshow(np.abs(NSA_1.k[...,0,1,0] - NSA_2.k[...,0,1,0]))
                plt.show(block=False)
                plt.pause(15)

                # Compress kspaces
                NSA_1.scale(dtype=np.uint16)
                NSA_2.scale(dtype=np.uint16)

                # Get mask
                maskim = mask.to_img()
                maskim = np.abs(maskim) > 0.25*np.abs(maskim).max()

                # Export kspaces and masks
                save_pyobject([NSA_1,NSA_2],'inputs/kspaces/CI_{:03d}_{:02d}_{:02d}.pkl'.format(d,fn,rn))
                save_pyobject(maskim,'inputs/masks/CI_{:03d}_{:02d}_{:02d}.pkl'.format(d,fn,rn))


# DENSE images generation
def generate_dense(resolutions, frequencies, patients, ini=0, fin=0):

    # Create folder
    if not os.path.isdir('inputs/kspaces'):
        os.makedirs('inputs/kspaces',exist_ok=True)
    if not os.path.isdir('inputs/masks'):
        os.makedirs('inputs/masks',exist_ok=True)

    # Resolutions loop
    for (rn, r) in enumerate(resolutions):

        # Create image
        I = DENSEImage(FOV=np.array([0.1, 0.1, 0.008]),
                  center=np.array([0.0,0.0,0.0]),
                  resolution=r,
                  encoding_frequency=np.array([0,0,0]),
                  T1=[1e-10,1e-10,0.85],
                  M0=[0,0,1],
                  flip_angle=15*np.pi/180,
                  kspace_factor=2,
                  slice_thickness=0.008,
                  oversampling_factor=2,
                  phase_profiles=r[1])

        # Filter specifications
        I.filter       = filter
        I.filter_width = filter_width
        I.filter_lift  = filter_lift

        # Frequencies loop
        for (fn, f) in enumerate(frequencies):

            # Set encoding frequency
            I.encoding_frequency = f
         
            # Data loop
            for d in range(ini,fin):

                # Debug
                MPI_print('[DENSE] Generating freq. {:d}, data {:d}'.format(fn,d))

                # Load spins and parameter objects
                spins = load_pyobject('inputs/spins/spins_{:03d}.pkl'.format(d),sep_proc=True)
                param = load_pyobject('inputs/parameters/p_{:03d}.pkl'.format(d))

                # Create phantom
                phantom = Phantom(spins, param, patient=patients[d], z_motion=False)

                # Artifact
                artifact = None

                # Generate kspaces
                NSA_1, NSA_2, mask = I.generate(artifact, phantom, param, debug=False)

                # Compress kspaces
                NSA_1.scale(dtype=np.uint16)
                NSA_2.scale(dtype=np.uint16)

                # Get mask
                maskim = mask.to_img()
                maskim = np.abs(maskim) > 0.25*np.abs(maskim).max()

                # Export kspaces
                save_pyobject([NSA_1,NSA_2],'inputs/kspaces/DI_{:03d}_{:02d}_{:02d}.pkl'.format(d,fn,rn))
                save_pyobject(maskim,'inputs/masks/DI_{:03d}_{:02d}_{:02d}.pkl'.format(d,fn,rn))


# REFERENCE images generation
def generate_reference(resolutions, frequencies, patients, ini=0, fin=0):

    # Create folder
    if not os.path.isdir('inputs/kspaces'):
        os.makedirs('inputs/kspaces',exist_ok=True)
    if not os.path.isdir('inputs/masks'):
        os.makedirs('inputs/masks',exist_ok=True)

    # Resolutions loop
    for (rn, r) in enumerate(resolutions):

        # Create image
        I = EXACTImage(FOV=np.array([0.1, 0.1, 0.008]),
                    center=np.array([0.0,0.0,0.0]),
                    resolution=r,
                    encoding_frequency=np.array([100.0,100.0,0.0]),
                    kspace_factor=2,
                    slice_thickness=0.008,
                    oversampling_factor=2,
                    phase_profiles=r[1])        

        # Filter specifications
        I.filter       = filter
        I.filter_width = filter_width
        I.filter_lift  = filter_lift

        # Frequencies loop
        for (fn, f) in enumerate([np.array([100.0,100.0,0.0])]):
         
            # Data loop
            for d in range(ini,fin):

                # Debug
                MPI_print('[REFERENCE] Generating freq. {:d}, data {:d}'.format(fn,d))

                # Load spins and parameter objects
                spins = load_pyobject('inputs/spins/spins_{:03d}.pkl'.format(d),sep_proc=True)
                param = load_pyobject('inputs/parameters/p_{:03d}.pkl'.format(d))

                # Create phantom
                phantom = Phantom(spins, param, patient=patients[d], z_motion=False)

                # Artifact
                artifact = None

                # Generate kspaces
                NSA_1, mask = I.generate(artifact, phantom, param, debug=False)

                # Compress kspaces
                NSA_1.scale(dtype=np.uint16)

                # Get mask
                maskim = mask.to_img()
                maskim = np.abs(maskim) > 0.25*np.abs(maskim).max()

                # Export kspaces
                save_pyobject(NSA_1,'inputs/kspaces/EI_{:03d}_{:02d}_{:02d}.pkl'.format(d,fn,rn))
                save_pyobject(maskim,'inputs/masks/EI_{:03d}_{:02d}_{:02d}.pkl'.format(d,fn,rn))


# Add noise to CSPAMM
def cspamm_noisy_images(resolutions, frequencies, noise_levels, ini=0, fin=0):

    # Create folder
    if not os.path.isdir('inputs/noisy_images'):
        os.mkdir('inputs/noisy_images')

    # Resolutions loop
    for (rn, r) in enumerate(resolutions):

        # Frequencies loop
        for (fn, f) in enumerate(frequencies):
         
            # Data loop
            for d in range(ini,fin):

                # Load kspaces
                (NSA_1, NSA_2) = load_pyobject('inputs/kspaces/Ck_{:03d}_{:02d}_{:02d}.pkl'.format(d,fn,rn))

                # Rescale kspaces
                NSA_1.rescale()
                NSA_2.rescale()

                # Noise loop
                for (nn, n) in enumerate(noise_levels):

                    # Add noise to kspaces
                    NSA_1.k = add_cpx_noise(NSA_1.k, mask=NSA_1.k_msk, sigma=n)
                    NSA_2.k = add_cpx_noise(NSA_2.k, mask=NSA_2.k_msk, sigma=n)

                    # Get images and scale
                    I1 = scale_image(NSA_1.to_img(),mag=False,real=True,compl=True,dtype=np.uint16)
                    I2 = scale_image(NSA_2.to_img(),mag=False,real=True,compl=True,dtype=np.uint16)

                    # Export noisy kspaces
                    save_pyobject([I1,I2],'inputs/noisy_images/CI_{:03d}_{:02d}_{:02d}_{:02d}.pkl'.format(d,fn,rn,nn))


# Add noise to DENSE
def dense_noisy_images(resolutions, frequencies, noise_levels, ini=0, fin=0):

    # Create folder
    if not os.path.isdir('inputs/noisy_images'):
        os.mkdir('inputs/noisy_images')

    # Resolutions loop
    for (rn, r) in enumerate(resolutions):

        # Frequencies loop
        for (fn, f) in enumerate(frequencies):
         
            # Data loop
            for d in range(ini,fin):

                # Load kspaces
                (NSA_1, NSA_2) = load_pyobject('inputs/kspaces/Dk_{:03d}_{:02d}_{:02d}.pkl'.format(d,fn,rn))

                # Rescale kspaces
                NSA_1.rescale()
                NSA_2.rescale()

                # Noise loop
                for (nn, n) in enumerate(noise_levels):

                    # Add noise to kspaces
                    NSA_1.k = add_cpx_noise(NSA_1.k, mask=NSA_1.k_msk, sigma=n)
                    NSA_2.k = add_cpx_noise(NSA_2.k, mask=NSA_2.k_msk, sigma=n)

                    # Get images and scale
                    I1 = scale_image(NSA_1.to_img(),mag=False,real=True,compl=True,dtype=np.uint16)
                    I2 = scale_image(NSA_2.to_img(),mag=False,real=True,compl=True,dtype=np.uint16)

                    # Export noisy kspaces
                    save_pyobject([I1,I2],'inputs/noisy_images/DI_{:03d}_{:02d}_{:02d}_{:02d}.pkl'.format(d,fn,rn,nn))


# 3D CSPAMM generation
def generate_3D_cspamm(parameters, imaging_parameters, phantom, resolution, ke, std_noise=None):

    # Create folder
    if not os.path.isdir('inputs/3D_experiments'):
        os.makedirs('inputs/3D_experiments',exist_ok=True)

    # Generate images:
    for v, view in enumerate(['base','mid','apex']):

        # Create complementary CSPAMM image
        I = CSPAMMImage(FOV=imaging_parameters['FOV'][v],
                  center=imaging_parameters['center'][v],
                  resolution=resolution,
                  encoding_frequency=np.array([ke,ke,0]),
                  T1=imaging_parameters['T1'],
                  M0=imaging_parameters['M0'],
                  flip_angle=imaging_parameters['flip_angle'],
                  encoding_angle=90*np.pi/180,
                  off_resonance=imaging_parameters['phi'],
                  kspace_factor=2,
                  slice_thickness=imaging_parameters['slice_thickness'],
                  slice_offset=imaging_parameters['offset'][v],
                  oversampling_factor=2,
                  phase_profiles=135)

        # EPI acquisiton object
        epi = EPI(receiver_bw=receiver_bw,
                  echo_train_length=echo_train_length,
                  off_resonance=off_resonance,
                  acq_matrix=I.acq_matrix,
                  spatial_shift=spatial_shift)

        # Filter specifications
        I.filter       = filter3
        I.filter_width = filter_width3
        I.filter_lift  = filter_lift3

        # Generate images
        NSA_1, NSA_2, mask = I.generate(epi, phantom, parameters, debug=True)

        # Add noise to DENSE images
        if std_noise is not None:
            noise_1 = np.random.normal(0, std_noise, NSA_1.k.shape) \
                    + 1j*np.random.normal(0, std_noise, NSA_1.k.shape)
            noise_2 = np.random.normal(0, std_noise, NSA_2.k.shape) \
                    + 1j*np.random.normal(0, std_noise, NSA_2.k.shape)
            f_noise_1 = NSA_1.k_msk*itok(noise_1)
            f_noise_2 = NSA_2.k_msk*itok(noise_2)
            NSA_1.k += f_noise_1
            NSA_2.k += f_noise_2

        # Get images
        In1 = NSA_1.to_img()
        In2 = NSA_2.to_img() 

        # Get mask
        maskim = mask.to_img()
        maskim = np.abs(maskim) > 0.25*np.abs(maskim).max()

        # Export noisy kspaces
        save_pyobject([In1,In2,maskim],'inputs/3D_experiments/CI_{:s}.pkl'.format(view))

        # # Plot
        # I = In1 - In2
        # if MPI_rank==0:
        #     # multi_slice_viewer(np.abs(I[:,:,0,0,:]))
        #     # multi_slice_viewer(np.angle(I[:,:,0,0,:]))
        #     # multi_slice_viewer(np.abs(itok(I[:,:,0,0,:])))
        #     # multi_slice_viewer(np.abs(NSA_1.k_acq[:,:,0,0,:]),clim=[0, 2000])
        #     # multi_slice_viewer(np.abs(NSA_1.k[:,:,0,0,:]),clim=[0, 2000])

    return std_noise

# 3D CSPAMM generation
def generate_3D_dense(parameters, imaging_parameters, phantom, resolution, ke,
                      rel_std=0.01, std_noise=None):

    # Create folder
    if not os.path.isdir('inputs/3D_experiments'):
        os.makedirs('inputs/3D_experiments',exist_ok=True)

    # Generate images:
    for v, view in enumerate(['base','mid','apex']):

        # Create complementary CSPAMM image
        I = DENSEImage(FOV=imaging_parameters['FOV'][v],
                  center=imaging_parameters['center'][v],
                  resolution=resolution,
                  encoding_frequency=np.array([ke,ke,0]),
                  T1=imaging_parameters['T1'],
                  M0=imaging_parameters['M0'],
                  flip_angle=imaging_parameters['flip_angle'],
                  off_resonance=imaging_parameters['phi'],
                  kspace_factor=2,
                  slice_thickness=imaging_parameters['slice_thickness'],
                  slice_offset=imaging_parameters['offset'][v],
                  oversampling_factor=2,
                  phase_profiles=63)

        # EPI acquisiton objects
        epi = EPI(receiver_bw=receiver_bw,
                  echo_train_length=echo_train_length,
                  off_resonance=off_resonance,
                  acq_matrix=I.acq_matrix,
                  spatial_shift=spatial_shift)

        # Filter specifications
        I.filter       = filter3
        I.filter_width = filter_width3
        I.filter_lift  = filter_lift3

        # Generate images
        NSA_1, NSA_2, mask = I.generate(epi, phantom, parameters, debug=True)

        # Get noise standard deviation
        I_max_r0 = max([np.abs(ktoi(NSA_1.k)[...,0,0]).max(), np.abs(ktoi(NSA_1.k)[...,1,0]).max()])
        std_noise = I_max_r0*rel_std

        # Add noise to DENSE images
        noise_1 = np.random.normal(0, std_noise, NSA_1.k.shape) \
                + 1j*np.random.normal(0, std_noise, NSA_1.k.shape)
        noise_2 = np.random.normal(0, std_noise, NSA_2.k.shape) \
                + 1j*np.random.normal(0, std_noise, NSA_2.k.shape)
        f_noise_1 = NSA_1.k_msk*itok(noise_1)
        f_noise_2 = NSA_2.k_msk*itok(noise_2)
        NSA_1.k += f_noise_1
        NSA_2.k += f_noise_2

        # Get images
        In1 = NSA_1.to_img()
        In2 = NSA_2.to_img() 

        # Get mask
        maskim = mask.to_img()
        maskim = np.abs(maskim) > 0.25*np.abs(maskim).max()

        # Export noisy kspaces
        save_pyobject([In1,In2,maskim],'inputs/3D_experiments/DI_{:s}.pkl'.format(view))

        # # Plot
        # I = In1 - In2
        # if MPI_rank==0:

        #     # # Contrast
        #     # maxl = []
        #     # for i in range(20):
        #     #     L = np.abs(I[:,:,0,0,i]).flatten()
        #     #     maxl.append(L.max())
        #     #     print(i, L.max())
        #     # plt.plot(maxl)
        #     # plt.show()

        #     # multi_slice_viewer(np.abs(I[:,:,0,0,:]))
        #     # multi_slice_viewer(np.angle(I[:,:,0,0,:]))
        #     # multi_slice_viewer(np.abs(NSA_1.k_acq[:,:,0,0,:]),clim=[0, 2000])
        #     # multi_slice_viewer(np.abs(NSA_1.k[:,:,0,0,:]),clim=[0, 2000])


    return std_noise

# 3D reference
def generate_3D_ref(parameters, imaging_parameters, phantom, resolution, name='I'):

    # Create folder
    if not os.path.isdir('inputs/3D_experiments'):
        os.makedirs('inputs/3D_experiments',exist_ok=True)

    # Generate images:
    for v, view in enumerate(['base','mid','apex']):

        # Create reference image
        I = EXACTImage(FOV=imaging_parameters['FOV'][v],
                  center=imaging_parameters['center'][v],
                  resolution=resolution,
                  encoding_frequency=np.array([100,100,0]),
                  T1=imaging_parameters['T1'],
                  M0=imaging_parameters['M0'],
                  kspace_factor=2,
                  slice_thickness=imaging_parameters['slice_thickness'],
                  slice_offset=imaging_parameters['offset'][v],
                  oversampling_factor=2,
                  phase_profiles=resolution[1])

        # Filter specifications
        I.filter       = filter3
        I.filter_width = filter_width3
        I.filter_lift  = filter_lift3

        # Generate images
        NSA_1, mask = I.generate(None, phantom, parameters, debug=True)

        # Get images
        I = NSA_1.to_img()

        # Get mask
        maskim = mask.to_img()
        maskim = np.abs(maskim) > 0.25*np.abs(maskim).max()

        # Export noisy kspaces
        save_pyobject([I,maskim],'inputs/3D_experiments/{:s}_{:s}.pkl'.format(name,view))

        # # Plot
        # if MPI_rank==0:
        #     multi_slice_viewer(np.abs(I[:,:,0,0,:]))
        #     multi_slice_viewer(np.angle(I[:,:,0,0,:]))
        #     multi_slice_viewer(np.abs(itok(np.abs(I[:,:,0,0,:]))))


    return 1
