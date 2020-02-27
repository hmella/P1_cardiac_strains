import os

from PyMRStrain import *


# Parameters generation
def generate_phantoms(nb_samples,ini=0,fin=0):

    # Create folders
    if not os.path.isdir('inputs/'):
        os.mkdir('inputs/')
    if not os.path.isdir('inputs/parameters'):
        os.mkdir('inputs/parameters')
    if not os.path.isdir('inputs/spins'):
        os.mkdir('inputs/spins')

    for d in range(ini,fin):

        # Create parameters
        p = Parameters(time_steps=18)
        p.h = 0.008

        # Create spins
        s = Spins(Nb_samples=nb_samples,parameters=p)

        # Write parameters and spins
        save_pyobject(p,'inputs/parameters/p_{:03d}.pkl'.format(d))
        save_pyobject(s,'inputs/spins/spins_{:03d}.pkl'.format(d), sep_proc=True)


# CSPAMM images generation
def generate_cspamm(resolutions, frequencies, patients, ini=0, fin=0, noise_free=False):

    # Create folder
    if not os.path.isdir('inputs/kspaces'):
        os.mkdir('inputs/kspaces')

    # Resolutions loop
    for (rn, r) in enumerate(resolutions):

        # Create image
        I = CSPAMMImage(FOV=np.array([0.1, 0.1, 0.008]),
                  center=np.array([0.0,0.0,0.0]),
                  resolution=r,
                  encoding_frequency=np.array([0,0,0]),
                  T1=0.85,
                  flip_angle=15*np.pi/180,
                  encoding_angle=90*np.pi/180,
                  kspace_factor=15,
                  slice_thickness=0.008,
                  oversampling_factor=1,
                  phase_profiles=r[1])

        # Filter specifications
        I.filter = 'Riesz'
        I.filter_width = 0.8
        I.filter_lift = 0.3

        # Frequencies loop
        for (fn, f) in enumerate(frequencies):

            # Set encoding frequency
            I.encoding_frequency = f
         
            # Data loop
            for d in range(ini,fin):

                # Load spins and parameter objects
                spins = load_pyobject('inputs/spins/spins_{:03d}.pkl'.format(d),sep_proc=True)
                param = load_pyobject('inputs/parameters/p_{:03d}.pkl'.format(d))

                # Create phantom
                phantom = Phantom(spins, param, patient=patients[d], z_motion=False)

                # Artifact
                artifact = None

                # Generate kspaces
                NSA_1, NSA_2, REF, mask = I.generate(artifact, phantom, param, debug=False)

                # Export noise-free images
                if noise_free:
                    # Create folder
                    if not os.path.isdir('inputs/noise_free_images'):
                        os.mkdir('inputs/noise_free_images')

                    # Get images and scale
                    I1 = scale_image(NSA_1.to_img(),mag=False,real=True,compl=True)
                    I2 = scale_image(NSA_2.to_img(),mag=False,real=True,compl=True)

                    # Export images
                    save_pyobject([I1,I2],'inputs/noise_free_images/CI_{:03d}_{:02d}_{:02d}.pkl'.format(d,fn,rn))

                # Compress kspaces
                NSA_1.scale()
                NSA_2.scale()

                # Export kspaces
                save_pyobject([NSA_1,NSA_2],'inputs/kspaces/Ck_{:03d}_{:02d}_{:02d}.pkl'.format(d,fn,rn))


# DENSE images generation
def generate_dense(resolutions, frequencies, patients, ini=0, fin=0, noise_free=False):

    # Create folder
    if not os.path.isdir('inputs/kspaces'):
        os.mkdir('inputs/kspaces')

    # Resolutions loop
    for (rn, r) in enumerate(resolutions):

        # Create image
        I = DENSEImage(FOV=np.array([0.1, 0.1, 0.008]),
                  center=np.array([0.0,0.0,0.0]),
                  resolution=r,
                  encoding_frequency=np.array([0,0,0]),
                  T1=0.85,
                  flip_angle=15*np.pi/180,
                  kspace_factor=15,
                  slice_thickness=0.008,
                  oversampling_factor=1,
                  phase_profiles=r[1])

        # Filter specifications
        I.filter = 'Riesz'
        I.filter_width = 0.8
        I.filter_lift = 0.3

        # Frequencies loop
        for (fn, f) in enumerate(frequencies):

            # Set encoding frequency
            I.encoding_frequency = f
         
            # Data loop
            for d in range(ini,fin):

                # Load spins and parameter objects
                spins = load_pyobject('inputs/spins/spins_{:03d}.pkl'.format(d),sep_proc=True)
                param = load_pyobject('inputs/parameters/p_{:03d}.pkl'.format(d))

                # Create phantom
                phantom = Phantom(spins, param, patient=patients[d], z_motion=False)

                # Artifact
                artifact = None

                # Generate kspaces
                NSA_1, NSA_2, REF, mask = I.generate(artifact, phantom, param, debug=False)

                # Export noise-free images
                if noise_free:
                    # Create folder
                    if not os.path.isdir('inputs/noise_free_images'):
                        os.mkdir('inputs/noise_free_images')

                    # Get images and scale
                    I1 = scale_image(NSA_1.to_img(),mag=False,real=True,compl=True)
                    I2 = scale_image(NSA_2.to_img(),mag=False,real=True,compl=True)

                    # Plot test
                    if MPI_rank==0:
                      III1 = NSA_1.to_img()
                      III2 = NSA_2.to_img()
                      III = III1 - III2
                      multi_slice_viewer(np.abs(NSA_1.k_msk[...,0,0,:]))
                      multi_slice_viewer(np.abs(III[...,0,0,:]))

                    # Export images
                    save_pyobject([I1,I2],'inputs/noise_free_images/DI_{:03d}_{:02d}_{:02d}.pkl'.format(d,fn,rn))

                # Compress kspaces
                NSA_1.scale()
                NSA_2.scale()

                # Export kspaces
                save_pyobject([NSA_1,NSA_2],'inputs/kspaces/Dk_{:03d}_{:02d}_{:02d}.pkl'.format(d,fn,rn))
