import os

from PyMRStrain import *


# Parameters generation
def generate_phantoms(nb_samples,nb_data):

    # Create folders
    if not os.path.isdir('inputs/'):
        os.mkdir('inputs/')
    if not os.path.isdir('inputs/parameters'):
        os.mkdir('inputs/parameters')
    if not os.path.isdir('inputs/spins'):
        os.mkdir('inputs/spins')

    for d in range(nb_data):

        # Create parameters
        p = Parameters(time_steps=18)

        # Create spins
        s = Spins(Nb_samples=nb_samples,parameters=p)

        # Write parameters and spins
        save_pyobject(p,'inputs/parameters/p_{:03d}.pkl'.format(d))
        save_pyobject(s,'inputs/spins/spins_{:03d}.pkl'.format(d), sep_proc=True)


# CSPAMM images generation
def generate_cspamm(resolutions, frequencies, patients, nb_data):

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

        # Frequencies loop
        for (fn, f) in enumerate(frequencies):

            # Set encoding frequency
            I.encoding_frequency = f
         
            # Data loop
            for d in range(nb_data):

                # Load spins and parameter objects
                spins = load_pyobject('inputs/spins/spins_{:03d}.pkl'.format(d),sep_proc=True)
                param = load_pyobject('inputs/parameters/p_{:03d}.pkl'.format(d))

                # Create phantom
                phantom = Phantom(spins, param, patient=patients[d], z_motion=False)

                # Artifact
                artifact = None

                # Generate kspaces
                NSA_1, NSA_2, REF, mask = I.generate(artifact, phantom, param, debug=False)

                # Compress kspaces
                NSA_1.scale()
                NSA_2.scale()

                # Export kspaces
                save_pyobject([NSA_1,NSA_2],'inputs/kspaces/Ck_{:03d}_{:02d}_{:02d}.pkl'.format(d,fn,rn))


# DENSE images generation
def generate_dense(resolutions, frequencies, patients, nb_data):

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

        # Frequencies loop
        for (fn, f) in enumerate(frequencies):

            # Set encoding frequency
            I.encoding_frequency = f
         
            # Data loop
            for d in range(nb_data):

                # Load spins and parameter objects
                spins = load_pyobject('inputs/spins/spins_{:03d}.pkl'.format(d),sep_proc=True)
                param = load_pyobject('inputs/parameters/p_{:03d}.pkl'.format(d))

                # Create phantom
                phantom = Phantom(spins, param, patient=patients[d], z_motion=False)

                # Artifact
                artifact = None

                # Generate kspaces
                NSA_1, NSA_2, REF, mask = I.generate(artifact, phantom, param, debug=False)

                # Compress kspaces
                NSA_1.scale()
                NSA_2.scale()

                # Export kspaces
                save_pyobject([NSA_1,NSA_2],'inputs/kspaces/Dk0_{:03d}_{:02d}_{:02d}.pkl'.format(d,fn,rn))
