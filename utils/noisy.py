import numpy as np
from PyMRStrain import *

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
                    I1 = scale_image(NSA_1.to_img(),mag=False,real=True,compl=True)
                    I2 = scale_image(NSA_2.to_img(),mag=False,real=True,compl=True)

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
                    I1 = scale_image(NSA_1.to_img(),mag=False,real=True,compl=True)
                    I2 = scale_image(NSA_2.to_img(),mag=False,real=True,compl=True)

                    # Export noisy kspaces
                    save_pyobject([I1,I2],'inputs/noisy_images/DI_{:03d}_{:02d}_{:02d}_{:02d}.pkl'.format(d,fn,rn,nn))