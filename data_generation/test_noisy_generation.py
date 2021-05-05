import os

import numpy as np
from PyMRStrain.IO import load_pyobject, scale_image
from PyMRStrain.Math import itok, ktoi
from PyMRStrain.Noise import add_cpx_noise
from PyMRStrain.Plotter import multi_slice_viewer
from scipy.io import savemat

from utils.args_parser import *
from utils.im_parameters import (harp_spacings, noise_levels, patients_arr,
                                 resolutions, sinmod_spacings, tag_spacings)


# Input folders
folders = ['inputs/kspaces/']

# Output folders
output_folder = 'inputs/noisy_images/'
os.makedirs(output_folder,exist_ok=True)

# Number of data
nb_data = len(patients_arr)
nb_noise = len(noise_levels)
nb_res = len(resolutions)

# Add noise to synthetic MR images
for folder in folders:
    for d in range(args.initial_data, args.final_data):
        print('  Generating data {:d}'.format(d))

        for r in range(nb_res):

            print('[H] Pixel-spacing {:f}, tag-spacing {:f}'.format(0.1/resolutions[r][0],tag_spacings[harp_spacings[r]]))
            print('[S] Pixel-spacing {:f}, tag-spacing {:f}'.format(0.1/resolutions[r][0],tag_spacings[sinmod_spacings[r]]))
            # Filenames
            H_filename = 'CI_{:03d}_{:02d}_{:02d}.pkl'.format(d,harp_spacings[r],r)
            S_filename = 'CI_{:03d}_{:02d}_{:02d}.pkl'.format(d,sinmod_spacings[r],r)
            D_filename = 'DI_{:03d}_{:02d}_{:02d}.pkl'.format(d,0,r)

            # Load images
            [H1,H2] = load_pyobject(folder+H_filename)
            [S1,S2] = load_pyobject(folder+S_filename)
            [D1,D2] = load_pyobject(folder+D_filename)

            # Rescale
            H1.rescale(), H2.rescale()
            S1.rescale(), S2.rescale()
            D1.rescale(), D2.rescale()

            # Maximum intensity of the image and standard deviation
            # of the noise
            if r == 0:
                I_max_r0 = max([np.abs(ktoi(D1.k)[...,0,0]).max(), np.abs(ktoi(D1.k)[...,1,0]).max()])

            for n, nlevel in enumerate(noise_levels):

                # Standard deviation of the noise
                std = nlevel*I_max_r0

                # Generate noise in the image domain
                noise_1 = np.random.normal(0, std, D1.k.shape) \
                        + 1j*np.random.normal(0, std, D1.k.shape)
                noise_2 = np.random.normal(0, std, D2.k.shape) \
                        + 1j*np.random.normal(0, std, D2.k.shape)

                # Noise in the fourier domain
                f_noise_1 = D1.k_msk*itok(noise_1)
                f_noise_2 = D2.k_msk*itok(noise_2)

                # Noisy kspaces
                nH1, nH2 = (H1.k + f_noise_1), (H2.k + f_noise_2)
                nS1, nS2 = (S1.k + f_noise_1), (S2.k + f_noise_2)
                nD1, nD2 = (D1.k + f_noise_1), (D2.k + f_noise_2)

                # import matplotlib.pyplot as plt
                # fig, ax = plt.subplots(1,2)
                # ax[0].imshow(np.abs(nD1[...,0,0,9]))
                # ax[1].imshow(np.abs(nD1[...,0,0,9] - nD2[...,0,0,9]))
                # plt.show()

                # kspace to images
                IH1, IH2 = ktoi(nH1), ktoi(nH2)
                IS1, IS2 = ktoi(nS1), ktoi(nS2)
                ID1, ID2 = ktoi(nD1), ktoi(nD2)

                # Scale images
                IH = scale_image(IH1-IH2, mag=False,real=True,compl=True)
                IS = scale_image(IS1-IS2, mag=False,real=True,compl=True)
                ID = scale_image(ID1-ID2, mag=False,real=True,compl=True)

                # Output filenames
                H_filename = 'HI_{:03.0f}_{:02.0f}_{:02.0f}.pkl'.format(d,n,r)
                S_filename = 'SI_{:03.0f}_{:02.0f}_{:02.0f}.pkl'.format(d,n,r)
                D_filename = 'DI_{:03.0f}_{:02.0f}_{:02.0f}.pkl'.format(d,n,r)

                # Export data
                savemat(output_folder+H_filename[:-4]+'.mat',{'I': IH})
                savemat(output_folder+S_filename[:-4]+'.mat',{'I': IS})
                savemat(output_folder+D_filename[:-4]+'.mat',{'I': ID})


# Generate EXACT images
for d in range(args.initial_data, args.final_data):
    print('  Generating data {:d}'.format(d))
    for r in range(nb_res):

        # Filename
        filename = 'EI_{:03d}_{:02d}_{:02d}.pkl'.format(d,0,r)

        # Load image
        K = load_pyobject(folder+filename)

        # Rescale
        K.rescale()

        # Inverse FFT
        I = ktoi(K.k)

        # Scale image
        I = scale_image(I, mag=False,real=True,compl=True)

        # Export
        savemat(output_folder+filename[:-4]+'.mat',{'I': I})
