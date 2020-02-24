import numpy as np

from utils.args_parser import *
from utils.im_parameters import (cfrequencies, dfrequencies, noise_levels,
                                 resolutions)
from utils.noisy import cspamm_noisy_images, dense_noisy_images

##################
# INPUT PARAMETERS
##################

# Number of data and spins samples
nb_data = args.nb_data
nb_samples = args.nb_samples

# Patients array
patients = np.ones([nb_data,],dtype=np.bool)
patients[0:int(patients.size/2)] = False 


##################
# DATA GENERATION
##################

# Add noise to CSPAMM images
cspamm_noisy_images(resolutions, cfrequencies, noise_levels, nb_data)

# Add noise to DENSE images
dense_noisy_images(resolutions, dfrequencies, noise_levels, nb_data)