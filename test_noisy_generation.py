import numpy as np

from utils.args_parser import *
from utils.im_parameters import (cfrequencies, dfrequencies, noise_levels,
                                 patients_arr, resolutions)
from utils.generation import cspamm_noisy_images, dense_noisy_images

##################
# DATA GENERATION
##################

# Add noise to CSPAMM images
cspamm_noisy_images(resolutions, cfrequencies, noise_levels, ini=args.initial_data)

# Add noise to DENSE images
dense_noisy_images(resolutions, dfrequencies, noise_levels, ini=args.final_data)
