import numpy as np

from utils.args_parser import *
from utils.generation import generate_cspamm, generate_dense, generate_phantoms

##################
# INPUT PARAMETERS
##################

# Number of data and spins samples
nb_data = args.nb_data
nb_samples = args.nb_samples

# Patients array
patients = np.ones([nb_data,],dtype=np.bool)
patients[0:int(patients.size/2)] = False 

# Resolutions
resolutions = [np.array([33, 33, 1]),    # pixel_size = 3.0  [mm]
               np.array([40, 40, 1]),    # pixel_size = 2.5  [mm]
               np.array([50, 50, 1]),    # pixel_size = 2.0  [mm]
               np.array([67, 67, 1]),    # pixel_size = 1.5  [mm]
               np.array([100, 100, 1])]  # pixel_size = 1.0  [mm]

# Encoding frequencies
tag_spacings = [0.0080, 0.0100, 0.0120, 0.0140, 0.0160] # tag spacings [m]
cfrequencies = [np.array([2*np.pi/l,2*np.pi/l,0]) for l in tag_spacings]        # encoding frequency [rad/m]
dfrequencies = [0.12*1000*2*np.pi*np.array([1,1,0])]                      # encoding frequency [rad/m]

##################
# DATA GENERATION
##################

# Generate phantoms
if args.generate_params:
    generate_phantoms(nb_samples,nb_data)

# Generate CSPAMM
generate_cspamm(resolutions, cfrequencies, patients, nb_data)

# Generate DENSE
generate_cspamm(resolutions, dfrequencies, patients, nb_data)
