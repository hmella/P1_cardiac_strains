import numpy as np

from utils.args_parser import *
from utils.generation import generate_cspamm, generate_dense, generate_phantoms
from utils.im_parameters import cfrequencies, dfrequencies, resolutions

##################
# INPUT PARAMETERS
##################

# Patients array
patients = np.ones([args.nb_data,],dtype=np.bool)
patients[0:int(patients.size/2)] = False 


##################
# DATA GENERATION
##################

# Generate phantoms
if args.generate_parameters:
    generate_phantoms(args.nb_samples,ini=args.initial_data,fin=args.final_data)

# Generate CSPAMM
if args.cspamm:
    generate_cspamm(resolutions, cfrequencies, patients, ini=args.initial_data,
                    fin=args.final_data, noise_free=args.noise_free)

# Generate DENSE
if args.dense:
    generate_dense(resolutions, dfrequencies, patients, ini=args.initial_data,
                  fin=args.final_data, noise_free=args.noise_free)
