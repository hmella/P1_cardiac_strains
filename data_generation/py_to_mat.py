import os

import numpy as np
from PyMRStrain.IO import load_pyobject, scale_image
from PyMRStrain.Math import itok, ktoi
from PyMRStrain.Plotter import multi_slice_viewer
from scipy.io import savemat

def if_exact(filename):
    return filename[0:2] == 'EI'

# Input folders
folders = ['inputs/kspaces/','inputs/masks/']

# Output folders
os.makedirs('inputs/noise_free_images',exist_ok=True)

# Convert files
for folder in folders:
    for i, filename in enumerate(os.listdir(folder)):
        fname, ext = os.path.splitext(filename)
        if ext != '.mat':

            if np.mod(i, 50) == 0: print(i)

            if folder is not 'inputs/masks/':         

                # Load kspaces and rescale
                if if_exact(fname):
                    K1 = load_pyobject(folder+filename)
                    K1.rescale()
                    I = ktoi(K1.k)
                else:
                    [K1,K2] = load_pyobject(folder+filename)
                    K1.rescale()
                    K2.rescale()
                    I = ktoi(K1.k - K2.k)

                # Scale image
                I = scale_image(I,mag=False,real=True,compl=True)

                # Export matlab object
                savemat('inputs/noise_free_images/'+fname+'.mat',{'I': I})

            else:

                # Load mask
                I = load_pyobject(folder+filename)

                savemat(folder+fname+'.mat',{'M': I})
