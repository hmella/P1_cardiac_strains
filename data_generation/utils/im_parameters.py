import numpy as np

# Resolutions
resolutions = [np.array([33, 33, 1]),    # pixel_size = 3.0  [mm]
               np.array([40, 40, 1]),    # pixel_size = 2.5  [mm]
               np.array([50, 50, 1]),    # pixel_size = 2.0  [mm]
               np.array([67, 67, 1]),    # pixel_size = 1.5  [mm]
               np.array([100, 100, 1])]  # pixel_size = 1.0  [mm]

# Encoding frequencies
tag_spacings = [0.0080, 0.0100, 0.0120, 0.0140, 0.0160]                   # tag spacings [m]
cfrequencies = [np.array([2*np.pi/l,2*np.pi/l,0]) for l in tag_spacings]  # encoding frequency [rad/m]
dfrequencies = [0.12*1000*2*np.pi*np.array([1,1,0])]                      # encoding frequency [rad/m]

# Noise levels
noise_levels = np.array([1e-30, 2e-02, 4e-02, 6e-02])

# Patients array
patients_arr = np.ones([100,],dtype=np.bool)
patients_arr[0:int(patients_arr.size/2)] = False

# Filter specs
filter       = None
filter_width = 0.8
filter_lift  = 0.3