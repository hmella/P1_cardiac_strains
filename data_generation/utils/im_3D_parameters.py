import numpy as np

# Noise levels
noise_levels = np.array([1e-30, 3.535e-02, 4.350e-02, 5.500e-02, 7.725e-02])

# DENSE imaging parameters
keD = 0.12*(1000*2*np.pi)     # encoding frequency [rad/m]
resD = np.array([129,129,1])  # resolution

# SPAMM encoding frequency
# keC = 2*np.pi/0.016           # encoding frequency [rad/m]
keC = 2*np.pi/0.008           # encoding frequency [rad/m]
resC = np.array([257,257,1])  # resolution

# Shared imaging parameters
FOV = [np.array([0.35, 0.35, 0.030]),    # Base
        np.array([0.35, 0.35, 0.025]),   # Mid
        np.array([0.35, 0.35, 0.020])]   # Apex
center = [np.array([0.0, 0.0, 0.025]),   # Base
          np.array([0.0, 0.0, 0.000]),   # Mid
          np.array([0.0, 0.0, -0.025])]  # Apex
offset = [0.012, 0.006, 0.000]
slice_thickness = 0.008

# T1 and M0
T1 = np.array([1e-10,1e-10,0.85])
M0 = np.array([0,0,1])

# Filter specs
filter3       = 'Riesz'
filter_width3 = 0.8
filter_lift3  = 0.0

receiver_bw=64*1000
echo_train_length=9
off_resonance=115
spatial_shift='top-down'