import numpy as np
import matplotlib.pyplot as plt
import time
import os
import h5py
t=np.zeros((25,50))
for l in np.arange(0,50):
    os.system('./test')
    f = h5py.File("time.h5", "r")
    t[:,l]=f['fields/t']
np.save("timedata",t)
