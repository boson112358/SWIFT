# coding: utf-8
import pandas as pd
import h5py

f1 = h5py.File('uniformBox.hdf5','r+')

a_group_key = list(f1.keys())[1]
ds_vel = f1[a_group_key][('Velocities')]
ds_coo = f1[a_group_key][('Coordinates')]


import matplotlib.pyplot as plt
plt.plot(ds_coo[:,1], ds_vel[:,1])

plt.xlabel('Y')
plt.ylabel('v_y')
plt.savefig('velocity.png')
