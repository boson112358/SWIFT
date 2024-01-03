# coding: utf-8
import h5py

f = h5py.File('uniformBox_0001.hdf5','r+')

a_group_key = list(f.keys())[8]

ds_gradient_v=f[a_group_key][('Gradient_v')]
import matplotlib.pyplot as plt
ds_coo = f[a_group_key][('Coordinates')]
fig, axs = plt.subplots(ncols=3, nrows=3, layout='constrained', figsize=(10,10))
        
index = 0
for row in range(3):
    for col in range(3):
        print(index)
        axs[row, col].annotate(f'[{row}, {col}]', (0.5, 0.5),
                               transform=axs[row, col].transAxes,
                               ha='center', va='center', fontsize=18,
                               color='darkgrey')
        axs[row,col].plot(ds_coo[:,row], ds_gradient_v[:,index])
        axs[row,col].set_ylim([0, 2])
        index +=1
        
fig.suptitle('Gradient_v after one timestep[axis,grad_v]')
plt.savefig('gradient.png')
