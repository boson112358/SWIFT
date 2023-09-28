#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from glob import glob
import os

fig, ax = plt.subplots(1, 2, sharex=True, figsize=(10, 5))

addr_book = sorted(glob('*/*/*.txt'))

Npart = 32*32*32/1
h = 1.3 * Npart**(-1/3)

print(Npart,h)

def plot_statistics_file(the_addr, the_name):
    the_statistics=np.transpose(np.loadtxt(the_addr))
    
    Time = np.array(the_statistics[1])
    #E_kin = np.array(the_statistics[13])
    #E_int = np.array(the_statistics[14])
    E_mag = np.array(the_statistics[34])
    divB = np.array(the_statistics[35])
    #E_tot=E_kin+E_int+E_mag
    B = np.sqrt(2*E_mag)

    
    B0 = B[0]

    alpha = find_growth_rate(Time, np.log(B/B0), nlast = 10)[0]
    print(f'growth rate is {alpha}')
    largest_grad = B/h
    R = divB/divB[0]
    ax[0].plot(Time, B/B0, label=the_name)
    ax[1].plot(Time, R, label=the_name)

def find_growth_rate(the_time, B_field, nlast = 3):
    l = len(B_field)
    B_field_cut = B_field[l-1-nlast:-1]
    time_cut = the_time[l-1-nlast:-1]
    #print(time_cut,B_field_cut)
    res = np.polyfit(time_cut,B_field_cut,1)
    return res

for i in range(len(addr_book)):
    addr = addr_book[i]
    name = os.path.basename(os.path.dirname(addr))
    print(name)
    plot_statistics_file(addr,name)
#ax.plot(Time, E_int/Emag_eq,label="E_int/Emag_eq")
#ax.plot(Time, E_mag/Emag_eq,label="E_mag/Emag_eq")
#ax.plot(Time, one ,label="1")
ax[0].set_xlabel("Time [s]")
ax[1].set_xlabel("Time [s]")
ax[0].set_ylabel("<B_rms>/<B_rms(0)>")
ax[1].set_ylabel("divB/divB(0)")
ax[0].legend(loc="best")
ax[1].legend(loc="best")
ax[0].set_yscale("log")
ax[1].set_yscale("log")
ax[0].set_ylim([1e-1, 1e3])
#ax[1].plot(Time, E_tot/E_tot[0],label="E(t)")
#ax[1].plot(Time, E_tot/E_tot,label="Energy conservation",color="Black",linestyle="dashed")
#ax[1].set_xlabel("Time [s]")
#ax[1].set_ylabel("E_tot(t) / E_tot(0)")
#ax[1].legend(loc="best")
#plt.tight_layout()


#ax[1].plot(Time, E_mag/Emag_eq,label="Relative_magnetic_energy")
#ax[1].plot(Time, one ,color="Black",linestyle="dashed")
#ax[1].plot(Time, 0.001*np.exp(Time * Rm * v0),color="Blue",linestyle="dashed")
#ax[1].set_xlabel("Time [s]")
#ax[1].set_ylabel("E(t)/Eeq")
#ax[1].set_yscale("log")
#ax[1].legend(loc="best")
#plt.tight_layout()
plt.savefig("B_compare.png", dpi=100)
