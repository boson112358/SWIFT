#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

the_statistics=np.transpose(np.loadtxt("statistics.txt"))

Time = np.array(the_statistics[1])
E_kin = np.array(the_statistics[13])
E_int = np.array(the_statistics[14])
E_mag = np.array(the_statistics[34])
E_tot=E_kin+E_int+E_mag


Rm=6
eta=0.01
v0=Rm*2*np.pi*eta
rho0=1
Emag_eq=1*rho0*(v0)**2/2
one=np.ones(len(E_mag))

tmax=16
mask=Time<=tmax
Time=Time[mask]
E_kin=E_kin[mask]
E_int=E_int[mask]
E_mag=E_mag[mask]
one=one[mask]

B = np.sqrt(2*E_mag)
Beq = np.sqrt(2*Emag_eq)

def find_growth_rate(the_time, B_field, nlast = 3):
    l = len(B_field)
    B_field_cut = B_field[l-1-nlast:-1]
    time_cut = the_time[l-1-nlast:-1]
    #print(time_cut,B_field_cut)
    res = np.polyfit(time_cut,B_field_cut,1)
    return res

B0 = B[0]

alpha = find_growth_rate(Time, np.log10(B/B0), nlast = 3)[0]
print(f'growth rate is {alpha}')
#print(np.max(E_kin)/Emag_eq)
fig, ax = plt.subplots(1, 1, sharex=True, figsize=(5, 5))
ax.plot(Time, B/B0,label="<B_rms>/B_eq")
#ax.plot(Time, E_int/Emag_eq,label="E_int/Emag_eq")
#ax.plot(Time, E_mag/Emag_eq,label="E_mag/Emag_eq")
#ax.plot(Time, one ,label="1")
ax.set_xlabel("Time [s]")
ax.set_ylabel("<B_rms>/B_eq")
ax.legend(loc="best")
ax.set_yscale("log")
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
plt.savefig("E_change.png", dpi=100)
