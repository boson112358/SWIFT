#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

the_statistics=np.transpose(np.loadtxt("statistics.txt"))

Time = the_statistics[1]
E_kin = the_statistics[13]
E_int = the_statistics[14]
E_mag = the_statistics[34]
E_tot=E_kin+E_int+E_mag

fig, ax = plt.subplots(1, 3, sharex=True, figsize=(10, 5))
ax[0].plot(Time, E_kin/E_tot,label="E_kin")
ax[0].plot(Time, E_int/E_tot,label="E_int")
ax[0].plot(Time, E_mag/E_tot,label="E_mag")
ax[0].plot(Time, E_tot/E_tot,label="E_tot",color="Black",linestyle="dashed")
ax[0].set_xlabel("Time [s]")
ax[0].set_ylabel("E/ E_tot")
ax[0].legend(loc="best")

ax[1].plot(Time, E_tot/E_tot[0],label="E(t)")
ax[1].plot(Time, E_tot/E_tot,label="Energy conservation",color="Black",linestyle="dashed")
ax[1].set_xlabel("Time [s]")
ax[1].set_ylabel("E_tot(t) / E_tot(0)")
ax[1].legend(loc="best")
plt.tight_layout()

ax[2].plot(Time, E_mag,label="Emag")
ax[2].plot(Time, E_mag,label="Magnetic energy",color="Black",linestyle="dashed")
ax[2].set_xlabel("Time [s]")
ax[2].set_ylabel("E_mag(t)")
ax[2].set_yscale("log")
ax[2].legend(loc="best")
plt.tight_layout()
plt.savefig("E_change.png", dpi=100)
