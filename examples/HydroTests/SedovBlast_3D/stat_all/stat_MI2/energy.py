import numpy as np
a_mi2_full = np.loadtxt('statistics_MI2_full.txt')
kin1 = a_mi2_full[:,13]
inter_u1 = a_mi2_full[:,14]

a_mi2_noac = np.loadtxt('statistics_MI2_noac.txt')
kin2 = a_mi2_noac[:,13]
inter_u2 = a_mi2_noac[:,14]

a_mi2_noacav = np.loadtxt('statistics_MI2_noacav.txt')
kin3 = a_mi2_noacav[:,13]
inter_u3 = a_mi2_noacav[:,14]

a_sphenix = np.loadtxt('statistics_sphenix.txt')
kin4 = a_sphenix[:,13]
inter_u4 = a_sphenix[:,14]

a_minimal = np.loadtxt('statistics_minimal.txt')
kin5 = a_minimal[:,13]
inter_u5 = a_minimal[:,14]

a_mi2_noav = np.loadtxt('statistics_MI2_noav.txt')
kin6 = a_mi2_noav[:,13]
inter_u6 = a_mi2_noav[:,14]

a_mi2_matthieu = np.loadtxt('statistics_MI2_matthieu.txt')
kin7 = a_mi2_matthieu[:,13]
inter_u7 = a_mi2_matthieu[:,14]


import matplotlib.pyplot as plt
plt.figure()
time = a_mi2_full[:,1]

plt.plot(time, inter_u1 + kin1, label='MI2_full')
plt.plot(time, inter_u2 + kin2, label='MI2_noac')
plt.plot(time, inter_u3 + kin3, label='MI2_noacav')
plt.plot(time, inter_u6 + kin6, label='MI2_noav')
plt.plot(time, inter_u4 + kin4, linewidth=4.0, label='sphenix')
plt.plot(time, inter_u5 + kin5, label='minimal')
plt.plot(time, inter_u7 + kin7, label='MI2_matthieu')

plt.ylim(0.996,1.008)
plt.xlabel('time/s')
plt.ylabel('energy')
plt.legend()
plt.savefig('Energy.png')
