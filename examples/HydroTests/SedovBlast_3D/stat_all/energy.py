import numpy as np
a_mi1_full = np.loadtxt('statistics_MI1_full.txt')
kin1 = a_mi1_full[:,13]
inter_u1 = a_mi1_full[:,14]

a_mi1_noac = np.loadtxt('statistics_MI1_noac.txt')
kin2 = a_mi1_noac[:,13]
inter_u2 = a_mi1_noac[:,14]

a_mi1_noacav = np.loadtxt('statistics_MI1_noacav.txt')
kin3 = a_mi1_noacav[:,13]
inter_u3 = a_mi1_noacav[:,14]

a_sphenix = np.loadtxt('statistics_sphenix.txt')
kin4 = a_sphenix[:,13]
inter_u4 = a_sphenix[:,14]

a_minimal = np.loadtxt('statistics_minimal.txt')
kin5 = a_minimal[:,13]
inter_u5 = a_minimal[:,14]

a_mi1_noav = np.loadtxt('statistics_MI1_noav.txt')
kin6 = a_mi1_noav[:,13]
inter_u6 = a_mi1_noav[:,14]


a_mi1_matthieu = np.loadtxt('statistics_MI1_matthieu.txt')
kin7 = a_mi1_matthieu[:,13]
inter_u7 = a_mi1_matthieu[:,14]


import matplotlib.pyplot as plt
plt.figure()
time = a_mi1_full[:,1]

plt.plot(time, inter_u1 + kin1, label='MI1_full')
plt.plot(time, inter_u2 + kin2, label='MI1_noac')
plt.plot(time, inter_u3 + kin3, label='MI1_noacav')
plt.plot(time, inter_u6 + kin6, label='MI1_noav')
plt.plot(time, inter_u4 + kin4, linewidth=4.0, label='sphenix')
plt.plot(time, inter_u5 + kin5, label='minimal')
plt.plot(time, inter_u7 + kin7, label='MI1_matthieu')
#plt.plot(time, inter_u8 + kin8, label='MI1_matthieu_noacav')

plt.ylim(0.996, 1.008)
plt.xlabel('time/s')
plt.ylabel('energy')
plt.legend()
plt.savefig('Energy.png')
