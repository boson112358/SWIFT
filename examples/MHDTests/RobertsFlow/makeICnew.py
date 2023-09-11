#############################

import h5py
import numpy as np
import matplotlib.pyplot as plt

# Parameters

rho = 1.0
cs2 = 3025.0 #3025.0
L = 1.0
k = 2 * np.pi / L
V0 = 0.12 #16 * np.pi

gamma = 5.0 / 3.0
u0 = cs2 / (gamma * (gamma - 1))

fileOutputName = "RobertsFlow.hdf5"

###---------------------------###

glass = h5py.File("glassCube_32.hdf5", "r")

pos = glass["/PartType0/Coordinates"][:, :]
h = glass["/PartType0/SmoothingLength"][:]

N = len(h)
vol = L ** 3

###---------------------------###

v = np.zeros((N, 3))
B = np.zeros((N, 3))
psi = np.zeros(N)
ids = np.linspace(1, N, N)
m = np.ones(N) * rho * vol / N
u = np.ones(N) * u0

psi[:] = (V0 / k) * np.cos(k * pos[:,0]) * np.cos(k * pos[:,1]) 

v[:, 0] = V0 * np.cos(k * pos[:, 0]) * np.sin(k * pos[:, 1]) # - was here
v[:, 1] = - V0 * np.sin(k * pos[:, 0]) * np.cos(k * pos[:, 1]) # + was here
v[:, 2] = k * np.sqrt(2) * psi[:]

#B[:, 0] = k * k * np.sin(k * pos[:, 1]) * np.sin(k * pos[:, 2])
B[:,0] = k * k * np.cos(k * pos[:, 2])
B[:,1] = k * k * np.sin(k * pos[:, 2])

###---------------------------###

# File
fileOutput = h5py.File(fileOutputName, "w")

# Header
grp = fileOutput.create_group("/Header")
