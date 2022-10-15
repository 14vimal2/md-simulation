import para
import numpy as np
import time
from matplotlib import rcParams
from fns import *
from particles import Particles

rcParams["figure.figsize"] = (6,5)
para.device.use()
t1 = time.time()


P = Particles()
P.set_arrays()
P.init()
# coarseGraining(P,10,True,'Velocity_fields_before.png')
# addVortices(P,para.eta,True, para.nvorticesrt)
# coarseGraining(P,10,True,'Velocity_fields_after.png')
ti = 0
n = int(para.tf/para.dt)
# PE_arr = np.zeros(n)
# KE_arr = np.zeros(n)
print('size of the box', para.L)
print('dimension', para.dim)

# times = np.arange(n-1000,n,50)

while (ti < n):
    # must compute for running simulation
    compute_interacting_pairs(P)
    compute_force(P)

    # compute for sample average
    # compute_PE(P)
    # compute_KE(P)
    # compute_T(P)
    # PE_arr[ti], KE_arr[ti] = P.PE, P.KE
    time_advance_single_step(P,para.dt)
    ti += 1

print('total collision count',P.collisionnumber)
t2 = time.time()
print('time take in simulation', t2 - t1)
# coarseGraining(P,10,True,'Velocity_fields_after_simulation.png')
np.save('mb_positions_{}_dim{}_.npy'.format(para.N, para.dim),P.r)
np.save('mb_velocity_{}_dim{}_.npy'.format(para.N, para.dim),P.v)
# np.save('PE_arr_{}.npy'.format(para.N),PE_arr)
