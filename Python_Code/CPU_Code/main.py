import para
import numpy as np
import time
from matplotlib import rcParams
from fns import * #compute_force, compute_interacting_pairs, compute_PE, time_advance_single_step, compute_KE, compute_T
from particles import Particles
rcParams["figure.figsize"] = (6,5)
t1 = time.time()

P = Particles()
P.set_arrays()
P.init()
ti = 0
n = int(para.tf/para.dt)
PE_arr = np.zeros(n)
KE_arr = np.zeros(n)
print('size of the box', para.L)
print('dimension', para.dim)

# showVelocityDist(P,bins=10)
# times = np.arange(0,n,n//10)
# print('time\tTemperature\t \tKinetic Energy\t \tPotential Energy\t\tTotal Energy',)
while (ti < n):
    # must compute for running simulation
    compute_interacting_pairs(P)
    compute_force(P)

    # compute for sample average
    compute_PE(P)
    compute_KE(P)
    # compute_T(P)
    PE_arr[ti], KE_arr[ti] = P.PE, P.KE 
    time_advance_single_step(P,para.dt)
    ti += 1
    # if ti in times:
    #     showVelocityDist(P, tag=ti)
        # print(ti*para.dt,'\t',P.T,'\t',P.KE,'\t',P.PE,'\t', P.KE + P.PE)




print('total collision count',P.collisionnumber)
xarr = np.arange(0,para.tf, para.dt)


# showVelocityDist(P,bins=25)
# to see energy plot
t2 = time.time()
print('time take in simulation',t2 - t1)
showEnergyWithTime(xarr, KE_arr, PE_arr)
# plt.show()