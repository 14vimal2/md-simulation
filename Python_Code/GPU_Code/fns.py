import numpy as np
from particles import Particles
import para
# import matplotlib.pyplot as plt
import cupy as cp

para.device.use()

def compute_interacting_pairs(P: Particles) -> None:
    """find interacting pairs by computing pairwise nearest image distance between all particles, 
    updates (N x N) d2matrix and all the k (= no. of interacting pairs) pairs ,i.e. iarr, jarr and r2arr for which d2matrix is less than rc2"""
    P.d2matrix = P.d2matrix*0
    for j in range(para.dim):
        temparr = cp.abs(P.r[:,j].reshape(para.N,1) - P.r[:,j].reshape(1,para.N))
        boolarr = temparr <= para.L / 2
        temparr = temparr * boolarr + (para.L - temparr) * (~ boolarr)
        P.d2matrix += temparr**2
    P.iarr, P.jarr = cp.where(P.d2matrix < para.rc**2)
    k = P.jarr > P.iarr
    P.iarr, P.jarr = P.iarr[k], P.jarr[k]
    # print(P.d2matrix)
    # print(P.iarr, P.jarr)
    # print('no. of interacting pairs',len(P.iarr))
    P.r2arr = P.d2matrix[P.iarr, P.jarr]
    P.collisionnumber += len(P.r2arr)
    # print(P.r2arr)
    # if P.r2arr.any() == 0 and len(P.r2arr != 0):        
    #     # print(np.where())
    #     # print(P.iarr,P.jarr,P.r2arr)
    #     raise ValueError('distance became 0')

def compute_force(P: Particles) -> None:
    """Compute force between all interacting pairs"""
    P.F = 0 * P.F
    rij = P.r[P.iarr] - P.r[P.jarr]
    rij = rij + para.L * (rij < -para.L/2) - para.L * (rij > para.L/2)
    r2i = P.r2arr**(-1)
    r6i = r2i**3
    ff = 48*r2i*r6i*(r6i - 1/2)
    for l in range(len(r6i)):
        P.F[P.iarr[l]] += rij[l] * ff[l]
        P.F[P.jarr[l]] -= rij[l] * ff[l]

def compute_PE(P: Particles) -> None:
    """Compute Potential energy of particles"""
    r6i = cp.asnumpy( P.r2arr**(-3))
    # print('r6i',r6i)
    P.PE = 4 * np.sum(r6i * (r6i - 1) ) - para.PE_cut * len(r6i)


def time_advance_single_step(P: Particles, dt: float) -> None:
    """Advances the particle by dt for given algorithm, e.g.- verlet algorithm"""
    rn, P.v = time_advance_verlet(P,dt)
    P.rp, P.r = P.r, rn

def time_advance_verlet(P: Particles, dt: float):
    """return Positions and velocity array using verlet formulas"""
    rdiff = P.r - P.rp 
    rn = 2 * P.r - P.rp + 2 * para.L * (rdiff < -para.L/2) - 2 * para.L * (rdiff > para.L/2) + P.F * dt**2
    rn = rn  + (rn < 0) * para.L - (rn > para.L) * para.L
    rdiff = rn - P.rp
    rdiff = rdiff + para.L * (rdiff < -para.L/2) - para.L * (rdiff > para.L/2)
    vn = rdiff/(2*dt)
    return rn, vn


def compute_KE(P: Particles) -> None:
    """update Kinetic Energy of Particles"""
    P.KE = np.sum(cp.asnumpy(cp.sum(P.v**2,axis=1)))/2

def compute_T(P: Particles) -> None:
    """update Temperatur of Particles"""
    P.T = np.sum(cp.asnumpy(cp.sum(P.v**2,axis=1)))/(para.dim * para.N)



def addVortices(P: Particles, eta: float, anticlockwise: bool = True, nvorticesrt: int = 1):
    endpoint = nvorticesrt*cp.pi
    x, y = P.r[:,0]/para.L * endpoint, P.r[:,1]/para.L * endpoint
    uvel = anticlockwise * cp.sin(x) * cp.cos(y) * eta
    vvel = -anticlockwise * cp.cos(x) * cp.sin(y) * eta
    P.v[:,0] += uvel
    P.v[:,1] += vvel
    P.rp = P.r - P.v * para.dt

