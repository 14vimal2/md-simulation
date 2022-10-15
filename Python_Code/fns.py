import numpy as np
from particles import Particles
import para
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.spatial import cKDTree
rcParams["figure.figsize"] = (6,5)



def compute_interacting_pairs(P: Particles) -> None:
    """find interacting pairs by computing pairwise nearest image distance
     between all particles, updates (N x N) d2matrix and all 
     the k (= no. of interacting pairs) pairs ,i.e. iarr, jarr and r2arr 
     for which d2matrix is less than rc2"""
    # matrix based
    P.d2matrix = P.d2matrix*0
    for j in range(para.dim):
        temparr = np.abs(P.r[:,j].reshape(para.N,1) - P.r[:,j].reshape(1,para.N))
        boolarr = temparr <= para.L / 2
        temparr = temparr * boolarr + (para.L - temparr) * (~ boolarr)
        P.d2matrix += temparr**2
    P.iarr, P.jarr = np.where(P.d2matrix < para.rc**2)
    k = P.jarr > P.iarr
    P.iarr, P.jarr = P.iarr[k], P.jarr[k]
    P.r2arr = P.d2matrix[P.iarr, P.jarr]

    # tree based
    # tree = cKDTree(P.r, boxsize=para.L)
    # colliding_indices = tree.query_pairs(para.rc,output_type = 'ndarray')
    # P.iarr, P.jarr = colliding_indices[:,0], colliding_indices[:,1]
    # dist_ = np.abs(P.r[P.iarr]-P.r[P.jarr])
    # boolarr = dist_ < para.L / 2
    # dist_ = dist_ * boolarr + (para.L - dist_) * (~ boolarr)
    # P.r2arr = np.sum(dist_**2, axis=1)


    P.collisionnumber += len(P.r2arr)
    # if P.r2arr.any() == 0 and len(P.r2arr != 0):        
    #     print(np.where(P.r2arr == 0))
    #     print('\n')
    #     print('\n',P.iarr,'\n',P.jarr,'\n',P.r2arr)
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
    r6i = P.r2arr**(-3)
    P.PE = 4 * np.sum(r6i * (r6i - 1) ) - para.PE_cut * len(r6i)


def time_advance_single_step(P: Particles, dt: float) -> None:
    """Advances the particle by dt for given algorithm, e.g.- verlet algorithm"""
    rn, P.v = time_advance_verlet(P,dt)

    P.rp, P.r = P.r, rn

def time_advance_verlet(P: Particles, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """return Positions and velocity array using verlet formulas"""
    rdiff = P.r - P.rp 
    rn = 2*P.r-P.rp+2*para.L*(rdiff<-para.L/2)-2*para.L*(rdiff>para.L/2)+P.F*dt**2
    rn = rn  + (rn < 0) * para.L - (rn > para.L) * para.L
    rdiff = rn - P.rp
    rdiff = rdiff + para.L * (rdiff < -para.L/2) - para.L * (rdiff > para.L/2)
    vn = rdiff/(2*dt)
    return rn, vn


def compute_KE(P: Particles) -> None:
    """update Kinetic Energy of Particles"""
    P.KE = np.sum(P.v**2)/2

def compute_T(P: Particles) -> None:
    """update Temperatur of Particles"""
    P.T = np.sum(P.v**2)/(para.dim * para.N)


def showVelocityDist(P: Particles, bins=10, tag=0) -> None:
    velocities2 = np.sum(P.v**2,axis=1)**(1/2)
    # velocities2 = velocities2[np.where(velocities2 < 15)]
    np.save('vel_{}_{}_{}.npy'.format(para.N, para.T, tag),velocities2)
    _, ax = plt.subplots()

    a, b, *_ = ax.hist(velocities2,bins=bins,alpha=0.1,color='y' )
    b = (b[0:-1] + b[1:])/2
    ax.plot(b,a, label= 'speed')
    ax.set_xlim(np.min(b),np.max(b))
    plt.title('Speed distribution of particles at {}'.format(tag*para.dt))
    ax.set_xlabel('speed')
    ax.set_ylabel('count')
    plt.savefig('veldist_{}.png'.format(tag))


def showEnergyWithTime(xarr, KE_arr,  PE_arr):
    fig, ax = plt.subplots()
    plt.title('Energy vs Time of Particles')
    ax.plot(xarr, PE_arr, label=r'$PE$')
    ax.plot(xarr, KE_arr, label=r'$KE$')
    ax.set_xlabel(r'time ($t^*=t\sqrt{ \frac{\epsilon}{m \sigma^2}  }$)', fontsize=10)
    ax.set_ylabel(r'energy ($E^*= \frac{E}{\epsilon}$)', fontsize=10)
    ax.set_xlim(0,xarr[-1])
    ax.set_ylim(np.min(PE_arr)-30,np.max(KE_arr)+ 30)
    plt.margins(20,40)
    E_mean = np.average(KE_arr + PE_arr)
    print('average envergy',E_mean)
    std_dev_E = np.sqrt(np.average((KE_arr+PE_arr-E_mean)**2))
    print('standard dev energy',std_dev_E)
    ax.plot(xarr, PE_arr + KE_arr,'k--', label=r'$KE + PE$')
    plt.legend(loc=7)
    plt.show()
