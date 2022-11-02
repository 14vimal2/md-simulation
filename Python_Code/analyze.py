from post_proc import *
import h5py
from scipy.spatial import cKDTree
i = 0

data_part = 21
t = 0

# f = h5py.File("%s/data_particles_%d_ensemble_%d_part_%d.hdf5"%(para.dir_name, para.N,i,data_part), "r")
# f = h5py.File("%s/data_particles_%d_ensemble_%d_part_%d.hdf5"%('simulation_data_test_1', para.N,i,data_part), "r")
# f = h5py.File('vortices_3/data_particles_1000_ensemble_1_part_%d.hdf5'%(data_part), "r")
f = h5py.File('data_particles_1000_ensemble_0_part_%d.hdf5'%(data_part), "r")

# print(f.keys())
pos = f['particles_positions_with_time']
vel = f['particles_velocities_with_time']
# print(val[:,:,0])


# showVelDist(vel[:,:,t])
# displayVelField(pos[:,:,t],vel[:,:,t])
# coarseGraining(pos[:,:,t],vel[:,:,t],10)
# print(pos[:,:,t])

tree = cKDTree(pos[:,:,t],boxsize=para.L)

# xv, yv = np.meshgrid(np.linspace(0,para.L,100),np.linspace(0,para.L,100))
# print(xv.shape)

no_of_gridpoints_1d = 30
def makeGridPoints(no_of_gridpoints_1d:int):
    rprime = np.empty((no_of_gridpoints_1d**2,para.dim),dtype=np.float32)
    index_ = np.arange(no_of_gridpoints_1d**2)
    for j in range(para.dim):
        rprime[:,j] = (index_//no_of_gridpoints_1d**j) % no_of_gridpoints_1d
    return ((rprime+0.5)/no_of_gridpoints_1d) * para.L

def gridPointsVelocities(rprime: np.ndarray, tree: cKDTree,v: np.ndarray, no_of_gridpoints_1d = no_of_gridpoints_1d, radius= 0.23):
    radius = radius * para.L
    N_ = no_of_gridpoints_1d**2
    vprime = np.empty((N_,para.dim),dtype=np.float32)
    for i in range(N_):
        particles_indices = tree.query_ball_point(rprime[i], radius)
        # print(particles_indices)
        vprime[i] = np.average(v[particles_indices],axis=0)
    return vprime

r_ = makeGridPoints(no_of_gridpoints_1d)
v_ = gridPointsVelocities(r_,tree,vel[:,:,t])
def getVorticesVel(no_of_vortices_sqrt, r):
    endpoint = no_of_vortices_sqrt*np.pi
    u = np.zeros(r.shape)
    x, y = r[:,0]/para.L * endpoint, r[:,1]/para.L * endpoint
    u[:,0] = np.sin(x)*np.cos(y)
    u[:,1] = -np.cos(x) * np.sin(y) 
    return u

# displayVelField(r_, getVorticesVel(3,r_))
displayVelField(r_,v_,show=True)
# particle_id = 7
# dd, ll = tree.query(pos[:,:,t][particle_id],k=10,p=2)
# coordinates = np.array([[250, 200]])
# dd, ll = tree.query(coordinates,k=20)
# print(dd, ll)
# plt.plot(pos[particle_id,0,t], pos[particle_id,1,t], 'bs')
# plt.plot(pos[:,0,t], pos[:,1,t], 'k.', lw=0.5)
# plt.plot(pos[:,0,t][ll], pos[:,1,t][ll], 'r.')
# plt.plot(coordinates[:,0], coordinates[:,1], 'bs')

# ll = tree.query_ball_point(coordinates, 50)[0]
# print(ll)
# plt.plot(pos[:,0,t], pos[:,1,t], 'k.', lw=0.5)
# plt.plot(pos[:,0,t][ll], pos[:,1,t][ll], 'r.')
# plt.show()
# plt.close()
