from post_proc import *
import h5py
from scipy.spatial import cKDTree
import random
import numpy as np

import para
# i = 0
# data_part = 0
# f = h5py.File("%s/data_particles_%d_ensemble_%d_part_%d.hdf5"%(para.dir_name, para.N,i,data_part), "r")
# vel = f['particles_velocities_with_time']
# print(vel.shape)
# showVelDist(vel[:,:,-1], filename= 'veldist.png')


def makeGridPoints(no_of_gridpoints_1d:int):
    rprime = np.empty((no_of_gridpoints_1d**2,para.dim),dtype=np.float32)
    index_ = np.arange(no_of_gridpoints_1d**2)
    for j in range(para.dim):
        rprime[:,j] = (index_//no_of_gridpoints_1d**j) % no_of_gridpoints_1d
    return ((rprime+0.5)/no_of_gridpoints_1d) * para.L
no_of_gridpoints_1d = 30

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





def ProcessEnsembleParts(dirname, ensemble, part, r_=r_):
    f = h5py.File("%s/data_particles_%d_ensemble_%d_part_%d.hdf5"%(dirname, para.N,ensemble,part),'r+')
    pos = f['particles_positions_with_time']
    vel = f['particles_velocities_with_time']
    no_of_snapshot_in_part = pos.shape[-1]
    vs = np.empty((no_of_gridpoints_1d**2, para.dim, no_of_snapshot_in_part))
    for i in range(no_of_snapshot_in_part):
        tree = cKDTree(pos[:,:,i], boxsize=para.L)
        vs[:,:,i] = gridPointsVelocities(r_,tree,vel[:,:,i])
    f.create_dataset('coarse_grained_positions', data=r_)
    f.create_dataset('coarse_grained_velocities_with_time',data=vs)
    # print(vs)
    
        # displayVelField(r_,v_,'ensemble_%d_%d.png' % (ensemble, part * no_of_snapshot_in_part + i),dirname=dirname)


dirs = ['simulation_data_no_v{}_new'.format(i**2) for i in range(1,5)]
dir_i = 3
for ensemble in range(para.no_of_ensemble):
    for part in range(200):
        ProcessEnsembleParts(dirs[dir_i], ensemble,part)
exit()
particle_id = 350
no_of_parts = 200

velocities_with_time = np.zeros((2, no_of_parts*100))
positions_with_time = np.zeros((2, no_of_parts*100))
# print(velocities_with_time)
dirs = ['simulation_data_not_v_new','simulation_data_with_v_new']
# for dir in dirs:
#     for part in range(no_of_parts):
#         f = h5py.File('%s/data_particles_1000_ensemble_0_part_%d.hdf5' % (dir,part), 'r')
#         pos = f['particles_positions_with_time'][particle_id,:,:]
#         vel = f['particles_velocities_with_time'][particle_id,:,:]
#         positions_with_time[:,part*100:(part+1)*100] = pos
#         velocities_with_time[:,part*100:(part+1)*100] = vel
#         f.close()
#     np.save('%s_vel.npy' % dir , velocities_with_time)
#     np.save('%s_pos.npy' % dir , positions_with_time)


# for dir in dirs:
#     for part in range(1):
#         f = h5py.File('%s/data_particles_1000_ensemble_0_part_%d.hdf5' % (dir,part), 'r')
#         pos = f['particles_positions_with_time'][:,:,0] #[particle_id,:,:]
#         vel = f['particles_velocities_with_time'][:,:,0] #[particle_id,:,:]
#         tree = cKDTree(pos, boxsize=para.L)
#         v_ = gridPointsVelocities(r_,tree,vel)
#         displayVelField(r_, v_, figname=None, show=True)

no_of_particles_in_plot = 100

particle_ids = random.sample(range(1000), no_of_particles_in_plot) #(0,1000,no_of_particles_in_plot) np.random.randint(0,1000,no_of_particles_in_plot)
particle_ids.sort()
types_ = ['no vortices', 'with vortices']
print('particle ids', particle_ids)
fig = plt.figure(figsize=(8, 4))
velocities_with_time_ = np.zeros((len(particle_ids), 2, no_of_parts*100))
freq_ = np.fft.fftfreq(velocities_with_time_.shape[-1], 1/200)
idx_ = np.argsort(freq_)
# lim_ = 0.2# np.max(freq_)
for dir_i in range(len(dirs)):
    # ax = fig.add_subplot(1,2,dir_i+1)
    # ax.set_title(types_[dir_i])
    for part in range(no_of_parts):
        f = h5py.File('%s/data_particles_1000_ensemble_0_part_%d.hdf5' % (dirs[dir_i],part), 'r')
        # if part == 0:
        #     vel_ = f['particles_velocities_with_time'][:,:,-1]
        vel = f['particles_velocities_with_time'][particle_ids,:,:]
        velocities_with_time_[:, :,part*100:(part+1)*100] = vel
    no_of_data = 20000
    plotdata = np.zeros(10001)
    for i in range(len(particle_ids)):
        # newdata =  np.sqrt(velocities_with_time_[i,0,-no_of_data:]**2 + velocities_with_time_[i,1,-no_of_data:]**2)
        velocities_x, velocities_y = velocities_with_time_[i,0,:], velocities_with_time_[i,1,:]
        ftx, fty = np.fft.rfft(velocities_x), np.fft.rfft(velocities_y)
        plotdata += np.abs(ftx) + np.abs(fty)

    plt.loglog(plotdata, label=dirs[dir_i])
plt.xlabel('Time', fontdict={'size':15})
plt.ylabel('$U_x(t)$', fontdict={'size':15})
# plt.title(dirs[dir_i], fontdict={'size':20})
plt.tight_layout()
plt.show()
plt.close()
        # ax.plot(np.abs(ft[1:]), label='$p_i = %d$'% particle_ids[i])
        # print('xft\n',x_ft[1:])
        # print('yft\n',y_ft[1:])
    # ax.set_xlim(-lim_,lim_)
    # ax.set_xlabel('k')
    # ax.set_ylabel('$|u_x|$')
    # ax.set_ylim(bottom=0)
    # ax.legend()


exit()
# plt.tight_layout()
# plt.show()
# plt.close()

















if True:
    raise ValueError('program halted')

vel_abs = np.sum(vel_**2, axis=1)
ft = np.fft.fft(vel_abs)
freq = np.fft.fftfreq(ft.shape[-1])
idx = np.argsort(freq)
plt.plot(freq[idx], ft[idx])
plt.xlim(-0.1, 0.1)
plt.show()
plt.close()




vel_not_vortex = np.load('%s_vel.npy' % dirs[0], 'r').transpose()
vel_with_vortex = np.load('%s_vel.npy' % dirs[1], 'r').transpose()
pos_not_vortex = np.load('%s_pos.npy' % dirs[0], 'r').transpose()
pos_with_vortex = np.load('%s_pos.npy' % dirs[1], 'r').transpose()

print(vel_not_vortex.shape)

# tree = cKDTree(pos_not_vortex, boxsize=para.L)
tree = cKDTree(pos_not_vortex, boxsize=para.L)

x = np.linspace(0, para.L, 100)
X, Y = np.meshgrid(x, x)
Zx, Zy = np.zeros((100,100)), np.zeros((100,100))

fig = plt.figure()

ax = fig.add_subplot(2,2,1)
for i in range(100):
    for j in range(100):
        particles_indices = tree.query_ball_point([X[i, j], Y[i, j]], 3)
        Zx[i,j], Zy[i,j] = np.sum(vel_not_vortex[particles_indices],axis=0)[0], np.sum(vel_not_vortex[particles_indices],axis=0)[1]
ft = np.fft.ifftshift(Zx)
ft = np.fft.fft2(ft)
ft = np.fft.fftshift(ft)
ax.imshow(np.abs(ft))
ax.set_title('No Vortex $|v_x(k)|$')

ax = fig.add_subplot(2,2,2)
ft = np.fft.ifftshift(Zy)
ft = np.fft.fft2(ft)
ft = np.fft.fftshift(ft)
ax.imshow(np.abs(ft))
ax.set_title('No Vortex $|v_y(k)|$')

tree = cKDTree(pos_not_vortex, boxsize=para.L)

ax = fig.add_subplot(2,2,3)
for i in range(100):
    for j in range(100):
        particles_indices = tree.query_ball_point([X[i, j], Y[i, j]], 3)
        Zx[i,j], Zy[i,j] = np.sum(vel_not_vortex[particles_indices],axis=0)[0], np.sum(vel_not_vortex[particles_indices],axis=0)[1]
ft = np.fft.ifftshift(Zx)
ft = np.fft.fft2(ft)
ft = np.fft.fftshift(ft)
ax.imshow(np.abs(ft))
ax.set_title('With Vortex $|v_x(k)|$')

ax = fig.add_subplot(2,2,4)
ft = np.fft.ifftshift(Zy)
ft = np.fft.fft2(ft)
ft = np.fft.fftshift(ft)
ax.imshow(np.abs(ft))
ax.set_title('With Vortex $|v_y(k)|$')
plt.tight_layout()

plt.show()
plt.close()

fig = plt.figure()

firstpoint, lastpoint = 0, -1
xvel_not_vortex, yvel_not_vortex = vel_not_vortex[:,0], vel_not_vortex[:,1]
kxvel_not_vortex, kyvel_not_vortex  = np.fft.fft(xvel_not_vortex), np.fft.fft(yvel_not_vortex)
print(xvel_not_vortex.shape)
freq = np.fft.fftfreq(xvel_not_vortex.shape[-1], 0.005)
print(freq.shape)
print(freq)
print(kxvel_not_vortex.real)
ax = fig.add_subplot(221)
n = 50
limit_value = 100 * 0.01 #5.0e-05
# selectedvalues =np.concatenate((np.arange(n), np.arange(-1,-n,-1)))
selectedvalues =np.arange(freq.shape[-1])
ax.plot(freq[selectedvalues], kxvel_not_vortex.real[selectedvalues], label='No Vortex $Re\{v_x(k)\}$')
is_ = np.where(kxvel_not_vortex.real == np.max(kxvel_not_vortex.real))
print('peaks in not vortex for x', is_[0], freq[is_], kxvel_not_vortex.real[is_])
# ax.plot(freq, kxvel_not_vortex.imag,'--', label='No Vortex $Im\{v_x(k)\}$')
# ax.set_xlim(-limit_value,limit_value)
plt.legend()
ax = fig.add_subplot(222)
ax.plot(freq[selectedvalues], kyvel_not_vortex.real[selectedvalues], label='No Vortex $Re\{v_y(k)\}$')
is_ = np.where(kyvel_not_vortex.real == np.max(kyvel_not_vortex.real))
print('peaks in not vortex for y', is_[0], freq[is_], kyvel_not_vortex.real[is_])
# ax.plot(freq, kyvel_not_vortex.imag,'--', label='No Vortex $Im\{v_y(k)\}$')
# ax.set_xlim(-limit_value,limit_value)
plt.legend()
xvel_with_vortex, yvel_with_vortex = vel_with_vortex[:,0], vel_with_vortex[:,1]
kxvel_with_vortex, kyvel_with_vortex = np.fft.fft(xvel_with_vortex), np.fft.fft(yvel_with_vortex)
ax = fig.add_subplot(223)
ax.plot(freq[selectedvalues], kxvel_with_vortex.real[selectedvalues], label='With Vortex $Re\{v_x(k)\}$')
# ax.plot(freq, kxvel_with_vortex.imag,'--', label='With Vortex $Im\{v_x(k)\}$')
is_ = np.where(kxvel_with_vortex.real == np.max(kxvel_with_vortex.real))
print('peaks in with vortex for x', is_[0], freq[is_], kxvel_with_vortex.real[is_])
# ax.set_xlim(-limit_value,limit_value)
plt.legend()
ax = fig.add_subplot(224)
ax.plot(freq[selectedvalues], kyvel_with_vortex.real[selectedvalues], label='With Vortex $Re\{v_y(k)\}$')
# ax.plot(freq, kyvel_with_vortex.imag,'--', label='With Vortex $Im\{v_y(k)\}$')
is_ = np.where(kyvel_with_vortex.real == np.max(kyvel_with_vortex.real))
print('peaks in with vortex for y', is_[0], freq[is_], kyvel_not_vortex.real[is_])
# ax.set_xlim(-limit_value,limit_value)
plt.legend()
plt.tight_layout()
plt.show()
plt.close()
# pos_not_vortex = np.load('%s_pos.npy' % dirs[0], 'r')
# t = pos_not_vortex.shape[-1]
# fig = plt.figure(figsize=(8,8))
# ax = fig.add_subplot(1, 1, 1)
# ax.scatter(pos_not_vortex[0,:], pos_not_vortex[1,:],s=0.5)
# ax.scatter(pos_not_vortex[0,0], pos_not_vortex[1,0])
# ax.set_xlim(0,para.L)
# ax.set_ylim(0,para.L)
# plt.savefig('pos_not_v.png')
# plt.show()
# plt.close()

# pos_not_vortex = np.load('%s_pos.npy' % dirs[1], 'r')
# t = pos_not_vortex.shape[-1]
# fig = plt.figure(figsize=(8,8))
# ax = fig.add_subplot(1, 1, 1)
# ax.scatter(pos_not_vortex[0,:], pos_not_vortex[1,:],s=0.5)
# ax.scatter(pos_not_vortex[0,0], pos_not_vortex[1,0])
# ax.set_xlim(0,para.L)
# ax.set_ylim(0,para.L)
# plt.savefig('pos_with_v.png')
# plt.show()
# plt.close()

# for part in range(no_of_parts):
#     f = h5py.File('%s/data_particles_1000_ensemble_0_part_%d.hdf5' % (dirs[0],part), 'r')
#     pos = f['particles_positions_with_time']
#     lenofparts = pos.shape[-1]
#     for t in range(0,lenofparts,50):
#         xpos, ypos = pos[:,0,t], pos[:,1,t]
#         plt.plot(xpos,ypos,'k.',lw=0.1)
#         plt.plot(xpos[particle_id], ypos[particle_id], 'r.', lw=0.1)
#         plt.xlim(0,para.L)
#         plt.ylim(0,para.L)
#         plt.savefig('pos%d.png' % ( (part*lenofparts) + t) )
#         plt.close()
        # plt.show()
        # break
    # break


if True:
    raise ValueError('finish')

velnotvortext = np.load('%s.npy' % dirs[1],'r')
number = no_of_parts*100
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(2, 2, 1)

ax.plot(np.arange(number),velnotvortext[0,:number],label='not vortex xvel')
ax.plot(np.arange(number),velnotvortext[1,:number],label='not vortex yvel')
plt.legend()
ax = fig.add_subplot(1,2,2)

fourier = np.fft.fft(velnotvortext[0,:])# .rfft2(velnotvortext)
# arr = r_
# fourier = np.fft.rfft2(np.sin(r_))
# print(r_)
print('from here')
print('shape',fourier.shape)
print(fourier)
freq = np.fft.fftfreq(fourier.shape[-1])
plt.plot(freq, fourier)
print(freq.shape)
print(freq)
# print(np.fft.rfft2(velnotvortext))

# ax.plot(freq,fourier[0,:],label='with vortex xvel')
# ax.plot(freq,fourier[1,:],label='with vortex yvel')
# ax.quiver(freq, freq, fourier[0,:], fourier[1,:])# plot(freq,fourier[1,:],label='with vortex yvel')
# ax.plot(np.arange(number),velwithvortext[1,:number],label='with vortex xvel')
plt.legend()
# # plt.plot(np.arange(69900),velwithvortext[0,:])
plt.show()
# f1 = h5py.File('data_particles_1000_ensemble_0_part_697.hdf5','r')
# f2 = h5py.File('data_particles_1000_ensemble_0_part_699.hdf5','r')
# # print(f1.keys())
# # pos = f['particles_positions_with_time']

# particle_id = 400
# vel1 = f1['particles_velocities_with_time'][particle_id,:,:]
# vel2 = f2['particles_velocities_with_time'][particle_id,:,:]
# print(vel1.shape)
# arr = np.concatenate((vel1,vel2), axis=1)
# print(arr.shape)
# print(arr)
# t = -1
# tree = cKDTree(pos[:,:,t], boxsize=para.L)
# v_ = gridPointsVelocities(r_, tree, vel[:,:,t])
# displayVelField(r_,v_,show=True)