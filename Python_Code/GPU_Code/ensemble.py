import para
import time
from fns import *
from particles import Particles
import h5py

para.device.use()

import logging
logging.basicConfig(filename='{}/ensemble.log'.format(para.dir_name), filemode='w', level=logging.DEBUG,format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')
try:
    t1 = time.time()
    ensemble = []
    files = []
    data_part_length = 100
    data_part = 0
    data_part_index = 0
    logging.info('Ensemble creation has started')
    # inputfiles = [] 
    for i in range(para.no_of_ensemble):
        # inputfiles.append(
        #     h5py.File("%s/data_particles_%d_ensemble_%d_part_%d.hdf5"%(para.dir_name, para.N,i,data_part-1), "r")
        # )
        files.append(
            h5py.File("%s/data_particles_%d_ensemble_%d_part_%d.hdf5"%(para.dir_name, para.N,i,data_part), "w")
            )
        ensemble.append(Particles())
        ensemble[i].set_arrays()
        ensemble[i].init(
            pos=np.load('{}/mb_positions_{}_{}_{}_{}.npy'.format('mb_dist_data',para.N,para.dim,para.T,i)),
            vel=np.load('{}/mb_velocities_{}_{}_{}_{}.npy'.format('mb_dist_data',para.N,para.dim,para.T,i))
        )
        # ensemble[i].init(
        #     pos=inputfiles[i]['particles_positions_with_time'][:,:,-1],
        #     vel=inputfiles[i]['particles_velocities_with_time'][:,:,-1]
        # )
        # inputfiles[i].close()
        if para.add_vortices:
            addVortices(ensemble[i],para.eta,True,para.nvorticesrt)
    logging.info('{} Ensembles are created'.format(para.no_of_ensemble))
    inputdata = open('para.py','r').read()
    logging.info(inputdata)
    simulation_program = open('ensemble.py', 'r').read()
    logging.info(simulation_program)
    ti = 0
    n = int(para.tf / para.dt)
    times = np.arange(0,n+n//para.no_of_snapshots,n//para.no_of_snapshots)
    logging.info('snapshot timings are {}'.format(times))
    positions_with_time = cp.empty((para.N,para.dim,data_part_length,para.no_of_ensemble),dtype=cp.float32)
    velocities_with_time = cp.empty((para.N,para.dim,data_part_length,para.no_of_ensemble),dtype=cp.float32)
    while (ti < n):
        for i in range(para.no_of_ensemble):            
            compute_interacting_pairs(ensemble[i])
            compute_force(ensemble[i])
            time_advance_single_step(ensemble[i], para.dt)
        if ti in times:
            for i in range(para.no_of_ensemble):            
                positions_with_time[:,:,data_part_index,i] = ensemble[i].r
                velocities_with_time[:,:,data_part_index,i] = ensemble[i].v
            data_part_index += 1
            if data_part_index == data_part_length:
                data_part_index = 0
                data_part += 1
                for i in range(para.no_of_ensemble):
                    files[i].create_dataset("particles_positions_with_time",data=cp.asnumpy(positions_with_time[:,:,:,i]))
                    files[i].create_dataset("particles_velocities_with_time",data=cp.asnumpy(velocities_with_time[:,:,:,i]))
                    files[i].close()
                files.clear()
                positions_with_time *= 0
                for i in range(para.no_of_ensemble):
                    files.append(h5py.File("%s/data_particles_%d_ensemble_%d_part_%d.hdf5"%(para.dir_name, para.N,i,data_part), "w"))
        ti += 1


    average_collision_no = 0
    for i in range(para.no_of_ensemble):
        average_collision_no += ensemble[i].collisionnumber
        files[i].close()

    logging.info('average_collision_number = {}'.format(average_collision_no/para.no_of_ensemble))
    t2 = time.time()
    logging.info('time take in simulation = {}'.format(t2 - t1))
except Exception as e:
    logging.error('Exception occured', exc_info=True)