
import os

# all physical quantities in LJ units

N = 1000                # No. of Particle

T =  0.306            # temperature
#vrms = 0
dt = 0.005                # timestep
tf = 20000                # final time
dim =  2                # dimension
rc  = 2.5                  # cutoff
rho = 3.04e-3              # density of particle


# for Vortices entry
add_vortices = True
eta = 1
nvorticesrt = 2


# no of ensemble
no_of_ensemble = 1

# no of snapshots 
# use int(tf/dt) for taking snapshot at each interval
no_of_snapshots = int(tf/dt)//200

# directory for data saving
dir_name = 'simulation_data_not_v_new'
if not os.path.exists(dir_name):
   os.makedirs(dir_name)

# cuda device
try:
    from cupy.cuda import Device
    no_of_devices = 7
    curr_device_index = 3 
    device = Device(curr_device_index)
except Exception as e:
    print(e)


#-----------Not to edit------------

L = (N/rho)**(1/dim)   # size of box
PE_cut = 4 * (1/rc**12 - 1/rc**6) # lj potential at rc

