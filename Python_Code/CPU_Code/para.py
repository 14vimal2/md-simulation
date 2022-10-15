

# all physical quantities in LJ units

N = 500                # No. of Particle

T =  0.306            # temperature
dt = 0.005                # timestep
tf = 10                # final time
dim =  2                # dimension
rc  = 2.5                  # cutoff
rho = 3.04e-3              # density of particle


#-----------Not to edit------------

L = (N/rho)**(1/dim)   # size of box
PE_cut = 4 * (1/rc**12 - 1/rc**6) # lj potential at rc

# print(L)

# print(1/rc**dim)
