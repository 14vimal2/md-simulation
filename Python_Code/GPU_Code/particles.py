from math import ceil
import numpy as np
import para
import cupy as cp

para.device.use()

class Particles:
    def __init__(self):
        self.r = []
        self.v = []
        self.rp = []
        self.vcm = []
        self.msv = None
        self.PE = None
        self.KE = None
        self.F = []
        self.T = None
        self.d2matrix = []
        self.iarr = []
        self.jarr = []
        self.r2arr = []
        self.collisionnumber = 0
    
    def set_arrays(self):
        self.r = cp.zeros((para.N, para.dim), dtype=cp.float32) 
        self.v = cp.zeros((para.N, para.dim), dtype= cp.float32)
        self.rp = cp.zeros((para.N, para.dim), dtype = cp.float32)
        self.vcm = cp.zeros((1, para.dim), dtype = cp.float32)
        self.F = cp.zeros((para.N, para.dim), dtype = np.float32)
        self.msv = 0
        self.PE = 0
        self.KE = 0
        self.T = 0
        self.d2matrix = cp.zeros((para.N, para.N), dtype = np.float32)


    def init(self, pos, vel):
        # if posFileName==None:
        # # set particle positions
        #     # set positions at lattice points
        #     sn = ceil(para.N**(1/ para.dim))
        #     # print('sn',sn)
        #     index_ = cp.arange(para.N)
        #     for j in range(para.dim):
        #         self.r[:,j] = (index_ // sn**j) % sn
        #     self.r = ((self.r + 1/2)/ sn )* para.L
        # else:
        #     self.r = cp.array(posFileName)
        # if velFileName==None:
        #     # set particle velocities
        #     velocities = cp.random.rand(para.N, dtype=cp.float32) - 1/2
        #     if (para.dim == 2):
        #         thetha = cp.random.rand(para.N) * 2 * np.pi
        #         self.v = cp.zeros((para.N, para.dim))
        #         xvel = (cp.sin(thetha) * velocities).reshape(para.N,1)
        #         yvel = (cp.cos(thetha) * velocities).reshape(para.N,1)
        #         self.v = cp.hstack((xvel,yvel))
        #     if (para.dim == 3):
        #         thetha = cp.random.rand(para.N) * np.pi
        #         phi = cp.random.rand(para.N) * 2 * np.pi
        #         xvel = (velocities * cp.sin(thetha) * cp.cos(phi)).reshape(para.N,1)
        #         yvel = (velocities * cp.sin(thetha) * cp.sin(phi)).reshape(para.N,1)
        #         zvel = (velocities * cp.cos(thetha)).reshape(para.N).reshape(para.N,1)
        #         self.v = cp.hstack((xvel, yvel, zvel))

        #     self.vcm = cp.sum(self.v, axis=0, dtype=cp.float32)/ para.N
        #     self.KE = np.sum(cp.asnumpy(cp.sum(self.v**2,axis=0, dtype=cp.float32)))
        #     scalefactor = (para.dim * para.T * para.N / self.KE) ** (1/2)
        #     self.v = (self.v - self.vcm) * scalefactor
        # else:
        #     self.v = cp.array(velFileName)
        self.r = cp.array(pos)
        self.v = cp.array(vel)
        print('max velocity', cp.max(cp.sqrt(cp.sum(self.v**2,axis=1))))

        # set particle previous position
        self.rp = self.r - self.v * para.dt
    
    
