from math import ceil
import numpy as np
import para
from scipy.stats import maxwell

if para._device_ == "gpu":
    import cupy as xp
    para.gpudevice.use()
else:
    xp = np


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
        self.r = xp.zeros((para.N, para.dim), dtype=xp.float32) 
        self.v = xp.zeros((para.N, para.dim), dtype= xp.float32)
        self.rp = xp.zeros((para.N, para.dim), dtype = xp.float32)
        self.vcm = xp.zeros((1, para.dim), dtype = xp.float32)
        self.F = xp.zeros((para.N, para.dim), dtype = xp.float32)
        self.msv = 0
        self.PE = 0
        self.KE = 0
        self.T = 0
        if para._device_ == "gpu":
            self.d2matrix = xp.zeros((para.N, para.N), dtype = xp.float32)


    def init(self):
        if True:
        # set particle positions
            # set positions at lattice points
            sn = ceil(para.N**(1/ para.dim))
            # print('sn',sn)
            index_ = xp.arange(para.N)
            for j in range(para.dim):
                self.r[:,j] = (index_ // sn**j) % sn
            self.r = ((self.r + 1/2)/ sn )* para.L
        else:
            self.r = pos
        if True:
            # set particle velocities
            vels = maxwell.rvs(loc=0,scale=np.sqrt(para.T), size=para.N)
            if para._device_ == "gpu":
                vels = xp.array(vels, dtype=xp.float32)
            velocities = vels - xp.average(vels)
            if (para.dim == 2):
                thetha = xp.random.rand(para.N) * 2 * np.pi
                self.v = xp.zeros((para.N, para.dim))
                xvel = (xp.sin(thetha) * velocities).reshape(para.N,1)
                yvel = (xp.cos(thetha) * velocities).reshape(para.N,1)
                self.v = xp.hstack((xvel,yvel))
            if (para.dim == 3):
                thetha = xp.random.rand(para.N) * np.pi
                phi = xp.random.rand(para.N) * 2 * np.pi
                xvel = (velocities * xp.sin(thetha) * xp.cos(phi)).reshape(para.N,1)
                yvel = (velocities * xp.sin(thetha) * xp.sin(phi)).reshape(para.N,1)
                zvel = (velocities * xp.cos(thetha)).reshape(para.N).reshape(para.N,1)
                self.v = xp.hstack((xvel, yvel, zvel))

            self.vcm = xp.sum(self.v, axis=0, dtype=xp.float32)/ para.N
            self.KE = xp.sum(xp.sum(self.v**2,axis=0, dtype=xp.float32))
            scalefactor = (para.dim * para.T * para.N / self.KE) ** (1/2)
            self.v = (self.v - self.vcm) * scalefactor
        else:
            self.v = xp.array(vel)
        print('max velocity', xp.max(xp.sqrt(xp.sum(self.v**2,axis=1))))

        # set particle previous position
        self.rp = self.r - self.v * para.dt
    
    
