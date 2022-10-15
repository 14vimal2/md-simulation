from math import ceil
import numpy as np
import para


class Particles:
    def __init__(self) -> None:
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
        self.r = np.zeros((para.N, para.dim))
        self.v = np.zeros((para.N, para.dim))
        self.rp = np.zeros((para.N, para.dim))
        self.vcm = np.zeros((1,para.dim))
        self.F = np.zeros((para.N, para.dim))
        self.msv = 0
        self.PE = 0
        self.KE = 0
        self.T = 0
        self.d2matrix = np.zeros((para.N, para.N))


    def init(self):
        # set positions at lattice points
        sn = ceil(para.N**(1/ para.dim))
        # print('sn',sn)
        index_ = np.arange(para.N)
        for j in range(para.dim):
            self.r[:,j] = (index_ // sn**j) % sn
        self.r = ((self.r + 1/2)/ sn )* para.L

        # set particle velocities
        velocities = np.random.rand(para.N) - 1/2

        if (para.dim == 2):
            thetha = np.random.rand(para.N) * 2 * np.pi
            self.v = np.zeros((para.N, para.dim))
            xvel = (np.sin(thetha) * velocities).reshape(para.N,1)
            yvel = (np.cos(thetha) * velocities).reshape(para.N,1)
            self.v = np.hstack((xvel,yvel))
        if (para.dim == 3):
            thetha = np.random.rand(para.N) * np.pi
            phi = np.random.rand(para.N) * 2 * np.pi
            xvel = (velocities * np.sin(thetha) * np.cos(phi)).reshape(para.N,1)
            yvel = (velocities * np.sin(thetha) * np.sin(phi)).reshape(para.N,1)
            zvel = (velocities * np.cos(thetha)).reshape(para.N).reshape(para.N,1)
            self.v = np.hstack((xvel, yvel, zvel))


        self.vcm = np.sum(self.v, axis=0)/ para.N
        self.KE = np.sum(self.v**2)
        scalefactor = (para.dim * para.T * para.N / self.KE) ** (1/2)
        self.v = (self.v - self.vcm) * scalefactor
        print('Vrms',(np.sum(self.v**2)/para.N)**(1/2))

        # set particle previous position
        self.rp = self.r - self.v * para.dt
    
    