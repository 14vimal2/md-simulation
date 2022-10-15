import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import matplotlib.path as path
import matplotlib.patches as patches

stepsize = 375
binsno = 40
vmax = 15


for step in range(3375,3380,stepsize):
    data = np.load('vel_900_2.84_370.npy')
    fig, ax = plt.subplots()
    velocities2 = data[np.where(data < 8)]
    # vmax = np.max(velocities2)
    # bins = np.linspace(0,vmax,binsno)
    a, b, *c = ax.hist(velocities2,bins=binsno,density=True, alpha=0.2,color='y')
    b = (b[0:-1] + b[1:])/2
    maxwell = stats.maxwell
    params = maxwell.fit(velocities2,floc=0)
    print(params)
    vmax = np.max(velocities2)
    x = np.linspace(0,vmax, 100)
    ax.plot(b,a, label= 'speed')
    ax.set_xlim(np.min(b),np.max(b))
    y = maxwell.pdf(x,*params)
    ratio = 1 # np.max(a)/np.max(y)
    ax.plot(x,y * ratio,'r--',lw=1,label='maxwell fit')
    plt.title('Speed distribution of particles')
    ax.set_xlabel('v')
    ax.set_ylabel('P(v)')
    plt.legend()
    plt.show()
    plt.close()
    # plt.savefig('veldist.png'.format(step))



# fig = plt.figure()
# speed_ax = fig.add_subplot()
# class Histogram:
#     """A class to draw a Matplotlib histogram as a collection of Patches."""

#     def __init__(self, data, xmax, nbars, density=False):
#         """Initialize the histogram from the data and requested bins."""
#         self.nbars = nbars
#         self.density = density
#         # self.bins = np.linspace(0, xmax, nbars)
#         self.hist, bins = np.histogram(data, binsno, density=density)
#         print(self.hist, bins)

#         # Drawing the histogram with Matplotlib patches owes a lot to
#         # https://matplotlib.org/3.1.1/gallery/animation/animated_histogram.html
#         # Get the corners of the rectangles for the histogram.
#         self.left = np.array(bins[:-1])
#         self.right = np.array(bins[1:])
#         self.bottom = np.zeros(len(self.left))
#         self.top = self.bottom + self.hist
#         nrects = len(self.left)
#         self.nverts = nrects * 5
#         self.verts = np.zeros((self.nverts, 2))
#         self.verts[0::5, 0] = self.left
#         self.verts[0::5, 1] = self.bottom
#         self.verts[1::5, 0] = self.left
#         self.verts[1::5, 1] = self.top
#         self.verts[2::5, 0] = self.right
#         self.verts[2::5, 1] = self.top
#         self.verts[3::5, 0] = self.right
#         self.verts[3::5, 1] = self.bottom
#         # self.patch = None
        

#     def draw(self, ax):
#         """Draw the histogram by adding appropriate patches to Axes ax."""
#         codes = np.ones(self.nverts, int) * path.Path.LINETO
#         codes[0::5] = path.Path.MOVETO
#         codes[4::5] = path.Path.CLOSEPOLY
#         barpath = path.Path(self.verts, codes)
#         self.patch = patches.PathPatch(barpath, fc='tab:green', ec='k', lw=0.5, alpha=0.5)
#         ax.add_patch(self.patch)

#     def update(self, data):
#         """Update the rectangle vertices using a new histogram from data."""
#         self.hist, bins = np.histogram(data, self.bins, density=self.density)
#         self.top = self.bottom + self.hist
#         self.verts[1::5, 1] = self.top
#         self.verts[2::5, 1] = self.top
# speeds = np.load('vel_20000_2.84_3500.npy')
# speed_hist = Histogram(speeds, np.max(speeds), 50, density=True)
# speed_hist.draw(speed_ax)
# plt.show()