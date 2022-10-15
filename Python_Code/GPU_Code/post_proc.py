import numpy as np
import para
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def pdf3D(x,a):
    return np.sqrt(2/np.pi) * (x**2) * np.exp(- x**2 / (2 * a**2)) / a**3

def pdf2D(x,a):
    return x * np.exp(- x**2 / (2 * a **2)) / a**2 * (x[1] - x[0])

def pdftheoretical3D(v,T):
    return (1/(2 * np.pi * T))**1.5 * 4 * np.pi * v**2 * np.exp(-v**2/ (2* T))

def pdftheoretical2D(v,T):
    return v/T * np.exp(-v**2 / (2 * T))*(v[1]-v[0])


def showVelDist(vel, color='b' ,bins=100, compare = False):
    velocities = np.sum(vel**2,axis=1)**(1/2)
    A, B = np.histogram(velocities,bins=bins,density=True)
    A /= np.sum(A)
    B = (B[0:-1] + B[1:])/2

    params, _ = curve_fit(pdf2D, B, A)
    print(params)
    x = np.linspace(0, np.max(B), 100)
    y = pdf2D(x,*params)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    if not compare:
        ax.plot(B,A,'r.',lw=0.5,label='simulation data')
        ax.plot(x,y,'k--',label='maxwell fit')
        Tprime = np.sum(velocities**2)/2/para.N
        print(para.T,Tprime)
        ax.plot(x,pdftheoretical2D(x,Tprime),'b',label='Theoretical prediction')
    else:
        ax.plot(B,A,color + '.',lw=0.5)
        ax.plot(x,y,color + '--')
        Tprime = np.sum(velocities**2)/2/para.N
        print(para.T,Tprime)
        ax.plot(x,pdftheoretical2D(x,Tprime),color,label='T = {}'.format(para.T))
    ax.set_xlabel('$v$',fontdict={"size":15})
    ax.set_ylabel('$P(v)$', fontdict={'size':15})
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0)
    ax.legend()
    plt.show()
    plt.close()


def getVel(prevPos,nextPos):
    rdiff = nextPos - prevPos
    rdiff = rdiff + para.L * (rdiff < -para.L/2) - para.L * (rdiff > para.L/2)
    return rdiff/para.dt/2.0

def displayVelField(r,v,figname=None,show=False, dirname=None):
    fig = plt.figure()
    ax = fig.add_subplot(111,aspect='equal')
    ax.quiver(r[:,0], r[:,1], v[:,0], v[:,1])
    veltickno = 6
    ax.set_xticks(np.linspace(0,para.L,veltickno))
    ax.set_xticklabels(np.round(np.linspace(0,para.L,veltickno,dtype=int)),fontsize=15)
    ax.set_yticks(np.linspace(0,para.L,veltickno))
    ax.set_yticklabels(np.round(np.linspace(0,para.L,veltickno,dtype=int)),fontsize=15)
    plt.xlim(0,para.L)
    plt.ylim(0,para.L)
    if figname !=None:
        plt.title(figname[:-4],fontsize=15)
        if dirname != None:
            plt.savefig(dirname + '/' + figname)
        else:
            plt.savefig(figname)

    if show:
        plt.show()
    plt.close()
    pass

def showEnergyWithTime(xarr, KE_arr,  PE_arr):
    plt.title('Energy vs Time of Particles')
    plt.plot(xarr, PE_arr, label=r'$PE$')
    plt.plot(xarr, KE_arr, label=r'$KE$')
    plt.xlabel(r'time ($t^*=t\sqrt{ \frac{\epsilon}{m \sigma^2}  }$)', fontsize=10)
    plt.ylabel(r'energy ($E^*= \frac{E}{\epsilon}$)', fontsize=10)
    plt.xlim(0,xarr[-1])
    plt.ylim(np.min(PE_arr)-30,np.max(KE_arr)+ 30)
    plt.margins(20,40)
    E_mean = np.average(KE_arr + PE_arr)
    print(E_mean)
    std_dev_E = np.sqrt(np.average((KE_arr+PE_arr-E_mean)**2))
    print(std_dev_E)
    plt.plot(xarr, PE_arr + KE_arr,'k--', label=r'$KE + PE$')
    plt.legend(loc=7)
    plt.savefig('energyvstime.png')

def coarseGraining(r, v, n: int,averaged: bool=True,filename=None):
    if n*n >= para.N:
        raise ValueError('n should be less than sqrt(N)')
    ln = para.L/n
    rprime = np.zeros((n**2,para.dim))
    index_ = np.arange(n**2)
    for j in range(para.dim):
        rprime[:,j] = (index_//n**j) % n
    rprime = ((rprime+0.5)/n) * para.L
    vprime = np.zeros((n**2,para.dim))
    counts = np.zeros(n**2)
    # r, v = cp.asnumpy(P.r), cp.asnumpy(P.v)
    for i in range(para.N):
        xind, yind = int(r[i,0]/ ln), int(r[i,1]/ln)
        vprime[yind*n + xind] += v[i]
        counts[yind*n + xind] += 1
    if averaged:
        vprime[:,0], vprime[:,1] = vprime[:,0]/counts, vprime[:,1]/counts
    displayVelField(rprime,vprime,filename)
    pass

