import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# defining relationship y = f(x)
def f(x, a, b, c):
    return a*x**2 + b * x + c

def f1(x, a, b):
    return a * x * np.log2(x) + b

xrealinputs = np.linspace(0,10,10) 
noise = 5*np.random.rand(10)
yreal_data = f(xrealinputs, 1, 2, 1.2) + noise


xrealinputs = [100, 500, 1000, 2500, 5000  ]
yreal_data1 = [
    (0.42107295989990234+0.3773813247680664+0.37511444091796875)/3,
    (0.8477673530578613+0.846865177154541+0.8530044555664062)/3,
    (1.419421672821045+1.4457647800445557+1.4076883792877197)/3,
    (2.9106132984161377+2.9690301418304443+2.9845783710479736)/3,
    (5.545851945877075+5.375452280044556+5.683040618896484)/3,

]
yreal_data = [
    (0.7769174575805664+0.7086718082427979+0.7318031787872314)/3,
    (44.59126806259155),
    (234.98265266418457),
    (1408.1233713626862),
    (5006.451021671295),

]




xfitinput = np.linspace(np.min(xrealinputs), np.max(xrealinputs), 100)
params, _ = curve_fit(f, xrealinputs, yreal_data)
print(params)
yfit_data = f(xfitinput, *params)
params, _ = curve_fit(f1, xrealinputs, yreal_data1)
print(params)
yfit_data1 = f1(xfitinput, *params)



plt.title('square matrix vs tree', fontsize=20)
plt.plot(xrealinputs,yreal_data,'ro',label='realdata $O(N^2)$')
plt.plot(xfitinput, yfit_data, 'k--',label='fitdata')
plt.plot(xrealinputs, yreal_data1, 'bs', label='realdata $O(NlogN)$')
plt.plot(xfitinput, yfit_data1, 'k--',label='fitdata')
plt.xlabel('x',fontsize=15)
plt.ylabel('$f(x)$',fontsize=15)
plt.legend()
plt.show()
plt.close()