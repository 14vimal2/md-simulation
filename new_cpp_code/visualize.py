import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from progress.bar import Bar 
# import regex

# make backgroud dark
plt.style.use('dark_background')

para_file_name = 'para.h'

f = open(para_file_name, 'r')

def present(txt : str):
    txt = txt.strip()
    if '#' in txt and txt[0] == '#':
        txt = txt[1:].strip()
        if txt[:6] == 'define':
            return 'NUMBER_OF_PARTICLES' in txt or \
                'DIMENSION' in txt or \
                'TIMESTEP' in txt or \
                'FINAL_TIME' in txt or \
                'DENSITY' in txt
    return False



def extractVar(txt : str):
    if '//' in txt:
        txt = txt[:txt.index('//')].strip()
    txt = txt[1:].strip()[6:]

    length = len(txt)
    ans = ''
    for i in range(1, length):
        if txt[i] == ' ' and txt[i-1] == ' ':
            continue
        else:
            ans += txt[i]
    return ans

    


lines = f.readlines()
lines = [ extractVar(txt) for txt in lines  if present(txt)]

variables = {}
for w in lines:
    ws = w.split(' ')
    variables[ws[0]] = ws[1]

N = int(variables['NUMBER_OF_PARTICLES'])
density = float(variables['DENSITY'])
dimension = int(variables['DIMENSION'])

L = (N/ density)**(1/dimension)
file_name = 'data.out'
data = np.loadtxt(file_name, delimiter=' ', skiprows=0)

# print(data.shape)

timesteps = data.shape[0]//N


jump = 1

# create animation

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot()


ax.plot([0, L], [0, 0], 'w', linewidth=2)
ax.plot([0, 0], [0, L], 'w', linewidth=2)
def init_func():
    pass


bar = Bar('Processing...', max=timesteps//jump)
def update(i):
    ax.clear()
    i = i*jump
    x = data[i*N: (i+1)*N, 0]
    y = data[i*N: (i+1)*N, 1]
    ax.set_axis_off()
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)


    bar.next()
    return ax.plot(x, y,'o')

ani = animation.FuncAnimation(fig, init_func=init_func, func=update, frames=timesteps // jump, interval=5)

# save animation as mp4
ani.save('movie.mp4', writer='ffmpeg')

bar.finish()
