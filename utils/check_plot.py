import enum
import numpy as np
import matplotlib.pyplot as plt
import sys as sys


fname = sys.argv[1]

idx, x, y = np.loadtxt(fname, unpack=True, usecols=(0,1,2))
idx = (np.asarray(idx, dtype=np.int32))
# print(idx)

def group_data(idx, x, y, grid):
    xg, yg = [], []
    cnt = 0;
    for i, id in enumerate(idx):
        if(id == grid):
            xg.append(x[i])
            yg.append(y[i])
            cnt = cnt + 1;
    return np.asarray(xg), np.asarray(yg), cnt

def plot_grouped(xg, yg, cnt, fig, 
                 ax, clr, smb, ms=10):
    for i in range(0, cnt):
        ii = np.asarray([0, (i+1)%cnt, (i+2)%cnt])
        xn_, yn_  = xg[ii], yg[ii] ;
        xn__, yn__  = xg[ii], yg[ii] ;
        for i in [1, 2]:
            dx = xn_[i] - xn_[0]
            if(dx >= 0.5): xn_[i] =  xn_[i] - 1  
            if(dx < -0.5): xn_[i] =  xn_[i] + 1  
            dy = yn_[i] - yn_[0]
            if(dy >= 0.5): yn_[i] =  yn_[i] - 1  
            if(dy < -0.5): yn_[i] =  yn_[i] + 1  
        ax.plot(xn_, yn_, '-', marker=smb, color=clr, mfc='none', ms=ms)

    ii = np.asarray([0, cnt-1, 1 ])
    xn_, yn_  = xg[ii], yg[ii] ;
    for i in [1, 2]:
        dx = xn_[i] - xn_[0]
        if(dx >= 0.5): xn_[i] =  xn_[i] - 1  
        if(dx < -0.5): xn_[i] =  xn_[i] + 1  
        dy = yn_[i] - yn_[0]
        if(dy >= 0.5): yn_[i] =  yn_[i] - 1  
        if(dy < -0.5): yn_[i] =  yn_[i] + 1  

    ax.plot(xn_, yn_, '-', marker=smb, color=clr, mfc='none', ms=ms)


fig, ax = plt.subplots()
x1, y1, cnt = group_data(idx, x, y, 1)
plot_grouped(x1, y1, cnt, fig, ax, 'tab:blue', 'o', ms = 20)
x1, y1, cnt = group_data(idx, x, y, 2)
plot_grouped(x1, y1, cnt, fig, ax, 'tab:orange', 's', ms = 15)
x1, y1, cnt = group_data(idx, x, y, 3)
plot_grouped(x1, y1, cnt, fig, ax, 'tab:green', 'P', ms = 10)
x1, y1, cnt = group_data(idx, x, y, 4)
plot_grouped(x1, y1, cnt, fig, ax, 'tab:red', 'd', ms = 5)

# x1, y1 = np.loadtxt('../Examples/cart_16.dat', unpack=True, usecols=(0,1))
# ax.plot(x1, y1, 'o', color = 'tab:blue', ms = 1)
# ax.plot(x1-1, y1, 'o', color = 'tab:blue',ms = 1)
# ax.plot(x1, y1-1, 'o', color = 'tab:blue',ms = 1)
# ax.plot(x1-1, y1-1, 'o', color = 'tab:blue',ms = 1)
# ax.plot(x1+1, y1-1, 'o', color = 'tab:blue',ms = 1)
# ax.plot(x1+1, y1, 'o', color = 'tab:blue',ms = 1)
# ax.plot(x1, y1+1, 'o', color = 'tab:blue',ms = 1)
# ax.plot(x1+1, y1+1, 'o', color = 'tab:blue',ms = 1)
# ax.plot(x1-1, y1+1, 'o', color = 'tab:blue',ms = 1)
# plt.grid(True)
# ofile = fname.replace('dat', 'png')
# fig.savefig(ofile)
plt.show()
