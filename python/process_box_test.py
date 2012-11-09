#!/usr/bin/env python
"""
Process the results of box_test.py and make some plots.
"""
import glob
import numpy as np
import matplotlib
from matplotlib import rc
import matplotlib.pyplot as plt
plt.rcParams.update({'backend':'pdf',
                     'text.fontsize':18,
                     'axes.labelsize': 18,
                     'legend.fontsize': 14,
                     'xtick.labelsize': 18,
                     'ytick.labelsize': 18,
                     'axes.linewidth':2})
plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

plt.ion()

data = {}

for i in 1000*(2**np.arange(10)):
    files = glob.glob('counts/'+str(i)+'_*.dat')
    if files:
        data[i] = []
        for f in files:
            temp = np.loadtxt(f,skiprows=1,dtype=np.float128)
            data[i].append(temp[:,2])
        data[i] = np.array(data[i])

xx = np.array(sorted(data.keys()))
sigma = {}
sigmabin = {}
sigma0 = []
sigma1 = []
sigma2 = []
sigma3 = []
for x in xx:
    sigma[x] = data[x].std(axis=0)/data[x].mean(axis=0)
    sigma0.append(sigma[x][0])
    sigma1.append(sigma[x][1])
    sigma2.append(sigma[x][2])
    sigma3.append(sigma[x][3])

plt.ylabel(r'$\sigma_{RR}/mean(RR)$')
plt.xlabel(r'$N_R$')
plt.loglog(xx,sigma0,'-^',lw=1.5,label='1st bin')
plt.loglog(xx,sigma1,'-v',lw=1.5,label='2nd bin')
plt.loglog(xx,sigma2,'-*',lw=1.5,label='3rd bin')
plt.loglog(xx,sigma3,'-s',lw=1.5,label='4th bin')
plt.loglog(xx,1/np.sqrt(xx/1e3)*0.025,lw=4)
plt.legend()
plt.xlim(9e2,3e5)
plt.ylim(1e-3,2e-1)

import pdb
pdb.set_trace()
