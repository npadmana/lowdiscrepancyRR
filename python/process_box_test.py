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
weightN = {}

#basedir = '../data/counts/'
basedir = '../data/deep2/counts/'

print 'bins vs. normalized counts'
for i in 1000*(2**np.arange(10)):
    files = glob.glob(basedir+str(i)+'_*.dat')
    if files:
        data[i] = []
        for f in files:
            if 'deep2' in basedir:
                # need to get the sum of weights
                ff = open(f,'r')
                weightN = float(ff.readline())
            else:
                ff = open(f,'r')
                weightN = 1
            temp = np.loadtxt(ff,skiprows=1,dtype=np.float64)
            data[i].append(temp[:,2]/weightN**2)
            bins = temp[:,0]
        data[i] = np.array(data[i])

xx = np.array(sorted(data.keys()))

sigma = {}
sigmabin = {}
mean = {}
sigma0 = []
sigma1 = []
sigma2 = []
sigma3 = []
sigma4 = []
sigmax = []
for x in xx:
    meanx = data[x].mean(axis=0)
    stdx = data[x].std(axis=0)
    sigma[x] = stdx/meanx
    mean[x] = meanx
    sigma0.append(sigma[x][0])
    sigma1.append(sigma[x][1])
    sigma2.append(sigma[x][2])
    sigma3.append(sigma[x][3])
    sigma4.append(sigma[x][4])

for i in xx:
    print i
    if 'deep2' in basedir:
        for b,x in zip(bins,mean[i]):
            print b,x
    else:
        for b,x in zip(bins,mean[i]):
            print b,x/i**2
    print '-------------'

    
plt.ylabel(r'$\sigma_{RR}/mean(RR)$')
plt.xlabel(r'$N_R$')
# some horizontal lines
plt.plot([1e3,1e6],[1e-2,1e-2],lw=0.5,color='grey')
plt.plot([1e3,1e6],[1e-3,1e-3],lw=0.5,color='grey')
plt.plot([1e3,1e6],[1e-4,1e-4],lw=0.5,color='grey')

# comparison lines
plt.loglog(xx,1/np.sqrt(xx/1e3)*0.025,'--',lw=3,label=r'$1/\sqrt{N_R}$')
plt.loglog(xx,1/(xx/1e3)*0.05,'--',lw=3,label=r'$1/N_R$')

plt.loglog(xx,sigma0,'-^',lw=1.5,label='1st bin')
plt.loglog(xx,sigma1,'-v',lw=1.5,label='2nd bin')
plt.loglog(xx,sigma2,'-*',lw=1.5,label='3rd bin')
plt.loglog(xx,sigma3,'-s',lw=1.5,label='4th bin')
plt.loglog(xx,sigma4,'-o',lw=1.5,label='5th bin')
plt.legend()
plt.xlim(9e2,6e5)
plt.ylim(1e-4,2e-1)
plt.subplots_adjust(bottom=0.14)
if 'deep2' in basedir:
    plt.savefig('../plots/sigma_RR-mr_pairs-deep2.pdf')
else:
    plt.savefig('../plots/sigma_RR-mr_pairs-box.pdf')

import pdb
pdb.set_trace()
