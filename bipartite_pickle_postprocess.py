# -*- coding: utf-8 -*-
"""
Created on Fri May 19 23:58:22 2017

@author: Alex
"""
from __future__ import division
import cPickle as pickle
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import gzip
from itertools import product

def draww(X,Y,Z,xlabel,ylabel,zlabel,title):
    fig = plt.figure()
    ax = fig.gca(projection='3d')

#    X=hit_range
#    Y = health_threshold
    X,Y = np.meshgrid(X,Y)
    surf = ax.plot_surface(X,Y,Z.as_matrix().T,cmap=cm.coolwar,rstride=1,cstride=1,linewidth=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.set_title(title)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()


def sum_up_old(lst):
    health=0
    iters=0
    avdeg = 0
    max_clust = 0
    for i,itm in enumerate(lst):
        health+=itm[itm.keys()[0]][0]
        iters+=itm[itm.keys()[0]][1]
        avdeg+=itm[itm.keys()[0]][2]
        max_clust+=itm[itm.keys()[0]][3]
    return health/(i+1.0), iters/(i+1.0), avdeg/(i+1.0), max_clust/(i+1.0)


def sum_up_new(lst):
    health = 0
    iters = 0
    avdeg = 0
    max_clust = 0
    crash_size = 0
    comp_size = 0
    for itm in lst:
        health += itm[0]
        iters += itm[1]
        avdeg += itm[2]
        max_clust += itm[3]
        crash_size += itm[4].fillna(0)
        comp_size += itm[5].fillna(method='ffill')
    length = len(lst)    
    return health/length, iters/length, avdeg/length, max_clust/length,crash_size/length, comp_size/length

'''
results2018-12-08_1000.pgz
results2018-12-07_2000.pgz
results2018-11-30_4000.pgz
results2018-09-28_8000.pgz
results2018-08-25.pgz
'''


def file_reader(filename):
    with gzip.open(filename,'r') as h:
        a = pickle.load(h)  # 2018-08-22, (a[(True,'er')][56][0][8]/4000).plot()
    er_health, er_iters, er_avdeg, er_max_clusts, er_crsh_sz, er_comp_sz = sum_up_new(a[True])
    return er_health, er_iters, er_avdeg, er_max_clusts, er_crsh_sz, er_comp_sz 


names = ['results2018-12-12_1000Parallel.pgz',
         'results2018-12-12_2000Parallel.pgz',
         'results2018-12-12_4000Parallel.pgz',
         'results2018-12-12_8000Parallel.pgz',
         'results2018-12-12_16000Parallel.pgz',
         'results2018-12-13_32000Parallel.pgz']
max_1 = []
max_2 = []
max_3 = []
plt.figure()
ratios = [1,2,4,8,16,32]
alpha = 0.1326325475290079
line = lambda x: x*0.1326325475290079 + 0.7894845636178769
for i,fil in enumerate(names):
    er_health, er_iters, er_avdeg, er_max_clusts, er_crsh_sz, er_comp_sz = file_reader(fil)
    (((1.0/ratios[i])**alpha)*er_iters[8]).plot()
    max_1.append(er_iters[8].max())
    max_2.append(er_iters[8][1-1/7])
    max_3.append(er_iters[8][1-1/7:1-1/7.75].mean())
    

with gzip.open('results2018-08-25.pgz','r') as h:
    a = pickle.load(h)  # 2018-08-22, (a[(True,'er')][56][0][8]/4000).plot()
# items before 2018-08-20 are lists of lists of dicts or something. use sum_up_old
er_health, er_iters, er_avdeg, er_max_clusts, er_crsh_sz, er_comp_sz = sum_up_new(a[(True,'er')])
sf_health, sf_iters, sf_avdeg, sf_max_clusts, sf_crsh_sz, sf_comp_sz = sum_up_new(a[True,'powpow'])
hh_health, hh_iters, hh_avdeg, hh_max_clusts, hh_crsh_sz, hh_comp_sz = sum_up_new(a[True,'man_half_half'])


# items from 2018-08-20 and on are dicts of lists



etall=0
efall = 0
ptall=0
pfall=0
mfall=0
mtall=0
for i in range(150):
    et = a[0][i][('er', True)]
#    nhet = et[0]
    avdget = et[2]
#    nhet.plot(title='et')

    pt = a[1][i][('powpow', True)]
#    nhpt = pt[0]
    avdgpt = pt[2]
#    nhpt.plot(title='pt')

    mt = a[2][i][('man_half_half', True)]
#    nhmt = mt[0]
    avdgmt = mt[2]
#    nhmt.plot(title='mt')

    ef = a[3][i][('er', False)]
#    nhef = ef[0]
    avdgef = ef[2]
#    nhef.plot(title='ef')

    pf = a[4][i][('powpow', False)]
#    nhpf = pf[0]
    avdgpf = pf[2]
    #nhpf.plot(title='pf')
    
    mf = a[5][i][('man_half_half', False)]
#    nhmf = mf[0]
    avdgmf = mf[2]
    #nhmf.plot(title='mf')
#plt.figure()
    mfall+=avdgmf[avdgmf>0].min()#.plot(style='.')
#    print avdgmf[avdgmf>0][17].min()
    mtall+=avdgmt[avdgmt>0].min()#.plot(style='.')
#    plt.title('Model');plt.xlabel('Initial average Degree');plt.ylabel('Cascade degree')
#    plt.figure()
    efall+=avdgef[avdgef>0].min()#.plot(style='.')
    etall+=avdget[avdget>0].min()#.plot(style='.')
#    plt.title('E-R');plt.xlabel('Initial average Degree');plt.ylabel('Cascade degree')
#    plt.figure()
    pfall+=avdgpf[avdgpf>0].min()#.plot(style='.')
    ptall+=avdgpt[avdgpt>0].min()#.plot(style='.')
plt.figure();(0.004*efall).plot(style='-o');(0.004*etall).plot(style='-.')
plt.figure();(0.004*pfall).plot(style='-o');(0.004*ptall).plot(style='-.')

plt.figure();(0.004*mfall).plot(style='-o');(0.004*mtall).plot(style='-.')
plt.title('Scale free');plt.xlabel('Initial average Degree');plt.ylabel('Cascade degree')