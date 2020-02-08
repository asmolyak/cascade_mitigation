# -*- coding: utf-8 -*-
"""
Created on Thu Feb 02 23:56:09 2017

@author: Alex
"""

#from multiprocessing import Pool, TimeoutError
from functools import partial
from scipy.special import binom
from scipy.stats import poisson, pareto
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import pandas as pd

def calc_min_holdings(average, power):
    return average*(power-1)/power


def draw(X,Y,Z,xlabel,ylabel,zlabel,title):
    fig = plt.figure()
    ax = fig.gca(projection='3d')

#    X=hit_range
#    Y = health_threshold
    X,Y = np.meshgrid(X,Y)
    surf = ax.plot_surface(X,Y,Z.as_matrix().T,cmap=cm.coolwarm,rstride=1,cstride=1,linewidth=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.set_title(title)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

def r(k, j, threshold):
    if j <= threshold*k:
        return 0 
    else:
        return 1


def calc_G(k_max, p_x, f, th, k_min=0):
    G = 0
    for k in range(k_min, k_max+1):
#        p = p_x(x=k+0.49)-p_x(x=k-0.5)
        p = p_x(k)
        for j in range(k_min, k+1):
            # print p,r(k,j,th),f,j,(1-f),k-j,binomkj[k,j]
            G += p * r(k, j, th) * (f**j) * ((1-f)**(k-j)) * binomkj[k, j]
            if p * r(k, j, th) * (f**j) * ((1-f)**(k-j)) * binomkj[k, j] < 0:
                print p, r(k, j, th), f, j, (1-f), k-j, binomkj[k, j]
    return G


def calc_H(k_max, k_mean, p_x, f, th, k_min=0):
    H = 0
    for k in range(k_min,k_max+1):
#        p = p_x(x=k+0.49)-p_x(x=k-0.5) # for scake free
        p = p_x(k)
        for j in range(k_min,k+1):
            H+= p * r(k,j,th) * j * (f**(j-1)) * ((1-f)**(k-j)) * binomkj[k,j] / k_mean
    return H


def calc_inner_loop(k_maxA, k_maxB, thA, thB, pA, pB, muA, k_min_A=0, k_min_B=0):
    muAcont = np.zeros(num_of_inner_iters)
    muBcont = np.zeros(num_of_inner_iters)
    for n in range(num_of_inner_iters):
        muB1 = calc_G(k_max=k_maxB, th=thA, p_x=pB, f=muA, k_min=k_min_B)
        muB = calc_H(k_max=k_maxB, th=thA, p_x=pB, f=muA, k_mean=k_meanB, k_min=k_min_B)
        muA = p*calc_H(k_max=k_maxA, th=thB, p_x=pA, f=muB, k_mean=k_meanA, k_min=k_min_A)
        muA1 = p*calc_G(k_max=k_maxA, th=thB, p_x=pA, f=muB, k_min=k_min_A)
        muAcont[n] = muA1
        muBcont[n] = muB1
        if muA1 != muA1:
            break
    return muAcont, muBcont
        
        
k_maxA = 50
k_maxB = 50
k_min_A = 0
k_min_B = 0
k_bin = max(k_maxA,k_maxB)

num_of_inner_iters = 200

binomkj = np.array(np.zeros([k_bin+1, k_bin+1]))
for k in range(k_bin+1):
    for j in range(k+1):
        binomkj[k,j] = binom(k, j)
          
p_range = [1-0.00025]  # np.arange(0.999, .9999, 0.0001)
threshold = 1-1.0/np.arange(4.0, 17.0, 0.25)  # np.arange(0.73,0.94,0.001)
  
agg = {}
agg_all = {}
# k_meanA = 3#12.7333
# k_meanB = 3#35.8125
# pA = partial(pareto.pdf,b=1.5)
for k_meanA in np.arange(5, 20):
    k_meanB = k_meanA
#    k_maxA = int(max(15, 1000))
#    k_maxB = int(max(15, 1000))
#    scale_a = calc_min_holdings(k_meanA,2.5)
#    scale_b = calc_min_holdings(k_meanB,2.5)
    pA = partial(poisson.pmf, mu=k_meanA)
    pB = partial(poisson.pmf, mu=k_meanB)
#    pA = partial(pareto.pdf, b=2.5, loc=0, scale=scale_a)
#    pB = partial(pareto.pdf, b=2.5, loc=0, scale=scale_b)
#    cA = partial(pareto.cdf, b=2.5, loc=0, scale=scale_a)
#    cB = partial(pareto.cdf, b=2.5, loc=0, scale=scale_b)

    x, y, z = 0, 0, 0
    muMat = np.ndarray(shape=(len(p_range), len(threshold), num_of_inner_iters))
    muBMat = np.ndarray(shape=(len(p_range), len(threshold), num_of_inner_iters))
    for p in p_range:  # [0.999]:#
        y = 0
        muA = p
        part_calc = partial(calc_inner_loop, k_maxA=k_maxA, k_maxB=k_maxB, pA=pA, pB=pB, muA=muA)
        for th in threshold:
            z = 0
            thA = th
            thB = th
            muAint, muBint = calc_inner_loop(k_maxA, k_maxB, thA, thB, pA, pB, muA, k_min_A=k_min_A, k_min_B=k_min_B)
            muMat[x, y, :] = muAint
            muBMat[x, y, :] = muBint
            y += 1
        x += 1
    out = np.array(np.zeros([len(p_range), len(threshold)]))
    outB = np.array(np.zeros([len(p_range), len(threshold)]))
    for i in range(len(p_range)):
        for j in range(len(threshold)):
            out[i, j] = muMat[i, j, -1]
            outB[i, j] = muBMat[i, j, -1]
    # generate fig 1 for paper, -1 for infinitesimal impact
    plt.plot(threshold,out[-1, :])
    agg_all[k_meanA] = out[-1, :]
    agg[k_meanA] = threshold[np.argmin(np.diff(out[-1, :]))]
    # draw(np.arange(1,0.9801,-0.0001),np.arange(0,1,0.01), pd.DataFrame(out),
    # 'impact','threshold','Asset health','Analytic solution')
    # draw(np.arange(1,0.9801,-0.0001),np.arange(0,1,0.01), 1000*(pd.DataFrame(outB+out)),
    # 'impact','threshold','Network health','Network Health')
