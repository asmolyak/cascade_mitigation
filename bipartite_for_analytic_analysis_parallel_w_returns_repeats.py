#!/usr/bin/env python2
# -*- coding: utf-8 -*-


from __future__ import division
from scipy.stats import powerlaw
from scipy import optimize
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import random
from matplotlib import cm
import pandas as pd
from collections import defaultdict
import cPickle as pickle
import gzip
from datetime import date
from itertools import product
import copy


def fragility(G, threshold, avg):
    """
    this function calculates a possible fragility thingamagig while looking at three parameters.
    first, "cent", having nonzero amount of fragile neighbors.
    second, "fragile", is 1 if the node itself is fragile and 0 otherwise
    third, "percolating" is the product of the two, i.e. it's nonzero for nodes that are both fragile and have fragile neighbors.
    """
    if threshold < 1:
        th = 1.0/(1-threshold)
    else:
        th = threshold
    test = {j: sum((G.degree(i) < th for i in G.neighbors(j))) for j in G.nodes()}
    test2 = {j: 0 if G.degree(j) > th+1 else 1 for j in G.nodes()}
    perc_nodes = nx.Graph()
    for node, data in G.nodes(data=True):
        data['cent'] = test[node]
        data['fragile'] = test2[node]
        data['percolating'] = test2[node] * test[node]
        if data['percolating'] > 1:
            perc_nodes.add_node(node)
    frag = nx.Graph()
    perc = nx.Graph()
    cent = nx.Graph()
    for e in G.edges():
        if G.node[e[0]]['fragile'] > 1 and G.node[e[1]]['fragile'] > 1:
            frag.add_edge(*e)
        if G.node[e[0]]['cent'] > 0 and G.node[e[1]]['cent'] > 0:
            cent.add_edge(*e)
        if G.node[e[0]]['fragile'] > 0 and G.node[e[1]]['fragile'] > 0 and \
                        G.node[e[0]]['cent'] > 1 and G.node[e[1]]['cent'] > 1 :
            perc.add_edge(*e)
    safe_perc = nx.Graph()
    safe_perc_val = nx.Graph()
    
    av = min((np.mean(dict(perc.degree()).values()))**0.5, 0.75*avg)
    for e in perc.edges():
        if (perc.degree(e[0]) > av-1 and perc.degree(e[1]) > av) or (perc.degree(e[0]) > av and perc.degree(e[1]) > av-1):
            safe_perc.add_edge(*e)
        if (G.node[e[0]]['percolating'] > 2 and G.node[e[1]]['percolating'] > 2) or \
            (G.node[e[0]]['percolating'] > 2 and G.node[e[1]]['percolating'] > 2):
            safe_perc_val.add_edge(*e)
#    if threshold>1:
#        print th, np.mean(G.degree().values()),np.mean(frag.degree().values()), np.mean(cent.degree().values()),  np.mean(perc.degree().values()), np.mean(safe_perc.degree().values())
    return frag, cent, perc, safe_perc,safe_perc_val,perc_nodes


def draww(X, Y, Z, xlabel, ylabel, zlabel, title):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
#    X=hit_range
#    Y = health_threshold
    X,Y = np.meshgrid(X,Y)
    surf = ax.plot_surface(X, Y, Z.as_matrix().T, cmap=cm.coolwarm, rstride=1, cstride=1, linewidth=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.set_title(title)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()


def average_degree(graph,size=None):
    if size:
        return sum([j for (i, j) in graph.degree()])/float(size)
    else:
        return sum([j for (i,j) in graph.degree()])/float(len(graph))


def calc_min_holdings(average, power):
    return average*(power-1)/power


def generate_distribution(p, min_holdings, num_of_banks, num_of_assets, fit=False, fixed_mean=False, average=16+2.0/3):
    if fixed_mean:
        min_holdings = calc_min_holdings(average,p)
#    np.random.seed(1)
    while True:
        rvss = 1/powerlaw.rvs(p, scale=1/min_holdings, size=num_of_banks)
        if max(rvss) < num_of_assets:
            break
    fixed_rvs = rvss.copy()
    round_rvs = map(round, fixed_rvs)
    if fit:
        hist, bins = np.histogram(rvss, bins=np.logspace(np.log10(min(rvss)), np.log10(max(rvss))))
        fitfunc = lambda p, x: p[0]*x**p[1]
        errfunc = lambda p, x, y: fitfunc(p, x) - y
        p0theo = (p-1)*(min_holdings**(p-1))
        p0 = [0.5, -3]
        p1, success = optimize.leastsq(errfunc,p0[:], args=(bins[:-1], hist))
#        plt.plot(bins[:-1],hist,'.',hold=True)
#        plt.hold()
#        plt.plot(bins[:-1], fitfunc(p1,bins[:-1]), hold=True)
        print p1, p0theo, np.mean(rvss)
    return round_rvs


def setup_network(num_of_banks, num_of_assets, round_rvs, kind='man_half_half', p=2.5, min_holdings=10,
                  single_net=False, project_network=False, av_deg=False, add_fedfunds=None):
    if not single_net:
        if kind == 'er':
            # np.random.seed(1)
            if av_deg:
                G = nx.bipartite.random_graph(num_of_banks, num_of_assets, av_deg/num_of_banks)  # ,seed=1
            else:
                G = nx.bipartite.random_graph(num_of_banks, num_of_assets, np.mean(round_rvs)/num_of_banks)  # ,seed=1
            if not nx.is_connected(G):
                G = max(nx.connected_component_subgraphs(G), key=len)
            bot, top = nx.bipartite.sets(G)
            assets = {t: 'a' + str(t-num_of_banks) for t in top}
            banks = {b: 'b' + str(b) for b in bot}
            z = assets.copy()
            z.update(banks)
            nx.relabel_nodes(G, z, copy=False)
            for node, data in G.nodes(data=True):
                if node.startswith('b'):
                    data['kind'] = 'bank'
                elif node.startswith('a'):
                    data['kind'] = 'asset'
#            print 'random: ', average_degree(G), len(G)
        elif kind == 'powpow':
            G = nx.bipartite.configuration_model(map(int, round_rvs), map(int, round_rvs), create_using=nx.Graph())
            bot, top = nx.bipartite.sets(G)
            assets = {t: 'a' + str(t-num_of_banks) for t in top}
            banks = {b: 'b' + str(b) for b in bot}
            z = assets.copy()
            z.update(banks)
            nx.relabel_nodes(G, z, copy=False)
            for node,data in G.nodes(data=True):
                if node.startswith('b'):
                    data['kind'] = 'bank'
                else:
                    data['kind'] = 'asset'
#            print 'powpow: ', average_degree(G)
        elif kind == 'man_half_half':
            banks = ['b'+str(i) for i in range(num_of_banks)]
            assets = ['a'+str(i) for i in range(num_of_assets)]
            G = nx.Graph()
            G.add_nodes_from(banks,kind='bank')
            G.add_nodes_from(assets,kind='asset')
#            np.random.seed(1)
            for i,j in enumerate(round_rvs):
                G.add_edges_from([('b'+str(i),'a'+str(k)) for k in np.random.choice(num_of_assets,int(j),replace=False)])
#            print 'power law: ', average_degree(G)
        elif kind == 'auto_half_and_half':
            temp_rvs = round_rvs[:-1]     
#            np.random.seed(1)
            while True:            
                pois = np.random.poisson(np.mean(round_rvs),num_of_assets)
                pd.Series(pois).hist()
                rem = sum(pois)-sum(temp_rvs)
                if rem < 0:
                    continue
                else:
                    if ((p-1)*(min_holdings**(p-1)))*(rem**(-p)) > 0.0001:
                        round_rvs[-1] = rem
                        break
                    else:
                        continue
            G = nx.bipartite.configuration_model(map(int, round_rvs), map(int, pois), create_using=nx.Graph())
            bot,top = nx.bipartite.sets(G)
            assets = {t: 'a' + str(t-num_of_banks) for t in top}
            banks = {b: 'b' + str(b) for b in bot}
            z = assets.copy()
            z.update(banks)
            nx.relabel_nodes(G,z,copy=False)
            for node,data in G.nodes(data=True):
                if node.startswith('b'):
                    data['kind'] = 'bank'
                else:
                    data['kind'] = 'asset'
#            print 'auto_h_h: ', average_degree(G)
        if project_network:
            G = nx.bipartite.project(G,[node for node,data in G.nodes(data=True) if data['kind']=='bank'])
            for node,data in G.nodes(data=True):
                data['init_deg'] = G.degree(node)
        else:
            if add_fedfunds:
                G.add_node('f01')
                new_edges = [('f01', 'b'+str(i)) for i in random.sample(range(num_of_banks), add_fedfunds)]
                G.add_edges_from(new_edges)
                G.node['f01']['kind'] = 'fed'

    else:
        if kind == 'er':
            # np.random.seed(1)
            G = nx.fast_gnp_random_graph(num_of_banks, np.mean(round_rvs)/num_of_banks)
        elif kind == 'man_half_half':
            if not nx.is_valid_degree_sequence(map(int,round_rvs)):
                round_rvs[0] += 1
            G = nx.configuration_model(map(int, round_rvs), create_using=nx.Graph())
        banks = {b: 'b' + str(b) for b in range(num_of_banks)}
        nx.relabel_nodes(G, banks, copy=False)
        for node, data in G.nodes(data=True):
            data['kind'] = 'bank'
    for node,data in G.nodes(data=True):
        data['init_deg'] = G.degree(node)
    return G


def run_core(params, hit_probability=0.000125):
    # fig = plt.figure()
    # gather_gc_by_avdeg = fig.gca()
    # fig2 = plt.figure()
    # gather_iters_by_avdeg = fig2.gca()
    fract_threshold, kind = params[0], params[1]
    num_of_banks = 4000
    num_of_assets = 4000
    if fract_threshold:
        health_threshold = 1-1.0/np.arange(4.0, 17.0, 0.25)  # np.arange(0.73,0.94,0.001)#
    else:
        health_threshold = np.arange(4.0, 17.0, 0.25)  #, th_base/((1+th_base)**fract_threshold)#np.arange(0,1,0.01)#
#        hit_range = [0.001]  # np.arange(0,.01,0.0001)
    print params, health_threshold
    av_deg_range = np.arange(6.0, 10.0, 0.5)  # np.logspace(np.log10(0.01),np.log10(18),100)
    power = 2.5
    min_holdings = 3
    res=[]
    for i in range(200):
        network_health = pd.DataFrame(columns=av_deg_range, index=health_threshold)
        iterations = network_health.copy()
        degree = network_health.copy()
        max_clust = network_health.copy()
        core_num = network_health.copy()
        bank_nodes = network_health.copy()
        asset_ndoes = network_health.copy()
        perc_frame = network_health.copy()
        safe_perc_frame = network_health.copy()
        safe_perc_val_frame = network_health.copy()
        cent_avdeg = network_health.copy()
        perc_avdeg = network_health.copy()
        safe_perc_avdeg = network_health.copy()
        safe_perc_bank = network_health.copy()
        safe_perc_asset = network_health.copy()
        perc_node_asset = network_health.copy()
    
        fail_counter = pd.DataFrame()
        size_counter = pd.DataFrame()
        for avg in av_deg_range:
            round_rvs = generate_distribution(power,min_holdings,num_of_banks=num_of_banks,\
                                              num_of_assets=num_of_assets,fit=False,fixed_mean=True,average=avg)
            single_net = False
            project_network = False
            link_temp = 0
        
            G = setup_network(num_of_banks=num_of_banks,num_of_assets=num_of_assets,\
                              round_rvs=round_rvs,p=power,kind=kind,min_holdings=min_holdings,\
                              single_net=single_net,project_network=project_network,av_deg=avg)
            output = {}
            G_copy = G.copy()
            for th in health_threshold:
                local_counter = []
                local_size = [len(G)]
                G = G_copy.copy()
#                random.seed(1)
                if project_network or single_net:
                    init_nodes_to_remove = ['b'+str(i) for i in random.sample(range(num_of_assets), 1)] #int(hit_probability*num_of_assets)
                else:
                    init_nodes_to_remove = ['a'+str(i) for i in random.sample(range(num_of_assets), 1)] # int(hit_probability*num_of_assets)
                link_temp = sum(dict(G.degree(init_nodes_to_remove)).values())
                G.remove_nodes_from(init_nodes_to_remove)
#                _, cent, perc, safe_perc,safe_perc_val,perc_nodes = fragility(G,threshold=th,avg=avg)
                local_counter.append(len(init_nodes_to_remove))
                counter = 0
                while True:
                    counter += 1
                    b_nodes_to_remove = set([node for node,data in G.nodes(data=True) if
                                             data['kind'] == 'bank' and
                                             ((data['init_deg'] == 0) or
                                              (1.0*G.degree(node)/(data['init_deg']**fract_threshold) <= th) or
                                                 ((not fract_threshold)*(G.degree(node) <= (th+1))*(np.random.rand() < (th%1))))])#and (node not in safe_perc)
                    if len(b_nodes_to_remove) == 0:
                        break
                    local_counter.append(len(b_nodes_to_remove))
                    link_temp += sum(dict(G.degree(b_nodes_to_remove)).values())
                    G.remove_nodes_from(b_nodes_to_remove)
                    try:
                        sz = len(max(nx.connected_components(G), key=len))
                    except:
                        sz = 0
                    local_size.append(sz)
    
                    if not project_network:
                        a_nodes_to_remove = set([node for node,data in G.nodes(data=True) if
                                                 data['kind']=='asset' and
                                                 ((data['init_deg']==0) or
                                                  (1.0*G.degree(node)/(data['init_deg']**fract_threshold) <= th) or
                                                     ((not fract_threshold)*(G.degree(node) <= (th+1))*(np.random.rand() < (th%1))))])#
                        if len(a_nodes_to_remove) == 0:
                            break
                        link_temp += sum(dict(G.degree(a_nodes_to_remove)).values())
                        G.remove_nodes_from(a_nodes_to_remove)
                try:
                    temp = average_degree(G, num_of_assets+num_of_banks)
                    temp2 = len(max(nx.connected_components(G), key=len))/(num_of_assets+num_of_banks)  # len(G)
                    temp3 = 0  # max(nx.core_number(G).values())
                except:
                    temp = 0
                    temp2 = 0
                    temp3 = 0
                temp4 = 0 # len([node for node, data in G.nodes(data=True) if data['kind'] == 'bank'])
                temp5 = 0 # len([node for node, data in G.nodes(data=True) if data['kind'] == 'asset'])
                temp6 = 0  # len(perc) if not np.isnan(len(perc)) else 0
                temp7 = 0  # len(safe_perc) if not np.isnan(len(safe_perc)) else 0
                temp8 = 0  # len(safe_perc_val) if not np.isnan(len(safe_perc_val)) else 0
                temp9 = 0  # np.mean(dict(cent.degree()).values())
                temp10 = 0  # np.mean(dict(perc.degree()).values())
                temp11 = 0  # np.mean(dict(safe_perc.degree()).values())
                temp12 = 0  # len([node for node in safe_perc.nodes() if node.startswith('b')])
                temp13 = 0  # len([node for node in safe_perc.nodes() if node.startswith('a')])
                temp14 = 0  # len([node for node in perc_nodes.nodes() if node.startswith('a')])
                # I removed the int() around avg because of fractional steps
                fail_counter = pd.concat([fail_counter,pd.Series(local_counter,name=(avg,th))],axis=1)
                size_counter = pd.concat([size_counter,pd.Series(local_size,name=(avg,th))],axis=1)
                output[th] = (counter,len(G),temp,temp2,temp3,temp4,temp5)
                # ,temp6,temp7,temp8,temp9,temp10,temp11,temp12,temp13,temp14
#            print 'done iterating'
            df = pd.DataFrame(output).T
            network_health[avg] = df[1]
            iterations[avg] = df[0]
            degree[avg] = df[2]
            max_clust[avg] = df[3]
            core_num[avg] = df[4]
            bank_nodes[avg] = df[5]
            asset_ndoes[avg] = df[6]
            # perc_frame[avg] = df[7]
            # safe_perc_frame[avg] = df[8]
            # safe_perc_val_frame[avg] = df[9]
            # cent_avdeg[avg] = df[10]
            # perc_avdeg[avg] = df[11]
            # safe_perc_avdeg[avg] = df[12]
            # safe_perc_bank[avg] = df[13]
            # safe_perc_asset[avg] = df[14]
            # perc_node_asset[avg] = df[15]
        res.append((network_health, iterations, degree, max_clust, fail_counter, size_counter))
        # , perc_frame, safe_perc_frame, safe_perc_val_frame, perc_node_asset
#    print 'done building df'
    print 'all done', params
    return params, res


if __name__ == "__main__":
    from multiprocessing.pool import Pool
    print "let's give it a go"
    fract_threshold = [True]#, False
    kind = ['er']  #, 'man_half_half', 'powpow'
    params = product(fract_threshold, kind)  # [(u,v) for u in fract_threshold for v in kind]
#    pool = Pool(6)
    results = {i:j for i,j in map(run_core, params)} # pool.
#    tmp = results[0][0][('er',False)][0]
#    tmp[tmp>0].min().plot()
#    pool.close()
#    pool.join()
#    rd = {}
    with gzip.GzipFile('results'+str(date.today())+'_8000Single.pgz', 'w') as h:
        pickle.dump(results, h, protocol=pickle.HIGHEST_PROTOCOL)
    print 'done with pool'
#    for i in results:
#        for k,v in i.iteritems():
#            print k
#            x=[]
#            y=[]
#            for col in v[2].columns:
#                x.append(col)
#                y.append( v[2][col][v[2][col]>0.1].min())
#            plt.plot(x,y)
