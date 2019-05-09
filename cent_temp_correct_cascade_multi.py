import networkx as nx
import numpy as np
import pandas as pd
from multiprocessing import Pool
from itertools import product, islice
from functools import partial
from operator import itemgetter
from scipy import stats as sts
import matplotlib.pyplot as plt
from parallel_cent_test import betweenness_centrality_parallel
from bipartite_for_analytic_analysis_parallel_w_returns_repeats import generate_distribution, setup_network


def sf_graph(num,power=2.5,avdeg=3.0):
    if avdeg:
        xmin = avdeg*(power-2)/(power-1)
        if xmin < 1:
            xmin = 1
    else:
        xmin=1
    rvs = sts.pareto.rvs(b=power-1,scale=xmin,size=num)
    rrvs = map(int, map(np.round, rvs))
    if sum(rrvs)%2:#not nx.is_valid_degree_sequence_erdos_gallai(rrvs):
        rrvs[0] += 1
    g = nx.configuration_model(rrvs,create_using=nx.Graph())
    return g


def set_threshold(G, threshold, kind, mean=0.5, std=0.2):
    if kind == 'simple' or kind == 'kcore' or  kind == 'articulation' or kind == 'core_number':
        return G
    elif kind == 'random':
        for n,d in G.nodes(data=True):
            d['th'] = mean+std*np.random.randn()
        return G


def clean_orphans(g,h=None):
    orphans = []
    for node,deg in g.degree():
        if deg == 0:
            orphans.append(node)
    if h:
        for node,deg in h.degree():
            if deg == 0:
                orphans.append(node)
    g.remove_nodes_from(orphans)
    if h:
        h.remove_nodes_from(orphans)
    if h:
        return g,h
    return g


def set_vulnerability(g, kind, threshold, factor=1.0):
    if kind == 'simple':
        test = {j: sum((g.degree(i) < factor/(1-threshold) for i in g.neighbors(j))) for j in g.nodes()}
        test2 = {j: 1 if (g.degree(j) < factor/(1-threshold) or threshold>1) else 0 for j in g.nodes()}
    elif kind == 'random':
        test = {j: sum((g.degree(i) < 1/(1-g.node[i]['th']) for i in g.neighbors(j))) for j in g.nodes()}
        test2 = {j: 1 if (g.degree(j) < 1/(1-g.node[j]['th']) or g.node[j]['th'] > 1) else 0 for j in g.nodes()}
    elif kind == 'kcore':
        test = {j: sum((g.degree(i) <= threshold+1 for i in g.neighbors(j))) for j in g.nodes()}
        test2 = {j: 1 if (g.degree(j) <= threshold+1) else 0 for j in g.nodes()}
    elif kind == 'articulation':
        all_max_nodes = max(nx.connected_components(g),key = len)
        all_max_graph = max(nx.connected_component_subgraphs(g),key = len)
        test2 = dict.fromkeys(g.nodes,False)
        test = dict.fromkeys(g.nodes,0)
        art_pt = nx.articulation_points(all_max_graph)
        for node in art_pt:
            temp = all_max_graph.copy()
            temp.remove_node(node)
            max_com_temp = max(nx.connected_components(temp),key = len)
            dlen = len(all_max_nodes) - len(max_com_temp)
            dnodes = set(all_max_nodes) - set(max_com_temp) - set([node])
            test[node] = dlen
            if dlen > 2:
                for k in dnodes:
                    test2[k] = True 
    elif kind == 'core_number':
        cn = nx.core_number(g)
        test = {j: sum((cn[i] < 3 for i in g.neighbors(j))) for j in g.nodes()}
        test2 = {j: 1 if (cn[j] < 3) else 0 for j in g.nodes()}
    seed_size = []
    unprotected_seed = []
    for node, data in g.nodes(data=True):
        data['init_deg'] = g.degree(node)
        data['cent'] = test[node]
        data['fragile'] = test2[node]
        data['percolating'] = test2[node] * test[node]
        if kind == 'random' and data['th'] > 1:
            seed_size.append(node)
            if data['percolating'] == 0:
                unprotected_seed.append(node)
    return g, seed_size, unprotected_seed


def get_protected_nodes(g, threshold, method, degree_min=2, avg=None):
    broadg = nx.subgraph(g, [n for n, d in g.nodes(True) if d['percolating'] >= degree_min])
    if not avg:
        avg = np.mean(dict(broadg.degree()).values())
    safe_perc = [n for n, d in broadg.degree() if d>(avg-1)]
#    print 'thresh: {0:.3f}, perc: {1:3d}, avdeg: {2:.2f}, safe_prec: {3:3d}'.format(threshold,len(broadg),avg,len(safe_perc)) 
    ln = len(safe_perc)
    lst = sorted(dict(g.degree()).iteritems(), key=itemgetter(1), reverse=True)
    if method == 'top_degree':
        safe_perc = zip(*lst)[0][:ln]
    elif method == 'bottom_deg':
        safe_perc = zip(*lst)[0][-ln:]
    elif method == 'random':
        safe_perc = np.random.choice(g,ln,replace=False)
    elif method == 'eigen_cent':
        centrality = nx.eigenvector_centrality_numpy(g)
        lst = sorted(centrality.iteritems(), key=itemgetter(1), reverse=True)
        safe_perc = zip(*lst)[0][:ln]
    elif method == 'between_cent':
        centrality = nx.betweenness_centrality(g,1000)#
        lst = sorted(centrality.iteritems(), key=itemgetter(1), reverse=True)
        safe_perc = zip(*lst)[0][:ln]
    elif method == 'degree':
        centrality = nx.degree_centrality(g)
        lst = sorted(centrality.iteritems(), key=itemgetter(1), reverse=True)
        safe_perc = zip(*lst)[0][:ln]
    elif method == 'closeness':
        centrality = nx.closeness_centrality(g)
        lst = sorted(centrality.iteritems(), key=itemgetter(1), reverse=True)
        safe_perc = zip(*lst)[0][:ln]
    elif method == 'katz':
        centrality = nx.katz_centrality_numpy(g)
        lst = sorted(centrality.iteritems(), key=itemgetter(1), reverse=True)
        safe_perc = zip(*lst)[0][:ln]
    return safe_perc


def cascade(g, safe_perc, kind, dont_save, threshold, num_of_init_nodes=1, h=None, safe_perch=None):
    remove = np.random.choice(g,num_of_init_nodes,replace=False)
    nodes_to_test = set([i for j in remove for i in g.neighbors(j)])
    if h:
        nodes_to_test_h = set([i for j in remove for i in h.neighbors(j)])
    g.remove_nodes_from(remove)
    if h:
        h.remove_nodes_from(remove)
    if kind == 'random':
        tested_nodes = set()
        while len(nodes_to_test)>0:
            node_to_test = nodes_to_test.pop()
            if node_to_test in g:
                if (node_to_test not in safe_perc) and (((g.degree(node_to_test) / float(g.nodes[node_to_test]['init_deg'])) < g.nodes[node_to_test]['th'])):
                    g.remove_node(node_to_test)
                elif (g.degree(node_to_test) / float(g.nodes[node_to_test]['init_deg']) < g.nodes[node_to_test]['th']) and (g.nodes[node_to_test] in safe_perc and np.random.rand()<dont_save):
                    g.remove_node(node_to_test)
    elif kind == 'simple':
        tested_nodes = set()
        while len(nodes_to_test)>0:
            node_to_test = nodes_to_test.pop()
            if node_to_test in g:
                if (node_to_test not in safe_perc) and ((g.degree(node_to_test)/float(g.nodes[node_to_test]['init_deg'])) < threshold):
                    tested_nodes |= set([node_to_test])
                    nodes_to_test |= set([n for n in g.neighbors(node_to_test) if n not in tested_nodes])
                    g.remove_node(node_to_test)
                elif (g.degree(node_to_test) / float(g.nodes[node_to_test]['init_deg']) < threshold) and ((node_to_test in safe_perc) and (np.random.rand() < dont_save)):
                    tested_nodes |= set([node_to_test])
                    nodes_to_test |= set([n for n in g.neighbors(node_to_test) if n not in tested_nodes])
                    g.remove_node(node_to_test)
    elif kind == 'kcore':
        tested_nodes = set()
        while len(nodes_to_test)>0:
            node_to_test = nodes_to_test.pop()
            if node_to_test in g:
                if (node_to_test not in safe_perc) and ((g.degree(node_to_test) < threshold) ):#or (g.degree(node_to_test) < (th+1) and np.random.rand()<(th%1))
                    tested_nodes |= set([node_to_test])
                    nodes_to_test |= set([n for n in g.neighbors(node_to_test) if n not in tested_nodes])
                    g.remove_node(node_to_test)
                elif (g.degree(node_to_test) < threshold) and ((node_to_test in safe_perc) and (np.random.rand() < dont_save)):
                    tested_nodes |= set([node_to_test])
                    nodes_to_test |= set([n for n in g.neighbors(node_to_test) if n not in tested_nodes])
                    g.remove_node(node_to_test)
    elif kind == 'articulation' or kind == 'core_number':
        max_clust_a = max(nx.connected_components(g),key = len)
        max_clust_b = max(nx.connected_components(h),key = len)
        conn_nodes_to_remove_a = []
        conn_nodes_to_remove_b = []
        for node_a,node_b in zip(g.nodes(),h.nodes()):
            if node_a not in max_clust_a:
                conn_nodes_to_remove_a+=list(g.neighbors(node_a))+[node_a]
                nodes_to_test_h|= set(h.neighbors(node_a))
            if node_b not in max_clust_b:
                conn_nodes_to_remove_b+=list(h.neighbors(node_b))+[node_b]
                nodes_to_test|= set(g.neighbors(node_b))
        conn_nodes_to_remove_a = list(set(conn_nodes_to_remove_a))
        conn_nodes_to_remove_b = list(set(conn_nodes_to_remove_b))
        g.remove_nodes_from(conn_nodes_to_remove_a)
        h.remove_nodes_from(conn_nodes_to_remove_a)
        
        h.remove_nodes_from(conn_nodes_to_remove_b)
        g.remove_nodes_from(conn_nodes_to_remove_b)
        nodes_to_test = nodes_to_test - set(remove)-set(conn_nodes_to_remove_a)
        nodes_to_test_h = nodes_to_test_h - set(remove)-set(conn_nodes_to_remove_b)
        while (len(nodes_to_test) or len(nodes_to_test_h)):
            #treat disconnected nodes in a, treat dependency nodes in b
            while len(nodes_to_test) > 0:
                max_clust_a = max(nx.connected_components(g),key = len)     
                node_a = nodes_to_test.pop()
                if (node_a in g) and ((node_a not in safe_perc) and (np.random.rand() < dont_save)) and (node_a not in max_clust_a):
                    nodes_to_test |= set(g.neighbors(node_a))
                    nodes_to_test_h |= set(h.neighbors(node_a))
                    try:
                        nodes_to_test_h.remove(node_a)
                    except:
                        pass
                    g.remove_node(node_a)
                    h.remove_node(node_a)
    
            #treat disconnected nodes in a, treat dependency ndoes in a
            while len(nodes_to_test_h) > 0:
                max_clust_b = max(nx.connected_components(h),key = len)            
                node_b = nodes_to_test_h.pop()
                if (node_b in h) and ((node_b not in safe_perch) and (np.random.rand() < dont_save)) and (node_b not in max_clust_b):
                    nodes_to_test |= set(g.neighbors(node_b))
                    nodes_to_test_h |= set(h.neighbors(node_b))
                    try:
                        nodes_to_test.remove(node_b)
                    except:
                        pass
                    g.remove_node(node_b)
                    h.remove_node(node_b)

    if h:
        return g,h
    return g


avdeg = 5.0
network_size = 5000
kind = 'articulation'

save=1
cont = {}
if kind == 'random':
    threshold = np.arange(0.74,0.8,0.002)
elif kind== 'simple':
    threshold = np.arange(0.8,0.92,0.01)
elif kind=='kcore':
    threshold = np.arange(5,13,0.5)
elif kind=='core_number' or kind=='articulation':
    threshold = np.arange(0.48,0.53,0.005)
network_health = pd.DataFrame(columns=np.arange(0,1.1,0.1),index=threshold)
safe_size = pd.DataFrame(columns=np.arange(0,1.1,0.1),index=threshold)
final_safe_size = pd.DataFrame(columns=np.arange(0,1.1,0.1),index=threshold)

# gc = g.copy()
tally = {}


def inner_core(th,mthd,dont_save):
#            round_rvs = generate_distribution(p=2.5,min_holdings=1,num_of_banks=5000,\
#                                              num_of_assets=5000,fit=False,fixed_mean=True,average=avdeg)
#            single_net = False
#            project_network = False
#        
#            g = setup_network(num_of_banks=5000,num_of_assets=5000,\
#                              round_rvs=round_rvs,p=2.5,kind='man_half_half',min_holdings=1,\
#                              single_net=single_net,project_network=project_network,av_deg=avdeg)
    g = nx.fast_gnp_random_graph(network_size,avdeg/(1.0*network_size))#sf_graph(network_size,power=2.5,avdeg=avdeg)#gc.copy() 
    h = nx.fast_gnp_random_graph(network_size,avdeg/(1.0*network_size))#sf_graph(network_size,power=2.5,avdeg=avdeg)#gc.copy() 
    g,h  = clean_orphans(g,h)
#            h = clean_orphans(h)
    g = set_threshold(g,th,kind,mean=th,std=0.1)
    h = set_threshold(h,th,kind,mean=th,std=0.1)
    g, seed_size, unprotected_seed = set_vulnerability(g,kind,th,factor=1)
    h, seed_size, unprotected_seed = set_vulnerability(h,kind,th,factor=1)
    safe_perc = get_protected_nodes(g,th,mthd)
    safe_perch = get_protected_nodes(h,th,mthd)
    g,h = cascade(g, safe_perc, kind, dont_save, th, num_of_init_nodes=int( th*network_size),h=h,safe_perch=safe_perch)
    return [len(g),len(safe_perc)]


def run(params):
    mthd = params[0]
    dont_save = params[1]
    tally={}
    pool = Pool(7)
    for th in threshold:
        tally[th]=[]
        pf = partial(inner_core, th=th,mthd=mthd,dont_save=dont_save)
        future_res = [pool.apply_async(pf) for _ in range(150)]
        res = [f.get() for f in future_res]
        tally[th] = res
#        for itr in range(150):
##            round_rvs = generate_distribution(p=2.5,min_holdings=1,num_of_banks=5000,\
##                                              num_of_assets=5000,fit=False,fixed_mean=True,average=avdeg)
##            single_net = False
##            project_network = False
##        
##            g = setup_network(num_of_banks=5000,num_of_assets=5000,\
##                              round_rvs=round_rvs,p=2.5,kind='man_half_half',min_holdings=1,\
#            g = nx.fast_gnp_random_graph(network_size,avdeg/(1.0*network_size))#sf_graph(network_size,power=2.5,avdeg=avdeg)#gc.copy() 
#            h = nx.fast_gnp_random_graph(network_size,avdeg/(1.0*network_size))#sf_graph(network_size,power=2.5,avdeg=avdeg)#gc.copy() 
#            g,h  = clean_orphans(g,h)
##            h = clean_orphans(h)
#            g = set_threshold(g,th,kind,mean=th,std=0.1)
#            h = set_threshold(h,th,kind,mean=th,std=0.1)
#            g, seed_size, unprotected_seed = set_vulnerability(g,kind,th,factor=1)
#            h, seed_size, unprotected_seed = set_vulnerability(h,kind,th,factor=1)##                              single_net=single_net,project_network=project_network,av_deg=avdeg)
#            safe_perc = get_protected_nodes(g,th,mthd)
#            safe_perch = get_protected_nodes(h,th,mthd)
#            g,h = cascade(g, safe_perc, kind, dont_save, th, num_of_init_nodes=int(th*network_size),h=h,safe_perch=safe_perch)
#            tally[th].append([len(g),len(safe_perc)])
    pool.close()
    return params,tally
            # print len(g)

if __name__ == '__main__':
#    p = Pool(2)
    protect = ['model']  #,'random','top_degree','bottom_deg' 'degree',,'cent','eigen_cent','closeness','katz','between_cent'
    dontsave = np.arange(0,1.1,0.1)#[0,1]  #0, ,0.3,0.6, 0.3, 0.6, 1
    final = {i: j for i,j in map(run, product(protect, dontsave))}  # p.
#    p.close()
    df_final = pd.DataFrame(columns=[i+'_'+str(np.round(j,1))+k for i, j, k in product(protect,dontsave,['_mean','_min','_max'])], index=final[('model',0)].keys())
    df_count = pd.DataFrame(columns=[i+'_'+str(np.round(j,1))+k for i, j, k in product(protect,dontsave,['_mean','_min','_max'])], index=final[('model',0)].keys())

    for out_key,out_val in final.iteritems():
        for k, v in out_val.iteritems():
            dfkey = out_key[0]+'_'+str(np.round(out_key[1],1))
            df_final[dfkey + '_mean'][k] = np.mean(zip(*v)[0]) / 1.0 / network_size /(1-k)
            df_count[dfkey + '_mean'][k] = np.mean(zip(*v)[1]) / 1.0 / network_size

            temp = []
            for i in range(3):
                temp.append(sum(zip(*v)[0][50 * i:50 * (i + 1)]) / 50.0 / network_size /(1-k))
            df_final[dfkey + '_min'][k] = min(temp) - df_final[dfkey + '_mean'][k]
            df_final[dfkey + '_max'][k] = max(temp) - df_final[dfkey + '_mean'][k]
    df_final.sort_index(inplace=True)
    df_count.sort_index(inplace=True)

    if save:
        df_final.to_csv('cent_temp_final_out_inter_'+str(pd.datetime.now())[:10]+'.csv')
        df_count.to_csv('cent_temp_count_out_inter_'+str(pd.datetime.now())[:10]+'.csv')

    fig = plt.figure(figsize=(16,12))
    plt.tick_params(labelsize=30)

#    plt.plot(df_final.index, df_final[df_final.columns[0]],'.-.',label='Proposed Protection')
#    plt.plot(df_final.index, df_final[df_final.columns[3]],'o-.',label='Random Nodes')
#    plt.plot(df_final.index, df_final[df_final.columns[6]],'^-.',label='Highest Degree')
#    plt.plot(df_final.index, df_final[df_final.columns[9]],'s-.',label='Lowest Degree')
#    
#    for protection strength:
    plt.plot(df_final.index, df_final[df_final.columns[0]],'.-.',label='1.0')
    plt.plot(df_final.index, df_final[df_final.columns[9]],'o-.',label='0.7')
    plt.plot(df_final.index, df_final[df_final.columns[18]],'^-.',label='0.4')
    plt.plot(df_final.index, df_final[df_final.columns[30]],'s-.',label='0.0')
#    
#    plt.errorbar(df_final.index, df_final[df_final.columns[0]], yerr=df_final[[df_final.columns[1], df_final.columns[2]]].abs().T.values,linestyle='-.',marker = '.',capsize=3,label='Proposed Protection')
#    plt.errorbar(df_final.index, df_final[df_final.columns[3]], yerr=df_final[[df_final.columns[4], df_final.columns[5]]].abs().T.values,linestyle='-.',marker = 'o',capsize=3,label='Closeness Centrality')
#    plt.errorbar(df_final.index, df_final[df_final.columns[6]], yerr=df_final[[df_final.columns[7], df_final.columns[8]]].abs().T.values,linestyle='-.',marker = '^',capsize=3,label='Katz Centrality')
#    plt.errorbar(df_final.index, df_final[df_final.columns[9]], yerr=df_final[[df_final.columns[10], df_final.columns[11]]].abs().T.values,linestyle='-.',marker = 's', capsize=3,label='hetweenness Centrality')
#
    plt.plot([0.48]*len(np.arange(0,1.0,0.1)),np.arange(0,1.0,0.1),'--',color='#1f77b4',linewidth=1)
    plt.plot([0.495]*len(np.arange(0,1.0,0.1)),np.arange(0,1.0,0.1),'--',color='#ff7f0e',linewidth=1)
    plt.plot([0.51]*len(np.arange(0,1.0,0.1)),np.arange(0,1.0,0.1),'--',color='#2ca02c',linewidth=1)
    plt.plot([0.53]*len(np.arange(0,1.0,0.1)),np.arange(0,1.0,0.1),'--',color='#d62728',linewidth=1)
    
    plt.legend(loc=0,fontsize=30)
    plt.xlabel('Threshold',fontsize=30)
    plt.ylabel('Probability of System survival',fontsize=30)
#    plt.title('Survival probability for the Model and centrality alternatives')

    axes1 = fig.add_subplot(111)
    axes2 = axes1.twinx()
    df_count['model_0.0_mean'].plot(ax=axes2,style='--r',ylim=[0.0,0.005],label='Protected Nodes')
    plt.ylabel('Fraction of nodes held safe',fontsize=30)
    plt.legend(loc=4,fontsize=30)
    plt.tick_params(labelsize=30)
    plt.show()
    
    
    
    # fig 3b
    fig = plt.figure(6,figsize=(16,12))
    plt.tick_params(labelsize=30)
#    objs = plt.plot(np.arange(0,1.1,0.1),df_final.loc[[0.85,0.865,0.885,0.895]][['model_{}_mean'.format(np.round(i,2)) for i in np.arange(0,1.1,0.1)]].T)
    plt.plot(np.arange(1.1,0,-0.1), df_final.loc[0.48][['model_{}_mean'.format(np.round(i,2)) for i in np.arange(0,1.1,0.1)]].T,'.-',label='0.85')
    plt.plot(np.arange(1.1,0,-0.1), df_final.loc[0.495][['model_{}_mean'.format(np.round(i,2)) for i in np.arange(0,1.1,0.1)]].T,'o-',label='0.865')
    plt.plot(np.arange(1.1,0,-0.1), df_final.loc[0.51][['model_{}_mean'.format(np.round(i,2)) for i in np.arange(0,1.1,0.1)]].T,'^-',label='0.885')
    plt.plot(np.arange(1.1,0,-0.1), df_final.loc[0.53][['model_{}_mean'.format(np.round(i,2)) for i in np.arange(0,1.1,0.1)]].T,'s-',label='0.895')


    plt.xlabel('Probability of Immunization',fontsize=30)
    plt.ylabel('Probability of System survival',fontsize=30)
    plt.legend(loc=0,fontsize=30)
    plt.tick_params(labelsize=30)
    '''
    fig = plt.figure(0)
    plt.errorbar(df_final.index, df_final['model_0_mean'], yerr=df_final[['model_0_min', 'model_0_max']].abs().T.values,
                 linestyle='--', capsize=3, label='Model')
    plt.errorbar(df_final.index, df_final['random_0_mean'],
                 yerr=df_final[['random_0_min', 'random_0_max']].abs().T.values, linestyle='--', capsize=3,
                 label='Random')
    plt.errorbar(df_final.index, df_final['top_degree_0_mean'],
                 yerr=df_final[['top_degree_0_min', 'top_degree_0_max']].abs().T.values, linestyle='--', capsize=3,
                 label='Top Degree')
    plt.errorbar(df_final.index, df_final['bottom_deg_0_mean'],
                 yerr=df_final[['bottom_deg_0_min', 'bottom_deg_0_max']].abs().T.values, linestyle='--', capsize=3,
                 label='Lowest Degree')
    plt.legend()
    plt.show()
    plt.xlabel('Fractional Threshold')
    plt.ylabel('Probability of System survival')
    plt.title('System survival probability for the Model and trivial alternatives')
    '''