from itertools import count
import pandas as pd
import networkx as nx
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Qt5Agg')
import os
import seaborn as sns


def cascade(netowrk, start_node, threshold, synchronous=True):
    if synchronous:
        nodes_to_check = set(netowrk.successors(start_node))  # All nodes case
        netowrk.remove_node(start_node)
        removed_nodes = []
        itr = 0
        while len(nodes_to_check) > 0:
            node_to_test = nodes_to_check.pop()
            if node_to_test in netowrk:
                itr += 1
                current_value = sum([netowrk[i][node_to_test]['weight'] for i in netowrk.predecessors(node_to_test)])
                if (current_value / netowrk.nodes[node_to_test]['money_in']) < (1 - th):
                    nodes_to_check |= set(netowrk.successors(node_to_test))
                    netowrk.remove_node(node_to_test)
                    removed_nodes.append(node_to_test)
        # print(init_node, th, len(the_net))
        net_remaining = len(netowrk)
        num_iters = itr
    else:
        nodes_to_remove = [start_node]
        netowrk.remove_nodes_from(nodes_to_remove)
        c = 0
        while len(nodes_to_remove) > 0:
            nodes_to_remove = [node for node in netowrk.nodes() if
                               'money_in' in netowrk.nodes[node].keys() and
                               (sum([netowrk[i][node]['weight'] for i in netowrk.predecessors(node)]) / netowrk.nodes[node]['money_in']) < (1 - threshold)]
            netowrk.remove_nodes_from(nodes_to_remove)
            c += 1
        print(c, len(netowrk))
        net_remaining = len(netowrk)
        num_iters = c
    return net_remaining, num_iters


def load_data(fname='D:\\Downloads\\Sufang\\ICIO2021_2018.csv', kind='OECD'):
    if kind == 'OECD':
        data = pd.read_csv(fname, index_col=0)
        main_data = data.iloc[:3195, :3195]
    elif kind == 'WIOD':
        tst = pd.read_excel(fname)
        innards = tst.loc[3:, 'Unnamed: 2':].copy()
        numbers = innards.iloc[2:2466, 2:2466].copy()
        row_col = ['_'.join(i) for i in zip(*innards.iloc[:2, 2:2466].values)]
        pd.DataFrame(numbers.values, index=row_col, columns=row_col)
        main_data = pd.DataFrame(numbers.values, index=row_col, columns=row_col)
    else:
        return 'Not Supported'
    main_data[main_data < 0.1] = 0
    main_data = (main_data - main_data*np.eye(main_data.shape[0]))
    return main_data

net_remaining = {}
num_iters = {}

# path = 'D:\\Downloads\\Sufang\\ICIO\\'
path = 'D:\\Downloads\\Sufang\\WIOD\\'

safe_nodes = {}
for fname in os.listdir(path)[:5]:
    # if not fname.endswith('xlsb'):
    if fname.endswith('xlsb'):
        year = fname.split('_')[0][-4:]

        main_data = load_data(path + fname, kind='WIOD')
        # main_data = load_data('D:\\Downloads\\Sufang\\WIOD\\WIOT2002_Nov16_ROW.xlsb', kind='WIOD')
    elif fname.endswith('csv'):
        year = fname.split('_')[1].split('.')[0]
        main_data = load_data(path + fname, kind='OECD')
    else:
        continue
    the_net = nx.from_pandas_adjacency(main_data.astype(float).T, create_using=nx.DiGraph)
    the_net = nx.subgraph(the_net, max(nx.weakly_connected_components(the_net), key=len))

    most_of_the_net = 0.98 * len(the_net)

    for n, d in the_net.nodes(data=True):
        country, industry = n.split('_')
        d['country'] = country
        d['industry'] = industry

    for n, d in the_net.nodes(data=True):
        for t in the_net.successors(n):
            if 'money_out' not in d.keys():
                d['money_out'] = the_net.edges[n, t]['weight']
            else:
                d['money_out'] += the_net.edges[n, t]['weight']

        for t in the_net.predecessors(n):
            d['local_trade'] = 0
            d['inter_trade'] = 0
            if 'money_in' not in d.keys():
                d['money_in'] = the_net.edges[t, n]['weight']
            else:
                d['money_in'] += the_net.edges[t, n]['weight']
            if the_net.nodes[t]['country'] == d['country']:
                d['local_trade'] += the_net.edges[t, n]['weight']
            else:
                d['inter_trade'] += the_net.edges[t, n]['weight']

    the_net_copy = the_net.copy()

    init_nodes = [n for n, d in the_net.nodes(data=True) if d['country'] in ['USA', 'CHN', 'DEU', 'CN1', 'CN2']]
    # init_nodes = the_net_copy.nodes
    # init_nodes = np.random.choice(the_net.nodes, 10)
    safe_nodes[year] = {}
    for init_node in init_nodes:
        the_net = the_net_copy.copy()
        current_country = the_net.nodes[init_node]['country']
        print(year, init_node)
        if init_node not in net_remaining.keys():
            net_remaining[init_node] = {}
            num_iters[init_node] = {}
        if year not in net_remaining[init_node].keys():
            net_remaining[init_node][year] = {}
            num_iters[init_node][year] = {}
        # safe_nodes[year][init_node] = {}

        for th in np.arange(0.03, 0.3, 0.01):
            the_net = the_net_copy.copy()

            # safe_nodes[year][init_node][th] = []
            if th not in safe_nodes[year].keys():
                subg = []
                for k, d in the_net_copy.nodes(data=True):
                    any_frag = any([(w['weight'] / d['money_in']) > th for _, _, w in the_net_copy.in_edges(k, data=True)])
                    num_frag = sum([(the_net_copy.edges[(k, node2)]['weight'] / the_net_copy.nodes[node2]['money_in']) > th for node2 in
                                                      the_net_copy.successors(k)])
                    if num_frag > 2 and any_frag is True:
                        subg.append(k)
                sub_net = nx.subgraph(the_net_copy, subg)
                avdeg = np.median([d for _, d in sub_net.degree])
                # avdeg = sum([d for _, d in sub_net.degree]) / len(sub_net)
                safe_nodes[year][th] = [n for n, d in sub_net.degree() if d > avdeg]


            nodes_to_check = set(n for n in the_net.successors(init_node) if n not in safe_nodes[year][th]) # All nodes case
            # nodes_to_check = set(n for n in the_net.successors(init_node) if the_net.nodes[n]['country'] != current_country)
            the_net.remove_node(init_node)
            removed_nodes = []
            itr = 0
            while len(nodes_to_check) > 0:
                node_to_test = nodes_to_check.pop()
                if node_to_test in the_net:
                    itr += 1
                    # current_country = the_net.nodes[node_to_test]['country']
                    # current_value = sum([the_net[node_to_test][i]['weight'] for i in the_net.successors(node_to_test)])
                    current_value = sum([the_net[i][node_to_test]['weight'] for i in the_net.predecessors(node_to_test)])
                    if (current_value / the_net.nodes[node_to_test]['money_in']) < (1 - th):
                        # if 'money_in' in the_net.nodes[node_to_test].keys() and (current_value / the_net.nodes[node_to_test]['money_in']) < (1 - th):
                        # nodes_to_check |= set(the_net.successors(node_to_test))
                        nodes_to_check |= set(n for n in the_net.successors(node_to_test) if n not in safe_nodes[year][th])
                        the_net.remove_node(node_to_test)
                        removed_nodes.append(node_to_test)
            # print(init_node, th, len(the_net))
            net_remaining[init_node][year][th] = len(the_net)
            num_iters[init_node][year][th] = itr
            # if th >= max(net_remaining_no_protect[init_node][year].keys()):
            #     break
            if len(the_net) > most_of_the_net:
                break

pc_agg = {}
for top_node, values in net_remaining.items():
    pc_agg[top_node] = {}
    for year, threshes in values.items():
        if len(threshes) > 1:
            pc_agg[top_node][year] = pd.Series(threshes).diff().idxmax()
            # print(top_node, year, pd.Series(threshes).diff().max(), )
dfagg = pd.DataFrame(pc_agg)
dfagg.dropna(axis=1, how='all', inplace=True)
dfagg[[c for c in dfagg.columns if 'CHN' in c]].max(axis=1).plot(label='CHN', legend=True, style='.-')
dfagg[[c for c in dfagg.columns if 'DEU' in c]].max(axis=1).plot(label='DEU', legend=True, style='.-')
dfagg[[c for c in dfagg.columns if 'USA' in c]].max(axis=1).plot(label='USA', legend=True, style='.-')
