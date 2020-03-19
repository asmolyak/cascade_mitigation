import numpy as np
import networkx as nx
from networkx.algorithms.community import girvan_newman, label_propagation_communities, greedy_modularity_communities
import community


def get_sizes(n, number_of_blocks=10, kind=None, a=1):
    interim = None
    if kind == 'linear':
        decay = a * np.arange(number_of_blocks)
        first_block = (n + sum(decay)) / number_of_blocks
        interim = [int(i) for i in np.array([first_block] * number_of_blocks) - decay]
    elif kind=='expon':
        if a > 1:
            print('can\'t have a coef>1 for exponential decay')
            return None
        decay = [a ** i for i in range(number_of_blocks)]
    elif kind == 'power':
        decay = np.array([1] + [i ** -a for i in range(2, number_of_blocks+1)])
    else:
        decay = np.array([1] * number_of_blocks)
    if interim is None:
        first_block = n/sum(decay)
        interim = [int(np.round(i)) for i in np.array([first_block] * number_of_blocks) * decay]
    if any(i<1 for i in interim):
        print('decay too steep')
        return None
    return np.array(interim)


def get_pees(number_of_blocks, number_of_diags=None, connectivity_decay=1.0):
    p = np.eye(number_of_blocks)
    if number_of_diags:
        for i in range(1, number_of_diags+1):
            p += np.diag([connectivity_decay**i]*(number_of_blocks-i), i)
            p += np.diag([connectivity_decay**i]*(number_of_blocks-i), -i)
    return p


def score_similarity(true_partition, community_partition):
    global_scores = {}
    for k, v in true_partition.items():
        scores = {}
        for node in v:
            if community_partition[node] in scores.keys():
                scores[community_partition[node]] += 1
            else:
                scores[community_partition[node]] = 1
        for num, val in enumerate(sorted(scores.values(),reverse=True)):
            if val > 1:
                if k in global_scores.keys():
                    global_scores[k] += num * val
                else:
                    global_scores[k] = num * val
    return global_scores


def get_block_degree(g, number_of_blocks):
    block_degree = {}
    for i in range(number_of_blocks):
        block_nodes = [n for n,d in g.nodes(data=True) if d['block']==i]
        subgraph = nx.subgraph(g, block_nodes)
        block_degree[i] = 2*len(subgraph.edges)/len(subgraph.nodes)
    return block_degree


def get_nodes_by_partition(g):
    true_part = {}
    for n, d in g.nodes(data=True):
        true_part[n] = d['block']
    flipped_tp = {}
    for k, v in true_part.items():
        if v in flipped_tp.keys():
            flipped_tp[v].append(k)
        else:
            flipped_tp[v] = [k]
    return flipped_tp


def get_network(network_size, average_degree, number_of_blocks, kind, number_of_diagonals, decay):
    p = get_pees(number_of_blocks=number_of_blocks, number_of_diags=number_of_diagonals, connectivity_decay=0.1)
    sizes = get_sizes(n=network_size, number_of_blocks=number_of_blocks, kind=kind, a=decay)
    resulting_degree = sizes@p@sizes.T/sum(sizes)
    factor = average_degree / resulting_degree
    p_deg = p*factor
    g = nx.stochastic_block_model(sizes=sizes,p=p_deg)
    return g, sizes


def test():
    network_size = 10000
    number_of_blocks = 20
    desired_degree = 8

    degrees = {}
    size_dict = {}
    score_dict = {}
    for a in np.arange(1, 0.8, -0.01):
         # np.array([4000,2500,2000,1500,500,50])
        g, sizes = get_network(network_size=network_size,average_degree=desired_degree,number_of_blocks=number_of_blocks,kind='expon',number_of_diagonals=2,decay=a)
        flipped_tp = get_nodes_by_partition(g)

        # gmc = greedy_modularity_communities(g)
        # lpa = label_propagation_communities(g)
        # gn = girvan_newman(g)
        part = community.best_partition(g)
        size_dict[a] = sizes
        degrees[a] = get_block_degree(g, number_of_blocks)
        score_dict[a] = score_similarity(flipped_tp, part, sizes)


if __name__=="__main__":
    test()