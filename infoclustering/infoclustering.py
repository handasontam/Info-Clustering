import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from itertools import combinations, chain
import logging
from networkx.algorithms.flow import shortest_augmenting_path
from sklearn.base import BaseEstimator, ClusterMixin
#from .data_structure import DisjointSetForest
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import os
from ctypes import *
curr_path = os.path.dirname(os.path.abspath(os.path.expanduser(__file__)))
lib1 = cdll.LoadLibrary(os.path.join(curr_path, 'ibfs_python3/libboost_python35.so.1.67.0'))
from .ibfs_python3 import ibfs_ext

import warnings

warnings.filterwarnings("ignore")

def max_threshold_clustering(full_linkage_matrix, min_cluster_size, num_cluster, n_data):
    """
        n-samples
        full_linkage_matrix has shape (n, 4): 
            column 1 and column 2 represent cluster id of u and v that merge together,
            column 3 represent the threshold (negative min cut),
            column 4 represent the merged cluster size
        Reduce the threshold to the largest value at which we have 
        k disjoint clusters of size at least m.
    """
    current_num_cluster = 0
    last_threshold = 0
    n = full_linkage_matrix.shape[0]
    for i, split in enumerate(reversed(full_linkage_matrix)):
        threshold = split[2]
        if threshold == last_threshold:
            # no need to do the cluster again if it has been done in previous iterations.
            continue
        label = fcluster(full_linkage_matrix, threshold, criterion='distance')
#         print(label)
        _, cluster_size = np.unique(label, return_counts=True)
        if (np.max(label) == num_cluster) and (np.min(cluster_size) >= min_cluster_size):
            cluster_above_threshold = np.append(full_linkage_matrix[n-i:, 0].flatten(), 
                                                full_linkage_matrix[n-i:, 1].flatten())
            cluster_above_threshold = cluster_above_threshold[cluster_above_threshold < n_data].astype(np.int)
            label[cluster_above_threshold] = -1
            
            return label
    return np.full(n_data, -1)


def print_psp(setofsetofset):
    for _set in setofsetofset:  # each partition
        __set = [set(___set) for ___set in _set]
        print(__set)

def print_partition(setofset):
    print([set(__set) for __set in setofset])

def Fset(it):
    # when creating a (set of set), the inner sets cannot be mutable, 
    # therefore, we need to use frozenset to make it immutable
    return frozenset(it)


def g_gamma(G, gamma, partition, capacity):
    # example: partition = [{1}, {2}, {3}],  Graph = Fig 2 in the paper,  gamma=0.5
    # output: 7 - 1.5 = 5.5
    def g_gamma_helper(_set):
        if _set:
            vertices = list(_set)
            all_in_edges = G.in_edges(nbunch=vertices,
                                      data=True)  # e.g. [(1, 2, {'weight': 1}), (1, 3, {'weight': 5}), (2, 3, {'weight': 1})]\
            return np.sum([attribute_dict[capacity] for v1, v2, attribute_dict in all_in_edges if
                           (v1 not in vertices) or (v2 not in vertices)]) - gamma
        else:
            return 0  # g_gamma({}) is 0

    g = np.sum([g_gamma_helper(_set) for _set in partition])  # sum all incut for all set in the partition
    return g


def split(G, partition_1, partition_2, psp=set(), in_cuts=set(), gamma=set(), verbose=False):
    ''' recursively find all the split in the Dilworth truncation '''

    # step1 find intersection point
    g_p1 = g_gamma(G, 0, partition_1, 'weight')
    g_p2 = g_gamma(G, 0, partition_2, 'weight')
    cardinality_p1 = len(partition_1)
    cardinality_p2 = len(partition_2)
    gamma_bar = (g_p1 - g_p2) / (cardinality_p1 - cardinality_p2)
    g_gamma_bar = (cardinality_p2 * g_p1 - cardinality_p1 * g_p2) / (cardinality_p2 - cardinality_p1)

    # step2 compute g_gamma(gamma_bar)
    partition = find_partition(G, gamma_bar)
    cut_value = g_gamma(G, gamma_bar, partition, 'weight')

    # Print debug information
    if verbose:
        print('g_0_p1         ', g_p1)
        print('g_0_p2         ', g_p2)
        print('partition_1    ')
        print_partition(partition_1)
        print('partition_2')
        print_partition(partition_2)
        print('cardinality_p1: ', cardinality_p1)
        print('cardinality_p2: ', cardinality_p2)
        print('gamma_bar       ', gamma_bar)
        print('g_gamma_bar     ', g_gamma_bar)
        print('partition of intersection: ')
        print_partition(partition)
        print('cut value of intersection: ', cut_value)
        print('\n\n')
    # step3 determine the stopping criterion
    if (partition_1 == partition) or (partition_2 == partition) or (cut_value >= g_gamma_bar):
        psp = psp.union({partition_1})
        psp = psp.union({partition_2})
        in_cuts = in_cuts.union({g_p1, g_p2})
        gamma = gamma.union({gamma_bar})
    elif cut_value < g_gamma_bar:
        m = len(partition)
        left_psp, left_incuts, left_gamma = split(G, partition_1, partition, psp, in_cuts, gamma, verbose)
        right_psp, right_incuts, right_gamma = split(G, partition, partition_2, psp, in_cuts, gamma, verbose)
        psp = psp.union(left_psp)
        psp = psp.union(right_psp)
        in_cuts = in_cuts.union(left_incuts)
        in_cuts = in_cuts.union(right_incuts)
        gamma = gamma.union(left_gamma)
        gamma = gamma.union(right_gamma)
    return psp, in_cuts, gamma


def sfm_minimize(G, x, j):
    '''
    # input: A weighted digraph D on vertex set V with capacity function c: V^2 -> R
    # output: the cut_value, and the sink partition in the max-flow min-cut.
    '''
    j_aug = ibfs_ext.IBFSGraph()
    num_edges = 0
    for edge in G.subgraph(range(j)).edges():
        if (edge[0] != edge[1]):
            num_edges += 1
    j_aug.initSize(j, num_edges)
    for i in range(j):
        if G.has_edge(i, j):
            j_aug.addNode(i, max(0, -x[i]), max(0, x[i]) + G[i][j]['weight'])
        else:
            j_aug.addNode(i, max(0, -x[i]), max(0, x[i]))
    for edge in G.subgraph(range(j)).edges(data=True):
        if (edge[0] != edge[1]):
            j_aug.addEdge(edge[0], edge[1], edge[2]['weight'], 0)
    j_aug.initGraph()
    mf = j_aug.computeMaxFlow()
    sink_label = [j]
    for v in range(j):
        if not j_aug.isNodeOnSrcSide(v):
            sink_label.append(v)
    return mf, frozenset(sink_label)


    # temp_DG = nx.DiGraph()
    # temp_DG.add_node('src')
    # added_edge_to_sink = []
    # for v in range(0, j):
    #     # if x[v] is negative, add a edge from source to v
    #     temp_DG.add_edges_from([('src', v, {'c_gamma': max(0, -x[v])})])
    #     if G.has_edge(v, j):  # if x[v] is positive, add a edge from v to sink instead
    #         temp_DG.add_edges_from([(v, j, {'c_gamma': max(0, x[v]) + G[v][j]['weight']})])
    #     else:
    #         temp_DG.add_edges_from([(v, j, {'c_gamma': max(0, x[v])})])
    #     if x[v] > 0:
    #         added_edge_to_sink += [(v, x[v])]
    #     # add the edge from the original graph
    #     temp_DG.add_edges_from(
    #         [(v, w, {'c_gamma': G[v][w]['weight']}) for w in range(0, j) if not (w == v) and (G.has_edge(v, w))])
    # cut_value, partition = nx.minimum_cut(temp_DG, _s='src', _t=j, capacity='c_gamma')
    # src_partition = partition[0]
    # sink_partition = partition[1]
    # temp_DG.remove_node('src')
    # for node, added_wgt in added_edge_to_sink:
    #     temp_DG[node][j]['c_gamma'] = temp_DG[node][j]['c_gamma'] - added_wgt
    # cut_value_without_source = g_gamma(temp_DG, 0, [sink_partition], 'c_gamma')
    # return cut_value, sink_partition, cut_value_without_source

def find_partition(G, gamma):
    # initialize
    num_V = G.number_of_nodes()
    P = {Fset({0})}
    x = np.zeros(num_V)
    x[0] = g_gamma(G=G, gamma=gamma, partition=[{0}], capacity='weight')

    # find partition
    for j in range(1, num_V):
        # run the max-flow min cut
        max_flow, B_star = sfm_minimize(G, x, j)
        # update the P*
        to_add = B_star
        for setp in P:
            if (setp & B_star):
                # if the intersection is not empty
                to_add = setp.union(to_add)
                # remove every C that intersects B_star from P_star
                P = P - {setp}
        P.add(Fset(to_add))

        # update x
        x[j] = g_gamma(G, gamma, [B_star], 'weight') - np.sum([x[_j] for _j in B_star if _j != j])
    return Fset(P)


def update_PStar(self, BStar):
    # merge PStar and BStar
    to_add = BStar
    for setp in self.P:
        if (setp & BStar):
            # if the intersection is not empty
            to_add = setp.union(to_add)
            # remove every C that intersects B_star from P_star
            self.P = self.P - {setp}
    self.P.add(Fset(to_add))


def info_clustering(G, verbose=False):
    """Apply Info-Clustering on a adjacency matrix:

    Parameters
    -----------
    G : networkx Graph

    cluster_selection_method : string, optional (default='eom')
        Options are:
            * ``max_threshold`` : 
                Reduce the threshold to the largest value at which we have 
                k disjoint clusters of size at least m.
            * ``best_first`` : 
                Start at the root of a dendrogram, search for clusters of size 
                at least m in a best first search manner: choose the largest cluster 
                to split first, until we have at least k disjoint clusters or have already 
                reached the maximum possible number of clusters. It is desirable 
                to implement this efficiently as it can be part of the clustering 
                algorithm if the knowledge of k and m is available. There may be a 
                simpler infeasibility condition or a way to return the maximum possible 
                k given m, or maximum possible m given k.

    """
    if G.is_directed():
        import sys
        print('G must be undirected graph')
        sys.exit()
    G = nx.convert_node_labels_to_integers(G, first_label=0, label_attribute='old_label')
    G = G.to_directed()
    for u, v in G.copy().edges():
        if u > v:
            # only preserve edge from small index to bigger index
            G.remove_edge(u, v)
            G.add_edges_from([(u, v, {'weight': 0})])  # otherwise, ibfs will trigger segmentation fault
    N = G.number_of_nodes()

    # start with the trivial partition and the singleton partition
    psp, in_cuts, gammas = split(G,
              partition_1=Fset([Fset(list(range(N)))]),  # trivial partition
              partition_2=Fset([Fset([i]) for i in range(N)]),  # singleton partition
              psp=set(), 
              gamma=set(), 
              verbose=verbose)

    original_psp = []
    for partition in psp:
        original_partition = []
        for community in partition:
            original_community = []
            for vertex in community:
                original_vertex = G.nodes[vertex]['old_label']
                original_community = original_community + [original_vertex]
            original_partition = original_partition + [original_community]
        original_psp = original_psp + [original_partition]
    psp = original_psp
        
    # construct the dendrogram matrix based on the Principle Sequence Partition
    #DSF = DisjointSetForest(N+1)
    linkage_matrix = np.zeros((N - 1, 4))
    i = 0
    sorted_psp = sorted(psp, key=len, reverse=True)
    sorted_in_cuts = sorted(in_cuts, reverse=True)
    if verbose:
        print('sorted_psp: ', sorted_psp)
        print('sorted in-cuts: ', sorted_in_cuts)
    # for partition, in_cut in zip(sorted_psp, sorted_in_cuts):
    #     for cluster in partition:
    #         if len(cluster) > 1:
    #             for v1, v2 in zip(list(cluster)[:-1], list(cluster)[1:]):
    #                 v1_root = DSF.fast_find(v1)
    #                 v2_root = DSF.fast_find(v2)
    #                 if v1_root == v2_root:
    #                     # already merged in previous lower threshold
    #                     continue
    #                 else:
    #                     linkage_matrix[i][0] = v1_root
    #                     linkage_matrix[i][1] = v2_root
    #                     linkage_matrix[i][2] = -in_cut  # negative in-cut
    #                     linkage_matrix[i][3] = DSF.size[v1_root] + DSF.size[v2_root]
    #                     DSF.union(v1_root, v2_root)
    #                     i += 1
    # linkage_matrix[:, 2] = linkage_matrix[:, 2] - min(linkage_matrix[:, 2])
    # if verbose:
    #     print(linkage_matrix)
    #     dendrogram(linkage_matrix)
    #     print('linkage matrix constructed success')
    gammas = list(gammas)
    gammas.append(-np.inf)
    gammas.append(np.inf)
    gammas = sorted(gammas)
    solutions = {}
    for gamma_1, gamma_2, cluster in zip(gammas[:-1], gammas[1:], sorted(psp, key=lambda x: len(x))):
        solutions[(gamma_1, gamma_2)] = cluster
    
    return linkage_matrix, solutions

class InfoClustering(BaseEstimator, ClusterMixin):
    """
    Info-Clustering Description:
    
    Parameters
    -----------
    n_jobs : int, optional (default = 1)
        The number of parallel jobs to run.
        If ``-1``, then the number of jobs is set to the number of CPU cores.

    verbose: bool, optional (default = False)
        If ``True``, debug message will be printed

    Attributes
    ----------
    affinity_matrix_ : array-like, shape (n_samples, n_samples)
        Affinity matrix used for clustering. Available only if after calling
        ``fit``.

    labels_ :
        Labels of each point
    """

    def __init__(self, G, n_jobs=1, verbose=False):
        self.G = G
        self.n_jobs = n_jobs
        self.verbose = verbose

    def fit(self):
        self.linkage_matrix, self.solutions = info_clustering(G=self.G, 
                                                verbose=self.verbose)
        if self.verbose:
            dendrogram(self.linkage_matrix)
        return self

