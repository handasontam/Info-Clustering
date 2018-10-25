import argparse
import networkx as nx
from infoclustering.data_preprocess import get_graph_from_data
from infoclustering.infoclustering import InfoClustering
 
 
def run(core, data_path, label_path, weighted, verbose=1):
    """
    Process the data into networkx DiGraph and run the algorithm
    :param core:
    :param data_path:
    :param directed:
    :param beta:
    :return:
    """
    # process the data
    G = get_graph_from_data(data_path=data_path, weighted=weighted)
    # run the algorithm
    IC = InfoClustering(n_clusters=2, verbose=False)
    IC.fit(G)
    print(IC.psp)
    print(IC.gammas)
 
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Info-Clustering')
    parser.add_argument('--cpu', type=int, required=True, help='number of cpu core to run the algorithm in parallel')
    parser.add_argument('--data', type=str, required=True, help='the file path of the graph data')
    parser.add_argument('--label', type=str, required=False, help='the file path of the label')
    # boolean argument deciding whether the graph is weighted/unweighted
    parser.add_argument('--weighted', dest='weighted', action='store_true')
    parser.add_argument('--unweighted', dest='weighted', action='store_false')
    parser.set_defaults(weighted=False)
 
    args = parser.parse_args()
 
    print(args)
    print(args.data)
    run(core=args.cpu, data_path=args.data, label_path=args.label, weighted=args.weighted)
