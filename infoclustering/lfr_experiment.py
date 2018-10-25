import argparse
import networkx as nx
from data_preprocess import get_graph_from_data
from infoclustering import InfoClustering
import os
import pickle
from joblib import Parallel, delayed
 
 
def run(core, experiment_path, label_path, weighted, verbose=1):
    def save_experiment(experiment_file):
        if experiment_file[-3:] == 'dat':
            graph_name = experiment_file[:-4]
            data_path = os.path.join(experiment_path, experiment_file)
        else:
            return 'not dat file'
        try:
            G = get_graph_from_data(data_path=data_path, weighted=weighted)
            print('Graph Loaded Success: {}'.format(graph_name))
        except NotADirectoryError:
            print('Not a directory: ', data_path)
            return 'not a directory'
        except FileNotFoundError:
            print('FileNotFoundError: ', data_path)
            return 'file not found'
        # run the algorithm
        IC = InfoClustering(n_clusters=2, verbose=False)
        IC.fit(G)
        with open(os.path.join(experiment_path, graph_name + '_infoclustering.solutions'), 'wb') as f:
            pickle.dump(IC.psp, f)
        return 'Good'

    x = Parallel(n_jobs=core, verbose = 2)(delayed(save_experiment)(experiment_file=exp_file) for exp_file in sorted(os.listdir(experiment_path), reverse=True))

    print(x)
     
 
 
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
    run(core=args.cpu, experiment_path=args.data, label_path=args.label, weighted=args.weighted)