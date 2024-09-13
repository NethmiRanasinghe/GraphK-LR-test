import pickle
from collections import defaultdict
from tqdm import tqdm
import numpy as np

def create_edges(initial_tool_results):

    read_id_idx = get_idx_maps(initial_tool_results + '/read_ids')

    # create and save the edge.npy file
    edges_nparr = alignments_to_edges(initial_tool_results + '/reads.alns', read_id_idx)

    np.save(initial_tool_results + '/edges.npy', edges_nparr)

    return edges_nparr


def get_idx_maps(read_ids_file_path):
    read_id_idx = {}
    
    with open(read_ids_file_path) as read_ids_file:
        for rid in tqdm(read_ids_file, desc="Read ids mapping"):
            rid = rid.split(" ")[0].strip()[1:]  #modified: rid = rid.strip()[1:] 
            read_id_idx[rid] = len(read_id_idx)
    
    return read_id_idx


def alignments_to_edges(alignments_file_path, read_id_idx):
    edges = []
    with open(alignments_file_path, "r") as af:
        for line in tqdm(af, desc="Alignments to edges"):
            u, v = line.strip().split('\t')
            if u != v:
                edges.append((read_id_idx[u], read_id_idx[v]))
    
    edges_nparr = np.array(edges, dtype=np.int32)
    return edges_nparr

def create_graph(edges):
    graph = defaultdict(list)
    for u, v in tqdm(edges, desc="Creating graph"):
        graph[u].append(v)
        graph[v].append(u)

    return graph

def run(exp_dir):
    edges = create_edges(exp_dir)
    
    #edges = np.load(exp_dir + 'edges.npy')

    # creating a dictionary with nodes as keys and neighbors as values
    graph = create_graph(edges)

    # Save the graph to a file
    with open(exp_dir + 'graph.pkl', 'wb') as f:
        pickle.dump(graph, f)



