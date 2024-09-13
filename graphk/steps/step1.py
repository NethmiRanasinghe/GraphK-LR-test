import numpy as np
from collections import defaultdict, deque
from tqdm import tqdm 
import pickle

def read_tsv_to_clusters(tsv_file, total_reads):
    
    # Initialize the clusters array with -1
    clusters = np.full(total_reads, -1, dtype=int)

    with open(tsv_file, 'r') as file:
        file.readline()
        for line in file:
            read_id, bin = line.strip().split('\t')
            read_index = int(read_id.split('_')[1]) - 1
            clusters[read_index] = int(bin)

    return clusters


def bfs_label_scores(graph, clusters, start_node):

    queue = deque([(start_node, 1)])
    visited = set()
    label_scores = defaultdict(float)
    
    while queue:
        node, depth = queue.popleft()
        visited.add(node)
        
        for neighbor in set(graph[node]):
            if neighbor not in visited:
                if clusters[neighbor] != -1:
                    score = (2 ** (-depth))
                    label_scores[clusters[neighbor]] += score
                else:
                    queue.append((neighbor, depth + 1))
                visited.add(neighbor)
    
    return label_scores

def get_misbinned(in_file, output_folder, graph):
 
    mis_binned = []

    # load the bins
    clusters = np.load(in_file)
 
    # Find the mislabeled nodes and compute label scores
    for node, node_label in tqdm(enumerate(clusters), total=len(clusters), desc="Checking nodes"):

        # possible labels from directly connected nodes
        possible_labels =  set(clusters[neighbor] for neighbor in graph[node])

        # len(possible_labels)==0 means isolated

        # obviously mis binned nodes
        if len(possible_labels)>0 and -1 not in possible_labels and node_label not in possible_labels :
            mis_binned.append(node)
            clusters[node] = -1
        
        elif node_label != -1:

            # having a direct neighbor with same label
            if node_label in possible_labels:

                label_scores = defaultdict(float)

                # Compute label scores considering only immediate neighbors
                for label in possible_labels:
                    for neighbor in set(graph[node]):
                        if clusters[neighbor] == label:
                            score = (2 ** (-1))
                            label_scores[label] += score

            # nodes with the same label might not be directly connected, but connected through a path with unlabelled nodes
            else: 

                # Compute label scores considering all neighbors depthwise
                label_scores = bfs_label_scores(graph, clusters, node)

                # not even the indirectly connected nodes are having the same label
                if node_label not in label_scores and len(label_scores) > 0:
                    mis_binned.append(node)
                    clusters[node] = -1
            
            # Get the key of the maximum value
            current_label_score  = label_scores[node_label]*1.5
            keys_greater_than_threshold = [key for key, value in label_scores.items() if value > current_label_score]
            if (len(keys_greater_than_threshold) > 0):
                clusters[node] = keys_greater_than_threshold[0] # relabel

    
    np.save(output_folder + '/relabelled_clusters.npy', clusters)
    np.save(output_folder + '/misbinned_reads.npy', list(mis_binned))

    return mis_binned

def run(exp_dir, in_file, out_dir):

    # Load the graph from the file
    with open(exp_dir + 'graph.pkl', 'rb') as f:
        graph = pickle.load(f)

    get_misbinned(in_file, out_dir, graph)
