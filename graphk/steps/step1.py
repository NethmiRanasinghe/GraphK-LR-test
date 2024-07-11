import numpy as np
from collections import defaultdict, deque
from tqdm import tqdm 
import subprocess

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

def get_misbinned(in_file, initial_tool_results, output_folder):
    print(output_folder)

    # Load the edges
    edges = create_edges(initial_tool_results)

    mis_binned = []

    result = subprocess.run(['wc', '-l', initial_tool_results+'/read_ids'], capture_output=True, text=True)
    total_reads = int(result.stdout.split()[0])
    print(total_reads)

    # load the bins
    # 1) semibin2
    # clusters = read_tsv_to_clusters(initial_tool_results + 'contig_bins.tsv', total_reads)

    # 2) oblr
    # clusters = np.load(in_file)['classes']

    # 3) lrb
    # clusters = np.loadtxt(initial_tool_results + 'bins.txt', dtype=int) 

    # 4) metabcc
    # clusters = bins_to_npy(initial_tool_results + 'final.txt')
    
    clusters = np.load(in_file)

    # creating a dictionary with nodes as keys and neighbors as values
    graph = defaultdict(list)
    for u, v in tqdm(edges, desc="Creating graph"):
        graph[u].append(v)
        graph[v].append(u)

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

def run(in_file, exp_dir, out_dir):

    get_misbinned(in_file, exp_dir, out_dir)
