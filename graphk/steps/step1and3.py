import numpy as np
from collections import defaultdict, deque
from tqdm import tqdm 

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

def create_graph(edges):
    graph = defaultdict(list)
    for u, v in tqdm(edges, desc="Creating graph"):
        graph[u].append(v)
        graph[v].append(u)

    return graph

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

# step 3

def get_bin(ref_vertex, classes):
    return classes[ref_vertex]

def find_marker_gene(read_id, marker_scores_dict):
    return marker_scores_dict.get(read_id, None)

def get_marker_genes_for_connected_vertices(connected_vertices, ambigous_vertices, marker_scores_dict):
    filtered_vertices = [vertex for vertex in connected_vertices if vertex not in ambigous_vertices]
    connected_vertices_info = {}
    marker_genes = [find_marker_gene(vertex, marker_scores_dict) for vertex in filtered_vertices]
    connected_vertices_info = {vertex: {'marker_gene': marker_gene} for vertex, marker_gene in zip(filtered_vertices, marker_genes)}
    return connected_vertices_info

def get_bins_for_connected_vertices(connected_vertices_info, ambigous_vertices, classes):
    filtered_vertices = {vertex: info for vertex, info in connected_vertices_info.items()
                         if vertex not in ambigous_vertices}
    
    bins = [get_bin(vertex, classes) for vertex in filtered_vertices.keys()]

    for vertex, info in filtered_vertices.items():
        info['bin'] = bins[list(filtered_vertices).index(vertex)]

    return connected_vertices_info

def update_bin(new_classes_file, vertices_info, classes):
    for ref_vertex in vertices_info:
        
        classes[ref_vertex] = vertices_info[ref_vertex]['bin']

    np.savez(new_classes_file, classes=classes)

def annotate_bins(out_dir, graph):

    vertices_file = out_dir + "/misbinned_reads.npy"
    marker_scores = out_dir + "/marker_scores.txt"
    classes_file = out_dir + '/relabelled_clusters.npy'

    marker_scores_dict = {(int(line.split()[0].split('_')[1])-1): line.split()[1] for line in open(marker_scores)}
    # edges_data = np.load(edges_file)
    ambigous_vertices = np.load(vertices_file)
    classes = np.load(classes_file)

    actually_updated_count = 0
    replaced_as_unlabelled = 0
    relabel_with_different_bin = 0
    mis_binned_without_markers = 0
    vertices_info = {}

    for ambigous_vertex in tqdm(ambigous_vertices, desc="Analyzing misbinned reads"):
        vertex = int(ambigous_vertex)
        
        old_bin = get_bin(vertex, classes)
        
        marker_gene = find_marker_gene(vertex, marker_scores_dict)
        

        if marker_gene:
            vertex_bin = -1

            connected_vertices = graph[vertex]
            connected_read_marker_info = get_marker_genes_for_connected_vertices(connected_vertices, ambigous_vertices, marker_scores_dict)
            connected_read_marker_info = get_bins_for_connected_vertices(connected_read_marker_info, ambigous_vertices, classes)

            marker_genes = [info['marker_gene'] for info in connected_read_marker_info.values()]
            bins = [info['bin'] for info in connected_read_marker_info.values()]

            if marker_gene in marker_genes:
                same_marker_gene_bins = [bin for bin, mg in zip(bins, marker_genes) if mg == marker_gene]
                if len(set(same_marker_gene_bins)) == 1:
                    vertex_bin = same_marker_gene_bins[0]

        else:
            vertex_bin = old_bin
            mis_binned_without_markers += 1

        if vertex_bin != -1 and vertex_bin != old_bin:
            relabel_with_different_bin += 1

        if vertex_bin != old_bin:
            actually_updated_count += 1

        if vertex_bin == -1:
            replaced_as_unlabelled += 1

        vertices_info[vertex] = {'bin': vertex_bin}
    
    updated_classes_file = out_dir + "/refined_classes.npz"
    update_bin(updated_classes_file, vertices_info, classes)


def run(in_file, exp_dir, out_dir):

    # Load the edges
    edges = create_edges(exp_dir)

    # creating a dictionary with nodes as keys and neighbors as values
    graph = create_graph(edges)

    get_misbinned(in_file, out_dir, graph)

    try:
        annotate_bins(out_dir, graph)
    except Exception as e:
        print(f"Error: {e}")