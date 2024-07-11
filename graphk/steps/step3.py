import numpy as np
import sys
from tqdm import tqdm
from functools import lru_cache
from collections import defaultdict


def get_bin(ref_vertex, classes):
    return classes[ref_vertex]

def get_connected_vertices(vertex, graph):
    return graph[vertex]

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
        #if (classes[ref_vertex] != vertices_info[ref_vertex]['bin']):
            #print("Mis match found. Label updating for node: ", ref_vertex)
        classes[ref_vertex] = vertices_info[ref_vertex]['bin']

    np.savez(new_classes_file, classes=classes)

def annotate_bins(exp_dir, out_dir):
    initial_tool_results = exp_dir

    vertices_file = out_dir + "/misbinned_reads.npy"
    marker_scores = out_dir + "/marker_scores.txt"
    edges_file = initial_tool_results + "/edges.npy"
    classes_file = out_dir + '/relabelled_clusters.npy'

    marker_scores_dict = {(int(line.split()[0].split('_')[1])-1): line.split()[1] for line in open(marker_scores)}
    edges_data = np.load(edges_file)
    ambigous_vertices = np.load(vertices_file)
    classes = np.load(classes_file)

    graph = defaultdict(list)
    for u, v in tqdm(edges_data, desc="Adding edges"):
        graph[u].append(v)
        graph[v].append(u)

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

            connected_vertices = get_connected_vertices(vertex, graph)
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

def run(exp_dir, out_dir):
    try:
        annotate_bins(exp_dir, out_dir)
    except Exception as e:
        print(f"Error step 3: {e}")

    