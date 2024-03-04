"""
File Name: map2pgg.py
Description: The function of this file is to realign the pan-genome using .sam files, record multiple alignment positions, and generate new .sam files.
Author: Xavier Xue
Python Version: 3.7
"""
import json
import time
import os
from utils import *
import warnings
import pandas as pd
import argparse
from collections import defaultdict
from datetime import datetime
from tqdm import trange
warnings.filterwarnings('ignore')
import random

parser = argparse.ArgumentParser()
parser.add_argument('-p', type=str, help="pan-genome file path", required=True)
parser.add_argument('-s', type=str, help="Path to the sam file.", required=True)
parser.add_argument('-n', type=str, help="Path to the NIPT_sample.csv", required=True)
parser.add_argument('-k', type=int, default=5, help="k-mer length", required=True)
args = parser.parse_args()


def generate_fragment_mapping(graph, k):
    fragment_mapping = {}
    stack = []
    for node in graph[0]:
        stack.append((node, [node]))

    while stack:
        node, current_fragment = stack.pop()
        if len(current_fragment) == k:
            fragment_mapping_key = ''.join([graph[2][node] for node in current_fragment])
            if fragment_mapping_key not in fragment_mapping:
                fragment_mapping[fragment_mapping_key] = [current_fragment.copy()]
            else:
                fragment_mapping[fragment_mapping_key].append(current_fragment.copy())
        else:
            neighbors = [edge[1] for edge in graph[1] if edge[0] == node]
            for neighbor in neighbors:
                stack.append((neighbor, current_fragment + [neighbor]))

    return fragment_mapping


def get_index_and_adjlist(pgg_data, k):
    indices = []
    adjacency_list = []
    for idx, pgg in enumerate(pgg_data):
        print(idx)
        index = generate_fragment_mapping(pgg, k)
        indices.append(index)

        adj = {}
        E = pgg[1]
        for edge in E:
            node, neighbor = edge
            if node in adj:
                adj[node].append(neighbor)
            else:
                adj[node] = [neighbor]
        adjacency_list.append(adj)

    return indices, adjacency_list


def calculate_path_probabilities(align_paths, chr_len_list):
    paths = defaultdict(lambda: {'prob': [0, 0], 'count': 0, 'w': 0, 'max_path': None})
    total_paths = len(align_paths)

    for align_path in align_paths:
        pos, path = align_path[0]
        chr_len = chr_len_list[pos][1]
        start_node = int(path[0][1:])
        w = sum(1 for node in path if node.startswith('v') and 1 <= int(node[1:]) <= chr_len)
        paths[(pos, start_node)]['count'] += 1
        if w > paths[(pos, start_node)]['w']:
            paths[(pos, start_node)]['w'] = w
            paths[(pos, start_node)]['max_path'] = path

    results = [
        [pos, start_node, path_info['count'] / total_paths, path_info['count'], path_info['max_path']]
        for (pos, start_node), path_info in paths.items()
    ]

    return results


def dfs(match_path, adj, sigma_upper, start_node, remaining_sequence, idx):
    stack = [(start_node, remaining_sequence, [start_node])]

    while stack:
        current_node, remaining_sequence, current_path = stack.pop()

        if not remaining_sequence:
            match_path.append((idx, current_path))
            continue

        next_char = remaining_sequence[0]
        remaining_sequence = remaining_sequence[1:]

        for neighbor in adj.get(current_node, []):
            if sigma_upper[neighbor] == next_char:
                stack.append((neighbor, remaining_sequence, current_path + [neighbor]))


def find_matching_paths(pgg_data, target_sequence, indices, adjacency_list, chr_len_list, k):
    align_paths_list = []
    for idx, pgg in enumerate(pgg_data):
        V, E, sigma = pgg
        target_sequence_upper = target_sequence.upper()
        sigma_upper = {v: sigma[v].upper() for v in V}
        index = indices[idx]
        adj = adjacency_list[idx]

        start_nodes = []
        match_paths = []
        for path in index.get(target_sequence_upper[:k], []):
            start_nodes.append(path[0])
        for start_node in start_nodes:
            match_path = []
            dfs(match_path, adj, sigma_upper, start_node, target_sequence_upper[1:], idx)
            if match_path:
                match_paths.append(match_path)
        if match_paths:
            align_paths_list.append(match_paths)

    align_paths = [item for sublist in align_paths_list for item in sublist]

    if align_paths:
        results = calculate_path_probabilities(align_paths, chr_len_list)
    else:
        results = []

    return results


def main():

    aligned_path = 'data/aligned/'
    if not os.path.exists(aligned_path):
        os.makedirs(aligned_path)
    pgg_data_path = args.p
    sam_path = args.s
    nipt_file_path = args.n
    k = args.k

    chr_len_path = 'data/chr_len.txt'
    chr_len_list = read_chr_len_file(chr_len_path)
    with open(chr_len_path, "r") as chr_file:
        chr_names = [line.split()[0] for line in chr_file]
    log_filename = "log/map2pgg_log.txt"  

    with open(pgg_data_path, 'r') as file:
        pgg_data = json.load(file)

    df = pd.read_csv(nipt_file_path, index_col=0)

    if not os.path.exists('data/indices.json') or not os.path.exists('data/adjacency_list.json'):
        indices, adjacency_list = get_index_and_adjlist(pgg_data, k)
        with open('data/indices.json', 'w') as indices_file:
            json.dump(indices, indices_file)
        with open('data/adjacency_list.json', 'w') as adjacency_list_file:
            json.dump(adjacency_list, adjacency_list_file)
    else:
        print("Both index.json and adjacency_list.json already exist.")
        with open('data/indices.json', 'r') as indices_file:
            indices = json.load(indices_file)
        with open('data/adjacency_list.json', 'r') as adjacency_list_file:
            adjacency_list = json.load(adjacency_list_file)

    for i in trange(df.shape[0]):
        mapping = f"sample_{i}"
        sample = df.loc[df['mapping'] == mapping]['nipt_files'].values[0]
        filename = os.path.join(sam_path, f'{sample}.sam')
        alignments = read_sam_file(filename)
        new_filename = os.path.join(aligned_path, f"new_sample_{i}.sam")
        with open(new_filename, 'w') as new_file:
            for idx, alignment in enumerate(alignments):
                try:
                    reads_id = f'reads{idx}'
                    rname = alignment['RNAME']
                    q_name = alignment['QNAME']
                    flag = alignment['FLAG']
                    mapq = alignment['MAPQ']
                    seq = alignment['SEQ']
                    qual = alignment['QUAL']

                    results = find_matching_paths(pgg_data, seq, indices, adjacency_list, chr_len_list, k)
                    print(results)
                    if results:
                        for result in results:
                            pos = result[0]
                            start_node = result[1]
                            prob = result[2]
                            count = result[3]
                            aligned_chr_name = f'aligned_{chr_len_list[pos][0]}'
                            new_file.write(f"{reads_id}\t{q_name}\t{flag}\t{mapq}\t{aligned_chr_name}\t{start_node}\t{prob}\t{count}\t{seq}\t{qual}\n")
                    elif rname in chr_names:
                        aligned_chr_name = f'aligned_{rname}'
                        new_file.write(
                            f"{reads_id}\t{q_name}\t{flag}\t{mapq}\t{aligned_chr_name}\t{start_node}\t{prob}\t{count}\t{seq}\t{qual}\n")
                    else:
                        continue
                except Exception as e:
                    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
                    error_message = f"{current_time}:Error occurred while reading {filename}, idx = {idx}, error: {e}"
                    print(error_message)
                    with open(log_filename, "a") as log_file:
                        log_file.write(error_message + "\n")
                    continue


if __name__ == '__main__':
    st = time.time()
    main()
    et = time.time()
    rt = et - st
    print(f"runtime: {rt}sec")