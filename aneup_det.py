"""
File Name: aneup_det.py
Description: The function of this file is to detect aneuploidy using pan-genomic graph.
Author: Xavier Xue
Date: Aug 13
Python Version: 3.7
"""

import random
import time
import os
from utils import *
import numpy as np
import json
import pandas as pd
from tqdm import trange
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-p', type=str, help="pan-genome file path", required=True)
parser.add_argument('-n', type=str, help="Path to the NIPT_sample.csv", required=True)
parser.add_argument('-l', type=float, default=-0.8, help="left threshold", required=True)
parser.add_argument('-r', type=float, default=1.5, help="right threshold", required=True)
args = parser.parse_args()

# Calculate the threshold for copy number on the left.
def left_threshold(aligned_file, chromo):

    aligned_chr_reads_ids = [element['reads_id'] for element in aligned_file if element['aligned_chr_name'] == chromo]
    filtered_elements = [element for element in aligned_file if element['reads_id'] in aligned_chr_reads_ids]

    reads_id_counts = {}
    for element in filtered_elements:
        reads_id = element['reads_id']
        if reads_id in reads_id_counts:
            reads_id_counts[reads_id] += 1
        else:
            reads_id_counts[reads_id] = 1

    reads_id_to_remove = set()

    for reads_id, count in reads_id_counts.items():
        if count > 1:
            aligned_chr_names = set()
            for element in filtered_elements:
                if element['reads_id'] == reads_id:
                    aligned_chr_names.add(element['aligned_chr_name'])
            if len(aligned_chr_names) > 1:
                reads_id_to_remove.add(reads_id)

    _filtered_elements = [element for element in filtered_elements if element['reads_id'] not in reads_id_to_remove]

    grouped_elements = {}
    for element in _filtered_elements:
        reads_id = element['reads_id']
        if reads_id in grouped_elements:
            grouped_elements[reads_id].append(element)
        else:
            grouped_elements[reads_id] = [element]

    filtered_elements_left = []
    for reads_id, elements in grouped_elements.items():
        sampled_element = random.choice(elements)
        filtered_elements_left.append(sampled_element)

    return filtered_elements_left

# Calculate the threshold for copy number on the right.
def right_threshold(aligned_file, chromo):

    aligned_chr_reads_ids = [element['reads_id'] for element in aligned_file if element['aligned_chr_name'] == chromo]
    filtered_elements = [element for element in aligned_file if element['reads_id'] in aligned_chr_reads_ids]
    filtered_elements_ = [element for element in filtered_elements if element['aligned_chr_name'] == chromo]

    grouped_elements = {}
    for element in filtered_elements_:
        reads_id = element['reads_id']
        if reads_id in grouped_elements:
            grouped_elements[reads_id].append(element)
        else:
            grouped_elements[reads_id] = [element]

    filtered_elements_right = []
    for reads_id, elements in grouped_elements.items():
        sampled_element = random.choice(elements)
        filtered_elements_right.append(sampled_element)

    return filtered_elements_right


def calculate_node_depth(filtered_elements, chr_len):
    node_depth = [0] * chr_len

    for item in filtered_elements:
        path = item['path']
        path_nodes = eval(path)
        for node in path_nodes:
            node_index = int(node[1:])
            if 1 <= node_index <= chr_len:
                node_depth[node_index - 1] += 1

    return node_depth


def z_score(pgg_depth, node_depth):

    mean_pgg = np.mean(pgg_depth)
    std_pgg = np.std(pgg_depth)
    mean_node = np.mean(node_depth)
    z = (mean_node - mean_pgg) / std_pgg

    return z


def tri_dec(mean_zscore_left, mean_zscore_right, threshold=(-1, 2)):
    if threshold[0] <= mean_zscore_left <= threshold[1] and threshold[0] <= mean_zscore_right <= threshold[1]:
        return 0
    elif mean_zscore_left < threshold[0]:
        return -1
    elif mean_zscore_right > threshold[1]:
        return 1


def main():
    aligned_path = 'data/aligned/'
    chr_len_path = 'data/chr_len.txt'
    chr_len_list = read_chr_len_file(chr_len_path)
    pgg_data_path = args.p
    threshold = (args.l, args.r)
    nipt_file_path = args.n
    df = pd.read_csv(nipt_file_path, index_col=0)
    with open(pgg_data_path, 'r') as file:
        pgg_data = json.load(file)
    detection_results_df = pd.DataFrame(columns=['sample'] + [item[0] for item in chr_len_list])

    for i in trange(df.shape[0]):
        time.sleep(0.01)
        aligned_file = read_aligned_sam_file(os.path.join(aligned_path, f"new_sample{i}.sam"))
        detection_results = []
        for idx, c in enumerate(chr_len_list):
            chr_len = c[1]
            filtered_elements_left = left_threshold(aligned_file, f'aligned_{c[0]}')
            filtered_elements_right = right_threshold(aligned_file, f'aligned_{c[0]}')
            node_depth_left = calculate_node_depth(filtered_elements_left, chr_len)
            node_depth_right = calculate_node_depth(filtered_elements_right, chr_len)

            pgg_depth = [int(pgg_data[idx][3]['v{}'.format(i)]) for i in range(1, c[1]+1)]
            pgg_depth = [value / df.shape[0] for value in pgg_depth]

            zscore_left = z_score(pgg_depth, node_depth_left)
            zscore_right = z_score(pgg_depth, node_depth_right)

            detection_result = tri_dec(zscore_left, zscore_right, threshold = threshold)
            detection_results.append(detection_result)

        detection_results_df.loc[i] = [f"new_sample{i}.sam"] + detection_results

    detection_results_df.to_csv('data/detection_results.csv', index=False)


if __name__ == '__main__':
    st = time.time()
    main()
    et = time.time()
    rt = et - st
    print(f"runtime: {rt}sec")
