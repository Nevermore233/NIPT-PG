"""
File Name: gen_pgg.py
Description: The purpose of this file is to generate a pan-genome from multiple .sam files.
Author: Xavier Xue
Python Version: 3.7
"""
import gc
import json
import warnings
import re
import os
import pandas as pd
from tqdm import trange
import time
from datetime import datetime
from utils import *
import argparse


warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser()
parser.add_argument('-r', type=str, help="Path to the reference genome file.", required=True)
parser.add_argument('-s', type=str, help="Path to the sam file.", required=True)
parser.add_argument('-n', type=str, help="Path to the NIPT_sample.csv", required=True)
parser.add_argument('-a', action='store_true', help="Add new samples to Pan-genome")
args = parser.parse_args()


def generate_graph(sequence):
    V = ['v{}'.format(i + 1) for i in range(len(sequence))]
    E = [(V[i], V[i + 1]) for i in range(len(sequence) - 1)]
    sigma = {V[i]: sequence[i] for i in range(len(sequence))}
    G = V, E, sigma
    return G


def add_missing_edges(G, position, count):
    V, E, sigma = G
    new_edge = ('v{}'.format(position - 1), 'v{}'.format(position + count))
    E.append(new_edge)
    G = V, E, sigma
    return G


def add_insertion_edges(G, position, count, seq, insert_base_pos):
    V, E, sigma = G
    new_node = {}
    for i in range(count):
        new_node[i] = 'v{}'.format(len(V) + i + 1)
        V.append(new_node[i])
        sigma[new_node[i]] = str(seq[insert_base_pos + i])
    if len(new_node) != 1:
        E = E + [(new_node[i], new_node[i + 1]) for i in range(len(new_node) - 1)]
    new_edge1 = ('v{}'.format(position-1), new_node[0])
    E.append(new_edge1)
    new_edge2 = (new_node[len(new_node) - 1], 'v{}'.format(position))
    E.append(new_edge2)
    G = V, E, sigma
    return G


def gen_pgg(G, position, sigar_seg, seq):
    insert_base_pos = 0

    if 'H' in sigar_seg[0]:
        sigar_seg.pop(0)

    if 'S' in sigar_seg[0]:
        count = int(re.findall(r'\d+', sigar_seg[0])[0])
        seq = seq[count:]
        sigar_seg.pop(0)

    for item in sigar_seg:
        count = int(re.findall(r'\d+', item)[0])
        if 'M' in item:
            position += count
            insert_base_pos += count
        if 'D' in item:
            G = add_missing_edges(G, position, count)
            position += count
        if 'I' in item:
            G = add_insertion_edges(G, position, count, seq, insert_base_pos)
            insert_base_pos += count

    return G


def remove_duplicate_edges(G):
    V, E, sigma = G
    unique_edges = list(set(E))
    G = V, unique_edges, sigma
    return G


def main():

    ref_seq_path = args.r
    sam_path = args.s
    nipt_file_path = args.n

    log_filename = "log/gen_pgg_log.txt" 
    df = pd.read_csv(nipt_file_path, index_col=0)
    pgg_data_path = 'data/pgg.json'

    with open("data/chr_len.txt", "r") as chr_file:
        chr_names = [line.split()[0] for line in chr_file]
    a = 0
    if a == 1:
        with open(pgg_data_path, 'r') as file:
            G = json.load(file)
    else:
        sequences = read_fasta_file(ref_seq_path)
        G = []
        for i in range(len(chr_names)):
            chr_name = chr_names[i]
            g_chr = generate_graph(sequences[chr_name])
            G.append(g_chr)
            del g_chr
            gc.collect()

    for i in trange(df.shape[0]):
        time.sleep(0.01)
        mapping = f"sample_{i}"
        sample = df.loc[df['mapping'] == mapping]['nipt_files'].values[0]
        filename = os.path.join(sam_path, f'{sample}.sam')
        if os.path.isfile(filename):
            try:
                alignments = read_sam_file(filename)
            except Exception as e:
                current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
                error_message = f"{current_time}:Error occurred while reading {filename}: {e}"
                print(error_message)
                # Write the error message to the log file
                with open(log_filename, "a") as log_file:
                    log_file.write(error_message + "\n")
                continue

            for idx, alignment in enumerate(alignments):
                r_name = alignment['RNAME']
                if r_name in chr_names:
                    position = alignment['POS']
                    sigar = alignment['CIGAR']
                    seq = alignment['SEQ']
                    index = chr_names.index(r_name)
                    if sigar != '*' and position != '*':
                        sigar_seg = re.findall(r'\d+[SHMDI]', sigar)
                        if len(sigar_seg) != 1:
                            G[index] = gen_pgg(G[index], position, sigar_seg, seq)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue
                gc.collect()
            del alignments
            gc.collect()

    print('start to remove duplicate edges:')
    for j in range(len(chr_names)):
        G[j] = remove_duplicate_edges(G[j])

    with open(pgg_data_path, 'w') as file:
        json.dump(G, file)
    print('finish！')


if __name__ == '__main__':
    st = time.time()
    main()
    et = time.time()
    rt = et - st
    print(f"runtime: {rt}sec")
