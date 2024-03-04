# NIPT-PG
NIPT-PG: Empowering Non-Invasive Prenatal Testing to learn from population genomics through an incremental pan-genomic approach

## Step 1. Package Dependency
first, install the NPIT-PG conda environment:
'''conda create -c NPIT-PG'''
'''conda activate NPIT-PG'''
then, in NPIT-PG environment, install the following package:
'''pip install pandas numpy tqdm argparse'''

## Step 2. generating pan-genome
usage:
python3 gen_pgg.py [-r REF.FA_FILE] [-s SAM_PATH] [-n NIPT_FILE]
optional arguments: 
-r path to the reference genome file (such as GRCh38.fa)
-s path to the folder containing the files to be tested
-n path to the nipt_files.csv
example:
python3 gen_pgg.py -r data/ref.fa -s data/sam/ -n data/nipt_files_ART-Random.csv
The content of the nipt_files.csv file is as illustrated in Table 1, documenting the file name mappings for each testing file. This practice aids in standardizing file management and enhances testing efficiency.

# Step 3. Sequence-to-graph alignment
usage:
python3 map2pgg.py [-p PGG_FILE] [-s SAM_PATH] [-n NIPT_FILE] [-k K_MER]
optional arguments: 
-p the pan-genome file path
-s path to the folder containing the files to be tested
-n path to the nipt_files.csv
-k k-mer length, default=5
example:
python3 map2pgg.py -p data/pgg.json -s data/sam/ -n data/nipt_files_ART-Random.csv -k 5

# Step 4. Z-score test based on multi-source aligned read
usage:
python3 aneup_det.py [-s SAM_PATH] [-g ALIGNED_SAM_PATH] [-n NIPT_FILE] 
[-l LEFT_THRESHOLD] 
[-r RIGHT_THRESHOLD]
[-c CONTROL SAMPLE]
optional arguments: 
-s path to the folder containing the files to be tested
-g path to the folder containing realigned samples
-n path to the nipt_files.csv
-l left threshold of z-score (default = -3)
-r right threshold of z-score (default = 3)
-c a file containing the mean and standard deviation of control samples (details in Table 2).
example:
python3 trisomy_detection.py -s data/sam/ -g data/aligned_sam/ -n data/nipt_files_ART-Random.csv -l -3 -r 3 -c data/control_group (AR).csv



