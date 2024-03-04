# NIPT-PG
NIPT-PG: Empowering Non-Invasive Prenatal Testing to learn from population genomics through an incremental pan-genomic approach

## Step 1. Package Dependency
first, install the NIPT-PG conda environment:
```bash
conda create -c NIPT-PG
```
```bash
conda activate NIPT-PG
```
then, in NPIT-PG environment, install the following package:
```bash
pip install pandas numpy tqdm argparse
```

## Step 2. generating pan-genome
### usage:
```bash
python3 gen_pgg.py [-r REF.FA_FILE] [-s SAM_PATH] [-n NIPT_FILE]
```
### optional arguments: 
• -r path to the reference genome file (such as GRCh38.fa)

• -s path to the folder containing the files to be tested

• -n path to the nipt_files.csv

### example:
```bash
python3 gen_pgg.py -r data/ref.fa -s data/sam/ -n data/nipt_files_ART-Random.csv
```
The content of the nipt_files.csv file is as illustrated in **Table 1**, documenting the file name mappings for each testing file. This practice aids in standardizing file management and enhances testing efficiency.

#### Table 1. Illustration of the nipt_files.csv file.

| id | nipt_files | mapping |
|--------|--------|--------|
| 0   | CL100050702_L02_91   | sample_0   |
| 1   | CL100025607_L02_22   | sample_1   |
| 2   | CL100035831_L01_15   | sample_2   |
| ...   |...   | ...   |

## Step 3. Sequence-to-graph alignment
### usage:
```bash
python3 map2pgg.py [-p PGG_FILE] [-s SAM_PATH] [-n NIPT_FILE] [-k K_MER]
```
### optional arguments: 
• -p the pan-genome file path

• -s path to the folder containing the files to be tested

• -n path to the nipt_files.csv

• -k k-mer length, default=5

### example:
```bash
python3 map2pgg.py -p data/pgg.json -s data/sam/ -n data/nipt_files_ART-Random.csv -k 5
```

## Step 4. Z-score test based on multi-source aligned read
### usage:
```bash
python3 aneup_det.py [-s SAM_PATH] [-g ALIGNED_SAM_PATH] [-n NIPT_FILE] 
[-l LEFT_THRESHOLD] 
[-r RIGHT_THRESHOLD]
[-c CONTROL SAMPLE]
```
### optional arguments: 
• -s path to the folder containing the files to be tested

• -g path to the folder containing realigned samples

• -n path to the nipt_files.csv

• -l left threshold of z-score (default = -3)

• -r right threshold of z-score (default = 3)

• -c a file containing the mean and standard deviation of control samples (details in **Table 2**).

### example:
```bash
python3 trisomy_detection.py -s data/sam/ -g data/aligned_sam/ -n data/nipt_files_ART-Random.csv -l -3 -r 3 -c data/control_group (AR).csv
```
#### Table 2. The mean and standard deviation of control samples.

| chr | mean | sd |
|--------|--------|--------|
| chr13   | 0.040201   | 0.009946   |
| chr18   | 0.041326   | 0.009731   |
| chr21   | 0.041521   | 0.009685   |



