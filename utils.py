import time

def read_fasta_file(filename):
    sequences = {}
    current_sequence = ""

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    sequences[sequence_name] = current_sequence
                    current_sequence = ""
                sequence_name = line[1:]
            else:
                current_sequence += line

        if current_sequence:
            sequences[sequence_name] = current_sequence

    return sequences


def read_fasta_file2(filename):
    with open(filename, 'r') as file:
        sequence_name = ""
        current_sequence = ""
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    yield sequence_name, current_sequence
                    current_sequence = ""
                sequence_name = line[1:]
            else:
                current_sequence += line
        if current_sequence:
            yield sequence_name, current_sequence



# .sam
def read_sam_file(filename):
    alignments = []

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if not line.startswith('@'):  # Ignore lines starting with '@' (header information)
                fields = line.split('\t')
                alignment = {
                    'QNAME': fields[0],
                    'FLAG': int(fields[1]),
                    'RNAME': fields[2],
                    'POS': int(fields[3]),
                    'MAPQ': int(fields[4]),
                    'CIGAR': fields[5],
                    'RNEXT': fields[6],
                    'PNEXT': int(fields[7]),
                    'TLEN': int(fields[8]),
                    'SEQ': fields[9],
                    'QUAL': fields[10],
                }
                alignments.append(alignment)

    return alignments


def read_aligned_sam_file(filename):
    alignments = []

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if not line.startswith('@'):  # Ignore lines starting with '@' (header information)
                fields = line.split('\t')
                alignment = {
                    'reads_id': int(fields[0].replace('reads', '')),
                    'q_name': fields[1],
                    'flag': int(fields[2]),
                    'mapq': int(fields[3]),
                    'aligned_chr_name': fields[4],
                    'pos': int(fields[5]),
                    'prob': float(fields[6]),
                    'count': int(fields[7]),
                    'seq': fields[8],
                    'qual': fields[9],
                }
                alignments.append(alignment)
    return alignments


def read_chr_len_file(file_path):
    result = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line:
                parts = line.split()
                if len(parts) == 2:
                    chrom, length = parts
                    result.append([chrom, int(length)])

    return result



