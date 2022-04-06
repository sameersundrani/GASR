"""
Preprocess 2 fastq files using Bio.SeqIO. For one fastq, takes about 60 seconds.

Example script call:
python preprocess.py 9_Mod32ngMS640pg_S9_L001_R1_001.fastq 9_Mod32ngMS640pg_S9_L001_R2_001.fastq
"""
import os
import sys
from Bio import SeqIO
from datetime import datetime

def trim_seqs(records, lenToRemove):
    """Removes first 3-10 bp from start of sequence.

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for record in records:
        yield record[lenToRemove:]

def main(input_files):

    saved_files = []
    for input in input_files:
        orig = SeqIO.parse(input, "fastq")
        trimmed_reads = trim_seqs(orig, 5)
        new_file_name = input[:-6] + '_' + "trimmed.fastq"
        count = SeqIO.write(trimmed_reads, new_file_name, "fastq")
        saved_files.append(new_file_name)
        print("Saved {} reads to {}".format(count, new_file_name))

    combine = 'cat ' + saved_files[0] + ' ' + saved_files[1] + ' > ' + input_files[0][:-6] + '_final_combined.fastq'
    os.system(combine)
    print("Saved trimmed and combined reads to {}".format(input_files[0][:-6] + '_final_combined.fastq'))

if __name__=='__main__':
    print('--Starting Sequence Trimming--')
    start = datetime.now()
    main([sys.argv[1], sys.argv[2]])
    print('--Ending Sequence Trimming--')
    end = datetime.now()
    print("--Total run time was {}--".format(str(end - start)))





