# FireLab

# Grep Assembler for Short Reads - *GASR*

### For inquiries, contact dejeong@stanford.edu or sundrani@stanford.edu

### Files
This repo contains a `preprocess.py` and `{}.py` for swift assembly of short read sequence datasets given a reference seed sequence.

### Usage

#### Package Requirments
Our python script relies on the [BioPython](https://biopython.org) toolkit suite to quickly `grep` large short read sequence files. 
Before running any of the scripts, install the required dependencies from `requirements.txt` by running `pip install -r requirements.txt` or create a separate `conda` environment and run `conda install --file requirements.txt`. 

#### Running the script

##### Script Arguments
```
usage: ShortSeqAssembler_v7_0406_SS.py [-h] [--input INPUT] [--seedseq SEEDSEQ] [--nmer NMER] [--maxround MAXROUND]
                                       [--outdir OUTDIR] [--process_multiple_fastq PROCESS_MULTIPLE_FASTQ]

options:
  -h, --help            show this help message and exit
  --input INPUT         input fastq sequence (preprocessed or not) to assemble
  --seedseq SEEDSEQ     seed sequence from where we will begin our search (e.g. GTGCTTCACCAACGTGTACGCCGACAGCTTCGTG)
  --nmer NMER           size of kmer to complete searching rounds (default = 25)
  --maxround MAXROUND   maximum number of searching rounds (default = 65)
  --outdir OUTDIR       path to directory to save output files (default is cwd)
  --process_multiple_fastq PROCESS_MULTIPLE_FASTQ
                        True if input files need to be preprocessed. Must input multiple fastq files, unprocesed separated
                        by a space (default = False)
```

##### Examples:

###### With already trimmed and combined .fastq
```
python ShortSeqAssembler_v7_0406_SS.py --input 1_Mod500pgMS0pg_S1_L001_R1_001_final_combined.fastq --outdir example_output --seedseq GGTTCGACAACCCCGTGCTGCCCTTCAACGACGGCGTGTACTTC
```
###### With untrimmed pairs of .fastq
```
python ShortSeqAssembler_v7_0406_SS.py --input '1_Mod500pgMS0pg_S1_L001_R1_001.fastq 1_Mod500pgMS0pg_S1_L001_R2_001.fastq' --outdir example_output --process_multiple_fastq True --seedseq GGTTCGACAACCCCGTGCTGCCCTTCAACGACGGCGTGTACTTC
```

### Output
The script should take no more than 2-10 minutes to run. We have tested with large (~1 gb) input files and a variety of test input seed sequences. In the output folder specified (or `cwd` if not), we also output two figures for scrutiny. 

The leftmost figure displays a metric for confidence in sequence position up- and downstream of your input seed sequence (shown as the black bar). The right rightmost figure displays the number of logFold counts of reads we see in the short read .fastq file for each sequence position up- and downstream of the seed sequence. 

<img src="https://user-images.githubusercontent.com/48189633/162048169-f4bf7dc2-8bbf-40f9-80c2-67235f28665e.png" width=450 align=left>
<img src="https://user-images.githubusercontent.com/48189633/162048154-c8a38c74-3675-4cdd-a4fe-a47469853a88.png" width=450 align=right>




