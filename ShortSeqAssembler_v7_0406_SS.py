"""
Fast Assembler (with) Sequences (of) Vaccine RNAseq
FASVrna
(EDIT)
"""
import sys
import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from Bio.Seq import Seq
from datetime import datetime
import argparse



def preprocess(input, trim=False):
    if trim:
        try:
            cmd = 'python preprocess.py ' + input
            os.system(cmd)
        except ValueError:
            print("Error with input format")

def right_extender(rightseqs):
    Position_list = [[0, 0, 0, 0, 0] for k in range(len(max(rightseqs, key=len)))]
    for rightseq in rightseqs:
        for i in range(len(rightseq)):
            Acount = 0
            Tcount = 0
            Ccount = 0
            Gcount = 0
            Ncount = 0

            if rightseq[i] == "A":
                Acount += 1
            elif rightseq[i] == "T":
                Tcount += 1
            elif rightseq[i] == "C":
                Ccount += 1
            elif rightseq[i] == "G":
                Gcount += 1
            else:
                Ncount += 1
            Position_list[i] = [x + y for x, y in zip(Position_list[i], [Acount, Tcount, Ccount, Gcount, Ncount])]
    return Position_list

def left_extender(leftseqs):
    Left_Position_list = [[0, 0, 0, 0, 0] for k in range(len(max(leftseqs, key=len)))]
    for leftseq in leftseqs:
        for i in range(len(leftseq) + 1):
            if i != 0:
                Acount = 0
                Tcount = 0
                Ccount = 0
                Gcount = 0
                Ncount = 0

                if leftseq[-i] == "A":
                    Acount += 1
                elif leftseq[-i] == "T":
                    Tcount += 1
                elif leftseq[-i] == "C":
                    Ccount += 1
                elif leftseq[-i] == "G":
                    Gcount += 1
                else:
                    Ncount += 1
                Left_Position_list[-i] = [x + y for x, y in zip(Left_Position_list[-i],
                                                                [Acount, Tcount, Ccount, Gcount, Ncount])]
            else:
                continue
    return Left_Position_list

def extension(file, nmer, searchseq, direction = 'R'):
    nucleotide = ["A", "T", "C", "G", "N"]
    with open(file) as f:
        seqlines = [seqline.strip() for seqline in f]
        seqs = []
        for seq in seqlines:
            index = seq.find(searchseq)
            if direction == 'R':
                rightseq = seq[index + nmer:]
                seqs.append(rightseq)
            else: #direction is 'L'
                leftseq = seq[:index]
                if index != 0:
                    seqs.append(leftseq)
                else:
                    continue

        Gainedseq = ""
        Probability_position = []

        if direction == 'R':
            Position_list = right_extender(seqs)
            for l in Position_list:
                probability = [(l[m] / sum(l)) * 100 for m in range(5)]
                Probability_position.append(probability)
                cand = max(probability)
                ind = probability.index(cand)
                if cand > 50:
                    Gainedseq = Gainedseq + nucleotide[ind]
                else:
                    Gainedseq = Gainedseq + "X"
            new_searchseq = Gainedseq[-(nmer + 5):-5]

            if "X" in Gainedseq[:-5]:
                print('-- reached a stopping point for extension --')
                raise ValueError

            return new_searchseq, Gainedseq[:-5], Probability_position[:-5], Position_list[:-5]

        else:  #direction == 'L'
            Left_Position_list = left_extender(seqs)
            for l in reversed(Left_Position_list):
                probability = [(l[m] / sum(l)) * 100 for m in range(5)]
                Probability_position.insert(0, probability)
                cand = max(probability)
                ind = probability.index(cand)
                if cand > 50:
                    Gainedseq = nucleotide[ind] + Gainedseq
                else:
                    Gainedseq = "X" + Gainedseq
            new_searchseq = Gainedseq[5:(nmer + 5)]

            if "X" in Gainedseq[5:]:
                print('-- reached a stopping point for extension --')
                raise ValueError

            return new_searchseq, Gainedseq[5:], Probability_position[5:], Left_Position_list[5:]

def right_extension(file, nmer, searchseq):

    with open(file) as f:
        seqlines = [ seqline.strip() for seqline in f ]
        rightseqs = []

        for seq in seqlines:
            index = seq.find(searchseq)
            rightseq = seq[index+nmer:]
            rightseqs.append(rightseq)
    
        Position_list = [ [0,0,0,0,0] for k in range(len(max(rightseqs, key=len)))]
    
        for rightseq in rightseqs:
            for i in range(len(rightseq)):
                Acount = 0
                Tcount = 0
                Ccount = 0
                Gcount = 0
                Ncount = 0
            
                if rightseq[i] == "A":
                    Acount += 1
                elif rightseq[i] == "T":
                    Tcount += 1
                elif rightseq[i] == "C":
                    Ccount += 1
                elif rightseq[i] == "G":
                    Gcount += 1
                else:
                    Ncount += 1
                Position_list[i] = [x+y for x,y in zip(Position_list[i], [Acount, Tcount, Ccount, Gcount, Ncount])]

        Percent_position = [[0,0,0,0,0] for k in range(len(max(rightseqs, key=len)))]
    
        Gainedseq = ""
        Probability_position = []
        for l in Position_list:
            probability = [(l[m] / sum(l))*100 for m in range(5)]
            Probability_position.append(probability)
            cand = max(probability)
            ind = probability.index(cand)
            if cand > 50:
                Gainedseq = Gainedseq + nucleotide[ind]
            else:
                Gainedseq = Gainedseq + "X"
        new_searchseq = Gainedseq[-(nmer+5):-5]

        if "X" in Gainedseq[:-5]:
            print('-- reached a stopping point for extension --')
            raise ValueError
        
    return new_searchseq, Gainedseq[:-5], Probability_position[:-5], Position_list[:-5]

def left_extension(file, nmer, Lsearchseq):
    
    with open(file) as f:
        seqlines = [ seqline.strip() for seqline in f ]
        leftseqs = []
    
        for seq in seqlines:
            index = seq.find(Lsearchseq)
            leftseq = seq[:index]
            if index != 0:
                leftseqs.append(leftseq)
            else:
                continue
        
        Left_Position_list = [ [0,0,0,0,0] for k in range(len(max(leftseqs, key=len)))]
    
        for leftseq in leftseqs:
            for i in range(len(leftseq)+1):
                if i != 0:
                    Acount = 0
                    Tcount = 0
                    Ccount = 0
                    Gcount = 0
                    Ncount = 0
                    
                    if leftseq[-i] == "A":
                        Acount += 1
                    elif leftseq[-i] == "T":
                        Tcount += 1
                    elif leftseq[-i] == "C":
                        Ccount += 1
                    elif leftseq[-i] == "G":
                        Gcount += 1
                    else:
                        Ncount += 1
                    Left_Position_list[-i] = [x+y for x,y in zip(Left_Position_list[-i], 
                                                                 [Acount, Tcount, Ccount, Gcount, Ncount])]
                else:
                    continue
        
        Left_Percent_position = [[0,0,0,0,0] for k in range(len(max(leftseqs, key=len)))]
    
        Gainedseq = ""
        Left_Probability_position = []
        for l in reversed(Left_Position_list):
            probability = [(l[m] / sum(l))*100 for m in range(5)]
            Left_Probability_position.insert(0, probability)
            cand = max(probability)
            ind = probability.index(cand)
        
            if cand > 50:
                Gainedseq = nucleotide[ind] + Gainedseq
            else:
                Gainedseq = "X" + Gainedseq
        new_searchseq = Gainedseq[5:(nmer+5)]
        
        if "X" in Gainedseq[5:]:
            print('-- reached a stopping point for extension --')
            raise ValueError
    
    return new_searchseq, Gainedseq[5:], Left_Probability_position[5:], Left_Position_list[5:]

def right_search(inputfile, seedseq, nmer, maxround):
    Rsearchseq = seedseq[-nmer:]
    Reads_position = pd.DataFrame(columns = ["A_counts", "T_counts", "C_counts", "G_counts", "N_counts"])
    nt_probability = pd.DataFrame(columns = ["A_percent", "T_percent", "C_percent", "G_percent", "N_percent"])

    extendedseq = ""
    search = True
    a = 1

    print("Extension from 3 prime of the seed sequence has been started.")

    while search == True:

        grepcmd = 'seqkit grep -s -P -p ' + Rsearchseq + ' ' + inputfile + ' | seqkit seq -s -o TEMP_seqlist.txt'
        os.system(grepcmd)
        revRsearchseq = str(Seq(Rsearchseq).reverse_complement())
        revgrepcmd = 'seqkit grep -s -P -p ' + revRsearchseq + ' ' + inputfile + ' | seqkit seq -s -o TEMP_rev_seqlist.txt'
        os.system(revgrepcmd)
        revrevfile = open("TEMP_revrev_seqlist.txt", "w")

        with open("TEMP_rev_seqlist.txt") as revlist:
            seqlines = [ seqline.strip() for seqline in revlist ]
            for seq in seqlines:
                seq = Seq(seq)
                revseq = seq.reverse_complement()
                revrevfile.write(str(revseq)+"\n")

        revrevfile.close()

        catcmd = 'cat TEMP_seqlist.txt TEMP_revrev_seqlist.txt > TEMP_total_seqlist.txt'
        os.system(catcmd)
        file = "TEMP_total_seqlist.txt"

        if a != maxround + 1:

            try:
                new, extended, nt_prob, Readdf = extension(file, nmer, Rsearchseq)
                if len(new) == nmer:
                    sconcat = datetime.now()
                    nt_probability = pd.concat([nt_probability, pd.DataFrame(nt_prob, columns=nt_probability.columns)], ignore_index=True)
                    Reads_position = pd.concat([Reads_position, pd.DataFrame(Readdf, columns=Reads_position.columns)], ignore_index=True)
                    extendedseq = extendedseq + extended
                    Rsearchseq = new
                    print("Round "+str(a)+" iteration for the extension to the right was completed.")
                    a += 1

                else:
                    print("Iteration for the right extension is over.")
                    break

            except ValueError:
                print("Iteration for the right extension is over.")
                break

        else:
            print("Iteration for the right extension has reached to the maximal number of rounds.")
            search = False

    print("Assembled sequence from 3 prime: ", extendedseq)
    print("Assembled sequence length: "+str(len(extendedseq))+"bp")

    return nt_probability, Reads_position, extendedseq

def left_search(inputfile, seedseq, nmer, maxround):
    Lsearchseq = seedseq[:nmer]
    LReads_position = pd.DataFrame(columns = ["A_counts", "T_counts", "C_counts", "G_counts", "N_counts"])
    Lnt_probability = pd.DataFrame(columns = ["A_percent", "T_percent", "C_percent", "G_percent", "N_percent"])

    Lextendedseq = ""
    search = True
    a = 1

    print("Extenstion from 5 prime of the seed sequence has been started.")

    while search == True:

        grepcmd = 'seqkit grep -s -P -p ' + Lsearchseq + ' ' + inputfile +' | seqkit seq -s -o TEMP_seqlist.txt'
        os.system(grepcmd)
        revLsearchseq = str(Seq(Lsearchseq).reverse_complement())
        revgrepcmd = 'seqkit grep -s -P -p ' + revLsearchseq + ' ' + inputfile + ' | seqkit seq -s -o TEMP_rev_seqlist.txt'
        os.system(revgrepcmd)
        revrevfile = open("TEMP_revrev_seqlist.txt", "w")

        with open("TEMP_rev_seqlist.txt") as revlist:
            seqlines = [ seqline.strip() for seqline in revlist ]
            for seq in seqlines:
                seq = Seq(seq)
                revseq = seq.reverse_complement()
                revrevfile.write(str(revseq)+"\n")

        revrevfile.close()

        catcmd = 'cat TEMP_seqlist.txt TEMP_revrev_seqlist.txt > TEMP_total_seqlist.txt'
        os.system(catcmd)
        file = "TEMP_total_seqlist.txt"


        if a != maxround + 1:

            try:
                Lnew, Lextended, Lnt_prob, LReaddf = extension(file, nmer, Lsearchseq, 'L')

                if len(Lnew) == nmer:
                    Lnt_probability = pd.concat([pd.DataFrame(Lnt_prob,  columns=Lnt_probability.columns), Lnt_probability], ignore_index=True)
                    LReads_position = pd.concat([pd.DataFrame(LReaddf, columns=LReads_position.columns), LReads_position], ignore_index=True)
                    Lextendedseq = Lextended + Lextendedseq
                    Lsearchseq = Lnew
                    print("Round "+str(a)+" iteration for the extension to the left was completed.")
                    a += 1

                else:
                    print("Iteration for the left extension is over.")
                    break

            except ValueError:
                print("Iteration for the left extension is over.")
                break

        else:
            print("Iteration for the left extension has reached to the maximal number of rounds.")
            search = False

    print("Assembled sequence at left: ", Lextendedseq)
    print("Assembled sequence length: "+str(len(Lextendedseq))+"bp")

    return Lnt_probability, LReads_position, Lextendedseq

def get_nucleotides(input_df, nucleotide):
    NT = []
    for index, rows in input_df.iterrows():
        count_list = [rows.A_counts, rows.T_counts, rows.C_counts, rows.G_counts, rows.N_counts]
        maxcount = max(count_list)
        ntindex = count_list.index(maxcount)
        NT.append(nucleotide[ntindex])
    return NT

def plot_figs(final_df, outdir):
    seeddf = final_df[final_df["A_counts"].isnull()]
    seeddf = seeddf.reset_index()
    seedstart = seeddf["Position"].iloc[0]
    seedend = seeddf["Position"].iloc[-1]

    ax = sns.lineplot(data=final_df, x="Position", y="log10TotalCounts")
    plt.plot([seedstart, seedend], [0.0, 0.0], color='black', linestyle='-', linewidth=2)

    plt.savefig(outdir + "/Coverage_basebybase_log10TotalCounts.png", format="png", dpi=1500, transparent=False, facecolor="white",
               bbox_inches="tight")

    plt.clf()

    ax2 = sns.lineplot(data=final_df, x="Position", y="MaxPercent")
    plt.plot([seedstart, seedend], [0.0, 0.0], color='black', linestyle='-', linewidth=2)

    plt.savefig(outdir + "/MaxPercentage_basebybase.png", format="png", dpi=1500, transparent=False, facecolor="white",
               bbox_inches="tight")

def main(inputfile, seedseq, nmer, maxround, outdir, process_multiple_fastq):

    Starttime=datetime.now()
    inputfile = str(inputfile)
    if process_multiple_fastq:
        preprocess(inputfile, trim=True)
        inputfile = inputfile.split()[0][:-6] + '_final_combined.fastq'

    print('Running with seed sequence= ' + seedseq + ' and with nmer = ' + str(nmer) +
          ' initiated at ' + str(Starttime))
    nucleotide = ["A", "T", "C", "G", "N"]
    num_to_print = 60 # For writing Fasta
    nt_probability, Reads_position, extendedseq = right_search(inputfile, seedseq, nmer, maxround)
    Lnt_probability, LReads_position, Lextendedseq = left_search(inputfile, seedseq, nmer, maxround)


    if not os.path.exists(outdir):
        # Create a new directory because it does not exist
        os.makedirs(outdir)
        print("Created a new output directory at {}!".format(outdir))

    print("Saving output files")

    nt_probability['MaxPercent'] = nt_probability[["A_percent", "T_percent", "C_percent", "G_percent", "N_percent"]].max(axis=1)
    Reads_position['MaxCounts'] = Reads_position[["A_counts", "T_counts", "C_counts", "G_counts", "N_counts"]].max(axis=1)
    NT = get_nucleotides(Reads_position, nucleotide)
    Reads_position["RepresentativeBase"] = NT
    concat = pd.concat([Reads_position,nt_probability], axis=1)
    concat.to_csv(outdir + "/Right_assembled_info.csv", sep='\t', index_label="Position")

    Lnt_probability['MaxPercent'] = Lnt_probability[["A_percent", "T_percent", "C_percent", "G_percent", "N_percent"]].max(axis=1)
    LReads_position['MaxCounts'] = LReads_position[["A_counts", "T_counts", "C_counts", "G_counts", "N_counts"]].max(axis=1)
    LNT = get_nucleotides(LReads_position, nucleotide)
    LReads_position["RepresentativeBase"] = LNT
    Lconcat = pd.concat([LReads_position,Lnt_probability], axis=1)
    Lconcat.to_csv(outdir + "/Left_assembled_info.csv", sep='\t', index_label="Position")

    final_assembled = Lextendedseq + seedseq + extendedseq

    print("Assembled sequence at both ends: ", final_assembled)
    print("Assembled sequence length: " + str(len(final_assembled)))

    final_output_file = open(outdir + "/Final_assembled.fasta", "w")
    final_output_file.write('>Final_assembled_'+ str(len(final_assembled))+"bp"+'\n')

    finalsplitseq = [final_assembled[i:i+num_to_print] for i in range(0, len(final_assembled), num_to_print)]

    for line in finalsplitseq:
        final_output_file.write(line+"\n")

    final_output_file.close()

    seed_df = pd.DataFrame(columns = ["RepresentativeBase", "A_counts", "T_counts", "C_counts", "G_counts",
                                      "N_counts", "MaxCounts", "A_percent", "T_percent", "C_percent", "G_percent",
                                      "N_percent", "MaxPercent"])
    seed_list = list(seedseq)
    seed_df["RepresentativeBase"] = seed_list

    concat = concat[["RepresentativeBase", "A_counts", "T_counts", "C_counts", "G_counts", "N_counts", "MaxCounts",
                     "A_percent", "T_percent", "C_percent", "G_percent", "N_percent", "MaxPercent"]]
    Lconcat = Lconcat[["RepresentativeBase", "A_counts", "T_counts", "C_counts", "G_counts", "N_counts", "MaxCounts",
                     "A_percent", "T_percent", "C_percent", "G_percent", "N_percent", "MaxPercent"]]

    final_df = pd.concat([Lconcat, seed_df], ignore_index=True)
    final_df = pd.concat([final_df, concat], ignore_index=True)

    final_df.reset_index(inplace=True)
    final_df = final_df.rename(columns = {'index':'Position'})

    final_df.to_csv(outdir + "/Final_assembled_info.csv", sep='\t', index=False)

    final_df["TotalCounts"] = final_df["A_counts"] + final_df["T_counts"] + final_df["C_counts"] + final_df["G_counts"] + final_df["N_counts"]
    final_df = final_df.transform(pd.to_numeric, errors='coerce')
    final_df["log10TotalCounts"] = np.log10(final_df["TotalCounts"])

    # Plot and Save Figures
    plot_figs(final_df, outdir)

    # Remove all intermediate files
    rmcmd = 'rm TEMP*'
    os.system(rmcmd)

    Endtime = datetime.now()
    print("Run has been completed at "+str(Endtime))
    print("Total run time: "+str(Endtime-Starttime))

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str,
                        default='/Users/sameersundrani/Desktop/FireLab/VaccineRNAPaperCode/1_Mod500pgMS0pg_S1_L001_R1_001_final_combined.fastq',
                        help='input fastq sequence (preprocessed or not) to assemble')
    parser.add_argument('--seedseq', type=str,
                        default='GTGCTTCACCAACGTGTACGCCGACAGCTTCGTG',
                        help='seed sequence from where we will begin our search (e.g. GTGCTTCACCAACGTGTACGCCGACAGCTTCGTG)')
    parser.add_argument('--nmer', type=int,
                        default=25,
                        help='size of kmer to complete searching rounds (default = 25)')
    parser.add_argument('--maxround', type=int,
                        default=65,
                        help='maximum number of searching rounds (default = 65)')
    parser.add_argument('--outdir', type=str,
                        default='.',
                        help='path to directory to save output files (default is cwd)')
    parser.add_argument('--process_multiple_fastq', type=bool,
                        default=False,
                        help='True if input files need to be preprocessed. Must input multiple fastq files, unprocesed separated by a space')

    args = parser.parse_args()
    main(args.input, args.seedseq, args.nmer, args.maxround, args.outdir, args.process_multiple_fastq)


