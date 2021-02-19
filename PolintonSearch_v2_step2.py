import os
import glob #add glob import for file lookups
import pandas as pd
import numpy as np
import itertools as it

directory = "/Users/sameersundrani/FirstPartResultsOct22_noPRO"
ends_with_genome = "_genomic.fna"
query = "/Users/sameersundrani/ConservedProteins_WBTransposon00000832.fa"

for root, dirs, files in os.walk(directory):
    if (len(dirs) > 0): # may be some empty directories
        myDirs = [directory + os.sep + str(file) for file in dirs] #get a list of folder paths that we created before

# print(myDirs)

candidateDirs = []
for potentialDir in myDirs:
    if (len([nameOfFile for nameOfFile in os.listdir(potentialDir) if (nameOfFile.endswith("polinton_noPRO.fa"))]) >= 1):
        candidateDirs.append(os.path.basename(os.path.normpath(potentialDir)))

#before we do the IRF let's save all our file names to a file
# with open('candidatePolintonList_Oct23.txt', 'w') as filehandle:
#     for item in candidateDirs:
#         filehandle.write('%s\n' % item)
#done saving file

# print([dir for dir in candidateDirs])
print("My candidates are: \n")
for name in candidateDirs:
    print(name) # Found 23 candidates in all 5 protein searches, Found __ in 4 protein (no PRO-1)
#
# breakpoint()

#done manually
# input on CrossOver: wineconsole cmd
# then: cd \Users\sameersundrani\FirstPartResultsOct18\
# inverted_repeat_finder = "irf307.dos.exe candidates_polinton.fa 2 3 5 80 10 40 500000 30000 -d -h -t4 74 -t5 493 -t7 30000"
# os.system(inverted_repeat_finder)
##transition back

##Step 5. IRF result -> Bed file

# have to change directories
print("\n Now I am parsing through the candidate files 1 by 1")
for candidate in candidateDirs:
    # if (candidate == 'GCA_002207785.1'): #this had a bad fasta file, need to check!
    #     continue
    #create and go to a working directory
    newpath = "/Users/sameersundrani/FirstPartResultsOct22_noPRO/" + candidate  #make new path with the genomic file folder name

    #change to the directory as necessary
    os.chdir(newpath)
    IRF_result = pd.read_csv("candidates_polinton_noPRO.fa.2.3.5.80.10.40.500000.30000.dat", sep=" ",
                             names = ["Left_start", "Left_end", "Left_length", "Right_start", "Right_end", "Right_length", "Loop", "Percent_matches", "Percent_indels", "Score", "Percent_AT", "Percent_GC", "Percent_ATpairs", "Percent_GCpairs", "Percent_GTpairs", "CenterTimes2", "AverageCenterTimes2", "LeftSeq", "RightSeq"])

    IRF_result.loc[(IRF_result["Left_start"] == "Sequence:"), "ContigName"] = IRF_result.loc[(IRF_result["Left_start"] == "Sequence:"), "Left_end"]
    IRF_result["ContigName"].ffill(inplace=True)
    IRF_result = IRF_result[pd.notnull(IRF_result["RightSeq"])].reset_index(drop=True)
    IRF_result = IRF_result[["ContigName", "Left_start", "Left_end", "Left_length", "Right_start", "Right_end", "Right_length", "Loop", "Percent_matches", "Percent_indels", "Score", "Percent_AT", "Percent_GC", "Percent_ATpairs", "Percent_GCpairs", "Percent_GTpairs", "CenterTimes2", "AverageCenterTimes2", "LeftSeq", "RightSeq"]]
    IRF_result.ContigName = IRF_result.ContigName.str.replace(' ', '')
    IRF_result2 = IRF_result.set_index("ContigName", drop=False)
    IRF_result2["Left_start"] = IRF_result2["Left_start"].astype(int)
    IRF_result2["Right_end"] = IRF_result2["Right_end"].astype(int)

    cand_IRF = IRF_result2[IRF_result2["Loop"] > 10000].reset_index(drop=True)

    cand_IRF.to_csv("IRFtables_candidate_Polinton_noPRO.csv", index = True, header = True)

    cand_Contigname = np.array(cand_IRF["ContigName"].tolist())
    cand_Leftstart = cand_IRF["Left_start"].tolist()
    cand_Leftstart2 = np.array(cand_Leftstart)
    cand_Leftstart2 = cand_Leftstart2 - 7 #indexing method
    cand_Rightend = cand_IRF["Right_end"].tolist()
    cand_Rightend2 = np.array(cand_Rightend)
    cand_Rightend2 = cand_Rightend2 + 6 #indexing method

    cand_BED = pd.DataFrame(data=[cand_Contigname, cand_Leftstart2, cand_Rightend2])
    num = cand_BED._get_numeric_data()
    num[num<0] = 0

    cand_BED.T.to_csv("IRFBED_candidate_Politon_noPRO.bed", sep="\t", index=False, header=False)

    ##Step 6. parsing candidate polinton sequences that have inverted repeat + 6 bp upstream and downstream

    seq_parsing2 = "seqkit subseq --bed IRFBED_candidate_Politon_noPRO.bed candidates_polinton_noPRO.fa > IRF_candidates_polinton_noPRO.fa"
    os.system(seq_parsing2)

    ##Step 7. testing whether candidate sequences have Target Site Duplication (TSD)

    seq_parsing_left = "cat IRF_candidates_polinton_noPRO.fa | seqkit subseq -r 1:6 > Left.txt"
    os.system(seq_parsing_left)
    seq_parsing_right = "cat IRF_candidates_polinton_noPRO.fa | seqkit subseq -r -6:-1 > Right.txt"
    os.system(seq_parsing_right)

    Left6bp = pd.read_csv("Left.txt", sep="\t", names=["Left_6bp"])
    Right6bp = pd.read_csv("Right.txt", sep="\t", names=["Right_6bp"])

    Left6bp.loc[(Left6bp["Left_6bp"].str.startswith(">")), "Contig"] = Left6bp.loc[(Left6bp["Left_6bp"].str.startswith(">")), "Left_6bp"]
    Left6bp["Contig"].ffill(inplace=True)
    Left6bp = Left6bp[~Left6bp.Left_6bp.str.contains(">")]
    Left6bp = Left6bp[["Contig", "Left_6bp"]]

    Right6bp.loc[(Right6bp["Right_6bp"].str.startswith(">")), "Contig"] = Right6bp.loc[(Right6bp["Right_6bp"].str.startswith(">")), "Right_6bp"]
    Right6bp["Contig"].ffill(inplace=True)
    Right6bp = Right6bp[~Right6bp.Right_6bp.str.contains(">")]
    Right6bp = Right6bp[["Contig", "Right_6bp"]]

    LeftRight6bp = Left6bp
    LeftRight6bp["Right_6bp"] = Right6bp["Right_6bp"]
    LeftRight6bp.Contig = LeftRight6bp.Contig.str.replace('>', '')

    LeftRight6bp["TSD_match"] = LeftRight6bp["Left_6bp"].str.upper() == LeftRight6bp["Right_6bp"].str.upper()
    LeftRight6bp.to_csv("Table_TSD_matched_Polinton-1_all_noPRO.csv", index=False)

    TSD_matched = LeftRight6bp[LeftRight6bp["TSD_match"] == True].reset_index(drop=True)
    TSD_matched_Contig = TSD_matched["Contig"].str.replace(" ", "")
    TSD_matched.to_csv("Table_TSD_matched_Polinton_noPRO.csv")
    TSD_matched_Contig.to_csv("Contigs_TSD_matched_Polinton_noPRO.csv", index=False)


    ##Step 8. parsing TSD-matched sequences

    seq_parsing_3 = "cat IRF_candidates_polinton_noPRO.fa | seqkit grep -f Contigs_TSD_matched_Polinton_noPRO.csv > TSDmatched_polinton_noPRO.fa"
    os.system(seq_parsing_3)


    ##Step 9. removing duplicate sequences

    rm_dup = "cat TSDmatched_polinton_noPRO.fa | seqkit rmdup -s -o rmdup_TSDmatched_polinton_noPRO.fa"
    os.system(rm_dup)

#DONE!

print("I finished all the steps!!!")

