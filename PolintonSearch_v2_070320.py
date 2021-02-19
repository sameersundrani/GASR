'''
Dae-Eun Jeong and Sameer Sundrani
Perform a combinatoric search for various polinton proteins
'''

import os
import glob #add glob import for file lookups
import pandas as pd
import numpy as np
import itertools as it
import csv
import argparse

#GLOBAL VARS
directory = "/Users/sameersundrani/nematode_all_Feb18/ncbi_dataset/data"
ends_with_genome = "_genomic.fna"
speciesname = "C.elegans_VC2010"
# genomefile = "c_elegans_test.fa" (genome original to test on from Dae-Eun)
# here is our query
query = "/Users/sameersundrani/ConservedProteins_WBTransposon00000832.fa" #C. briggsae polinton aa seq


##STEP 0 IS TO MAKE SURE PATH VARIABLE HAS THE NCBI  ->  export PATH=$PATH:$HOME/ncbi-blast-2.10.1+/bin
# view using echo $PATH

##Step 1. tblastn

##Need to loop through and get all of our proteins

#lets get the folder paths named myDirs
def getFolderPaths(directory):
    for root, dirs, files in os.walk(directory):
        if (len(dirs) > 0): # may be some empty directories
            myDirs = [directory + os.sep + str(file) for file in dirs] #get a list of folder paths
    return myDirs

#The following code generates our _genomic.fna files to run our Polinton search through

def mergeGenomeReads(myDirs):
    for curDir in myDirs:
        #remove mistakes (because I messed up the first time)
        # if any(name.endswith("_genome.fna") or name.endswith("combinedGenome_genomic.fna") for name in os.listdir(curDir)):
        #     for fileToRem in os.listdir(curDir):
        #         if (fileToRem.endswith("_genome.fna") or fileToRem.endswith("combinedGenome_genomic.fna")):
        #             pathToRem = curDir + os.sep + fileToRem
        #             os.remove(pathToRem) #rem the bad file
        #             print("I removed " + fileToRem)

        if not any(name.endswith(ends_with_genome) for name in os.listdir(curDir)):
            curListToCat = [nameOfFile for nameOfFile in os.listdir(curDir) if (nameOfFile.endswith(".fna"))]
            commandToExecute = "cat"
            for fnaFile in curListToCat:
                commandToExecute += " " + curDir + os.sep + fnaFile
            last_part = os.path.basename(os.path.normpath(curDir))
            commandToExecute += " > " + curDir + "/" + last_part + "_genomic.fna" #solves the issue, prints the correct path name
            # print(commandToExecute)
            os.system(commandToExecute) #uncomment to run this command

    #quick check to see if we have the genomic files in the folders
    for curDir in myDirs:
        #remove mistakes
        if (len([nameOfFile for nameOfFile in os.listdir(curDir) if (nameOfFile.endswith("_genomic.fna"))]) >= 1):
            print("I have more than 1 genomic file here \n")

    # breakpoint() #break here for testing


# before we look for tBlastN matches let's save all our file names to a file
def saveOrigAssemblies(myDirs, logPath):
    listAssemblyNames = []
    for curDir in myDirs:
        listAssemblyNames.append(os.path.basename(os.path.normpath(curDir)))
    assemblyNamePath = logPath + 'assemblyList_nematodes_NCBI.txt'
    with open(assemblyNamePath, 'w') as filehandle:
        for item in listAssemblyNames:
            filehandle.write('%s\n' % item)
# done saving file


## Loop through ALL of our files and find polinton candidates
def searchFiles(myDirs, logDir, polFolder, nProts): #return num of pol candidates
    #init our finds to see how many have a >0 size bed file
    PolintonPtns = ["PRO-1_WBTransposon00000832_CBp", "POLB-1_WBTransposon00000832_CBp",
                    "ATP-1_WBTransposon00000832_CBp", "PY-1_WBTransposon00000832_CBp",
                    "INT-1_WBTransposon00000832_CBp"]
    polProteinAbbrevs = ['PRO', 'POLB-1', 'ATP-1', 'PY-1', 'INT-1']
    listCombinations = list(it.combinations([0, 1, 2, 3, 4], nProts)) #all combs of length nProts by index
    totalPolintonCandidatesDict = {}

    # # try first without PRO-1, need to edit feb 18 to be modular
    # Ptn_1 = PolintonPtns[0]
    # Ptn_2 = PolintonPtns[1]
    # Ptn_3 = PolintonPtns[2]
    # Ptn_4 = PolintonPtns[3]
    # Ptn_5 = PolintonPtns[4]
    if not os.path.isdir(logDir + polFolder):
        os.mkdir(logDir + polFolder) #build results directory

    for protComb in listCombinations:
        pathComb = '_'.join([polProteinAbbrevs[i] for i in protComb])
        myCandidateFinds = 0
        os.mkdir(logDir + polFolder + os.sep + pathComb)
        for curDir in myDirs:
            #create and go to a working directory
            # listFiles = [name for name in os.listdir(curDir) if name.endswith("_genome.fna") or name.endswith("combinedGenome_genomic.fna")]
            genomefile = [nameOfFile for nameOfFile in os.listdir(curDir) if (nameOfFile.endswith("_genomic.fna"))][0] #get our genome fasta file
            newpath = logDir + polFolder + os.sep + pathComb + os.sep + os.path.basename(os.path.normpath(curDir)) #make new path with the genomic file folder name

            #make directories in my main folder
            # try:
            #     os.mkdir(newpath)
            # except OSError:
            #     print("Creation of the directory %s failed" % newpath)
            # else:
            #     print("Successfully created the directory %s " % newpath)
            os.mkdir(newpath)
            print("Successfully created the directory %s " % newpath)
            os.chdir(newpath)
            makeblast_cmd = "makeblastdb -dbtype nucl -out Genome -in " + curDir + os.sep + genomefile
            os.system(makeblast_cmd)

            tblastn_cmd = "tblastn -db Genome -num_threads 2 -outfmt 6 -out tBLASTn_result.txt -query " + query
            os.system(tblastn_cmd)


            ##Step 2. tblastn result -> Bed file

            tBLASTn_result = pd.read_csv("tBLASTn_result.txt", sep="\t",
                                         names = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])

            sorted_result = tBLASTn_result.sort_values(["sseqid", "sstart"], ascending=[True, True]).reset_index(drop=True)
            sorted_result = sorted_result.groupby("sseqid").filter(lambda x: len(x) >= 5).reset_index(drop=True)
            sstart_col = sorted_result["sstart"]
            sstartendmin = []

            for i in range(len(sstart_col)):
                sstartendmin.append(sorted_result.iloc[i, 8:10].min())

            sorted_result["sstartendmin"] = sstartendmin
            contig_count = sorted_result.groupby("sseqid").count().reset_index()


            Index_distance = []
            bed_final_combined = pd.DataFrame()

            for contig in contig_count["sseqid"]:
                contig = sorted_result[sorted_result["sseqid"] == contig].reset_index(drop=True)

                Ptn_col = contig.qseqid

                cand_index_total = []

                for j in range(len(contig["sstart"])):
                    distance = contig.iloc[j, -1] - np.array(contig["sstartendmin"])
                    cand_distance = []
                    cand_index = []
                    klist = distance[j:j+15].tolist()

                    for k in it.takewhile(lambda x: abs(x) < 20000, klist):
                        cand_distance.append(k)

                    if len(cand_distance) >= 5:
                        for cand in range(j, j+len(cand_distance)):
                            cand_index.append(cand)

                        cand_index_ptn = []
                        for l in cand_index:
                            cand_index_ptn.append(Ptn_col[l])
                        # commented out PRO-1, need to reindent as needed
                        if (all(PolintonPtns[index] in cand_index_ptn for index in protComb)): #returns true only if all in comb are found
                            cand_index_total.append(cand_index)
                        # if Ptn_1 in cand_index_ptn:
                        # if Ptn_2 in cand_index_ptn:
                        #     if Ptn_3 in cand_index_ptn:
                        #         if Ptn_4 in cand_index_ptn:
                        #             if Ptn_5 in cand_index_ptn:
                        #                 cand_index_total.append(cand_index)
                        #             else:
                        #                 continue
                        #         else:
                        #             continue
                        #     else:
                        #         continue
                        # else:
                        #     continue
                        # else:
                        #     continue
                        #end changes
                    else:
                        continue

                cand_sseqid = []
                cand_sstart = []
                cand_send = []

                if len(cand_index_total) != 0:
                    for ci in range(len(cand_index_total)):
                        df_cand = contig.iloc[cand_index_total[ci]]
                        cand_sseqid.append(df_cand.iloc[0]["sseqid"])
                        cand_sstart.append(df_cand[["sstart", "send"]].min().min())
                        cand_send.append(df_cand[["sstart", "send"]].max().max())
                        df_bed_preliminary = pd.DataFrame(np.column_stack([cand_sseqid, cand_sstart, cand_send]), columns=['sseqid', 'sstart', 'send'])

                    df_bed_preliminary['send'] = df_bed_preliminary['send'].astype(str).astype(int)
                    df_bed_preliminary['sstart'] = df_bed_preliminary['sstart'].astype(str).astype(int)
                    df_bed_preliminary["length"] = df_bed_preliminary['send'] - df_bed_preliminary["sstart"]

                    df_bed_selected = df_bed_preliminary[df_bed_preliminary["length"] < 30000] # arbitrary
                    df_bed_selected2 = df_bed_selected[df_bed_selected["length"] > 5000]
                    df_bed_selected2["sstart-10kb"] = df_bed_selected2["sstart"] - 10000 # arbitrary choice for IRF
                    df_bed_selected2["send+10kb"] = df_bed_selected2["send"] + 10000
                    num = df_bed_selected2._get_numeric_data()
                    num[num<0] = 0
                    df_bed_final = df_bed_selected2[["sseqid", "sstart-10kb", "send+10kb"]]

                    bed_final_combined = pd.concat([bed_final_combined, df_bed_final], axis=0)

                else:
                    continue
            bedFileName = 'BEDcandidate_Polinton_' + pathComb + '.bed'
            bed_final_combined.to_csv(bedFileName, sep="\t", index=False, header=False)

            ##Step 3. parsing candidate polinton sequences
            filesize = os.path.getsize(bedFileName)
            if (filesize != 0): #only search the ones that have file sizes that have things in them
                myCandidateFinds += 1
                if pathComb not in totalPolintonCandidatesDict:
                    totalPolintonCandidatesDict[pathComb] = [os.path.basename(os.path.normpath(curDir))]
                else:
                    totalPolintonCandidatesDict[pathComb].append(os.path.basename(os.path.normpath(curDir)))
                print("We have now found: ", myCandidateFinds, " Polinton Candidates") #length == 222
                seq_parsing1 = 'seqkit subseq --bed' + bedFileName + ' ' + curDir + os.sep + genomefile + ' > candidates_polinton_' + pathComb + '.fa'
                os.system(seq_parsing1)
    return totalPolintonCandidatesDict
##Step 4. Inverted Repeat Finder (IRF) (on new file)
##NEED TO DO WINECONSOLE CMD IN CROSSOVER FOR THIS TO WORK
##transition to CX here

# breakpoint() ###set breakpoint for transition

# ##transition back

def main(nProts, genomeDir, logPath, polFolder):
    myPaths = getFolderPaths(genomeDir)
    # # print(myPaths)
    # # mergeGenomeReads(myPaths) only needed the first time, which I've already done with the nematodes (date =  Feb182021)
    # saveOrigAssemblies(myPaths, logPath)
    # For everything
    for i in range(2, 6): #already did 1
        d = searchFiles(myPaths, logPath, polFolder, i)
        print(d)
        df = pd.DataFrame.from_dict(d, orient='index')
        df.transpose().to_csv(logPath + 'polintonCandidates_proteinCountIs_' + str(i) + '.csv')


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--nprots', type=int, default = 5, help='input number of proteins you want to search for ex. 4)')
    parser.add_argument('--dir', type=str, default = '/Users/sameersundrani/nematode_all_Feb18/ncbi_dataset/data',
                        help='fasta original directory')
    parser.add_argument('--logPath', type=str, default = '/Users/sameersundrani/', help='path to save files to')
    parser.add_argument('--polFolder', type=str, default='PolintonSearchResultsFeb18', help='folder name to new things to')

    args = parser.parse_args()
    main(args.nprots, args.dir, args.logPath, args.polFolder)
