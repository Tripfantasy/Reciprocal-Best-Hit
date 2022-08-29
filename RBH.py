#!/usr/bin/env python3.10
import argparse
import os 

def get_args():
    parser = argparse.ArgumentParser(description = "Program to capture RBH of BLAST Data for 2 Specimen")
    parser.add_argument("-bf","--blastfiles", nargs = 2, help ="Input blast file names", type = str, required = True)
    parser.add_argument("-mf","--martfiles", nargs = 2, help ="Input biomart file names", type = str, required = False)

    return parser.parse_args()
args = get_args()

# Interactive environment command: srun --account=bgmp --partition=bgmp --nodes=1 --ntasks-per-node=1 --time=2:00:00 --cpus-per-task=1 --pty bash
# command to run: ./RBH.py -mf zfish_mart_rbh.txt human_mart_rbh.txt -bf HZ_rbh_set ZH_rbh_set 
#H_Z_blast = "/projects/bgmp/shared/Bi623/PS7_RBH_Bi623/H_to_zfishdb.blastp"
#Z_H_blast = "/projects/bgmp/shared/Bi623/PS7_RBH_Bi623/Z_to_homodb.blastp"
#-bf /projects/bgmp/dmarro/bioinfo/Bi623/assignments/RBH/HZ_rbh_set /projects/bgmp/dmarro/bioinfo/Bi623/assignments/RBH/ZH_rbh_set
#-mf /projects/bgmp/dmarro/bioinfo/Bi623/assignments/RBH/human_mart_rbh.txt /projects/bgmp/dmarro/bioinfo/Bi623/assignments/RBH/zfish_mart_rbh.txt

#Open mart files and define fields
z_mart = {}
with open(args.martfiles[0],'rt') as Zbiomart:
    for line in Zbiomart:
        column = line.strip("\n")
        column = line.split("\t")
        GeneID = column[0].strip()
        Gene_Name = column[1].strip()
        ProteinID = column[2].strip()
        z_mart[ProteinID]=[GeneID,Gene_Name]

h_mart = {} 
with open(args.martfiles[1],'rt') as Hbiomart:
    for line in Hbiomart:
        column = line.strip("\n")
        column = line.split("\t")
        GeneID = column[0].strip()
        Gene_Name = column[1].strip()
        ProteinID = column[2].strip()
        h_mart[ProteinID]=[GeneID,Gene_Name]

#Uncertain list will contain duplicate best-hit values to be omitted from HZ_dict.
uncertain_list1 = []
HZ_dict = {}
with open(args.blastfiles[0],'rt') as hz_set:
    for line in hz_set:
        column = line.strip("\n")
        columns = column.split(" ")
        H_ProteinID = columns[0]
        Z_ProteinID = columns[1]
        evalue = float(columns[2])
#Filters best hits per gene and places in HZ_dict.
        if H_ProteinID not in HZ_dict:
            HZ_dict[H_ProteinID]=[Z_ProteinID,evalue]
            continue
        if HZ_dict[H_ProteinID][1] == evalue:
            uncertain_list1.append(H_ProteinID)
            del HZ_dict[H_ProteinID]
#Makes sure uncertain hits are omitted post-filtering.
for entry in uncertain_list1:
    if entry in HZ_dict:
        del HZ_dict[entry]
print("Human to Zebra: ", len(HZ_dict))

#Uncertain list will contain omitted duplicate rbh for the same genes. (Same as prior chunk, but for ZH)
uncertain_list2 = []
ZH_dict = {}
with open(args.blastfiles[1], 'rt') as zh_set:
    for line in zh_set:
        column = line.strip("\t")
        columns = column.split(" ")
        Z_ProteinID = columns[0]
        H_ProteinID = columns[1]
        evalue = float(columns[2])

        if Z_ProteinID not in ZH_dict:
            ZH_dict[Z_ProteinID]=[H_ProteinID,evalue]
            continue
        if ZH_dict[Z_ProteinID][1] == evalue:
            uncertain_list2.append(Z_ProteinID)
            del ZH_dict[Z_ProteinID]

for entry in uncertain_list2:
    if entry in ZH_dict:
        del ZH_dict[entry]
print("Zebra to Human: ", len(ZH_dict))

#Filter reciprocal hits. From each direction's dictionary (Z:Z, H:Z) into matches{}
matches = {}
for Z_ProteinID in ZH_dict:
    MatchedHProt = ZH_dict[Z_ProteinID][0]

    if MatchedHProt in HZ_dict:
        MatchedZProt = HZ_dict[MatchedHProt][0]
    else:
        continue
    
    if HZ_dict[MatchedHProt][0] == Z_ProteinID:
        matches[MatchedHProt]=MatchedZProt
print("Reciprocal Best Hits: ",len(matches))
#Combine biomart info for gene name with the matches' info. (ID-MATCHID)
summary=[]
for H_ProteinID in matches:
    Z_ProteinID = matches[H_ProteinID]
    #print(z_mart[ProteinID])
    summary.append([Z_ProteinID]+z_mart[Z_ProteinID]+[H_ProteinID]+h_mart[H_ProteinID])
#Writes to summary file (ID,NAME,ID,ID,NAME,ID)
with open('Human_Zebrafish_RBH.tsv','wt') as summary_file:
    summary_file.write("Human Gene ID\tGene Name\tProtein ID\tZebra Gene ID\tGene Name\tProtein ID\n")
    for entries in summary:
        summary_file.write(entries[4]+"\t"+entries[5]+"\t"+entries[3]+"\t"+entries[1]+"\t"+entries[2]+"\t"+entries[0]+"\n")

