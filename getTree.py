# import  
import sys, os, subprocess, time
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool
import settings
################################################################################
settings= settings.getSettings()

################################################################################
# functions

# function getDictFromFasta
# this cuntion read a fasta and return a dict {seqId : seq}
# @param,
# @ fastaFile, file to read
def getDictFromFasta(fastaFile):
    return {record.id[:record.id.find("_")]: record.seq for record in SeqIO.parse(fastaFile, "fasta")}


# function concatAllRBH
# this function read all RBH file into a dict and concatenate them 
# Finally it create a big concatenated fasta.
def concatAllRBH():
    if not os.path.isdir(settings["path"]["data"]) : os.mkdir(settings["path"]["data"])
    if not os.path.isdir(settings["path"]["trees"]) : os.mkdir(settings["path"]["trees"])
    t1 = time.time()
    dictList = []
    concatDict= {}
    with Pool() as p:
        dictList =p.map(getDictFromFasta, [settings["path"]["trimmedAlignments"]+ file for file in os.listdir(settings["path"]["trimmedAlignments"]) if file.endswith(".fasta")])
    for d in dictList: concatDict = concatDict | d
    concatDf = pd.DataFrame.from_dict({"seqId": concatDict.keys(), "seq": concatDict.values()})
    concatDf.to_csv(settings["path"]["trees"] + "concatenatedRBHAlignments.tab", sep= "\t", index=False, header=False)
    SeqIO.convert(settings["path"]["trees"] + "concatenatedRBHAlignments.tab", 'tab', settings["path"]["trees"] + "concatenatedRBHAlignments.fasta", 'fasta-2line')
    os.remove(settings["path"]["trees"] + "concatenatedRBHAlignments.tab")
    print("number of species in the concatenated alignment :", len(concatDf))
    print("number of proteins in the concatenated alignment :", len(dictList))
    print("length of the concatenated alignment :", len(concatDf.seq.tolist()[0]))
    t2 = time.time()
    print("concatenated alignment file generated in " + str((t2-t1)/60) + " sec")

def generateTree(concatenatedAlignmentFile):
    if os.path.isfile(concatenatedAlignmentFile):
        outputFile= concatenatedAlignmentFile.replace("fasta", "tree")

        cmd = ["raxmlHPC","-T","12","-f","a","-m","PROTGAMMAAUTO","-p","12345","-x","12345","-#","100","-s", concatenatedAlignmentFile,"-n", outputFile]

        child= subprocess.check_output(cmd , text=True)
        # print(child)
    else: print("concatenated file not found : ", concatenatedAlignmentFile)

if __name__ == '__main__':
    concatAllRBH()
    #generateTree()