# import  
import sys, os, subprocess, time
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool
import psutil
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
    t1 = time.time()
    dictList = []
    concatDict= {}
    with Pool() as p:
        dictList =p.map(getDictFromFasta, [settings["path"]["trimmedAlignments"]+ file for file in os.listdir(settings["path"]["trimmedAlignments"]) if file.endswith(".fasta")])
    for d in dictList: 
        for k,v in d.items():
            if k in concatDict.keys():
                concatDict[k] += v
            else: concatDict[k] = v
    concatDf = pd.DataFrame.from_dict({"seqId": concatDict.keys(), "seq": concatDict.values()})
    concatDf.to_csv(settings["path"]["trees"] + "concatenatedRBHAlignments.tab", sep= "\t", index=False, header=False)
    SeqIO.convert(settings["path"]["trees"] + "concatenatedRBHAlignments.tab", 'tab', settings["path"]["trees"] + "concatenatedRBHAlignments.fasta", 'fasta-2line')
    os.remove(settings["path"]["trees"] + "concatenatedRBHAlignments.tab")
    print("number of species in the concatenated alignment :", len(concatDf))
    print("number of proteins in the concatenated alignment :", len(dictList))
    print("length of the concatenated alignment :", len(list(concatDict.values())[0]))
    t2 = time.time()
    print("concatenated alignment file generated in " + str((t2-t1)/60) + " sec")
    return settings["path"]["trees"] + "concatenatedRBHAlignments.fasta"

def generateTree(concatenatedAlignmentFile):
    if os.path.isfile(concatenatedAlignmentFile):
        outputFile= concatenatedAlignmentFile.replace("fasta", "tree")

        cpu_affinity = len(psutil.Process().cpu_affinity())
        cpu_count = psutil.cpu_count()
        # Check if the process has CPU affinity set
        if cpu_affinity > 0:
        # Limit the CPU affinity count to the available CPU count
            cpu_affinity = min(cpu_affinity, cpu_count)
        else:
            cpu_affinity = cpu_count
        
# Convert the final count to a string
        cpu_affinity_str = str(cpu_affinity)
        print("number of cpu provided : ", cpu_affinity)

        cmd = ["raxmlHPC-PTHREADS", "-T", cpu_affinity_str, "-f", "a", "-m", "PROTGAMMAAUTO", "-p", "12345", "-x", "12345", "-#", "100", "-s", concatenatedAlignmentFile[concatenatedAlignmentFile.rfind("/")+1:], "-n", outputFile[outputFile.rfind("/")+1:]]
        if cpu_affinity  < 4 :
            print("not enough cpu provided to use the more efficiant PTHREADS version, will use the single raxmlHPC version ")
            cmd = ["raxmlHPC", "-f", "a", "-m", "PROTGAMMAAUTO", "-p", "12345", "-x", "12345", "-#", "100", "-s", concatenatedAlignmentFile[concatenatedAlignmentFile.rfind("/")+1:], "-n", outputFile[outputFile.rfind("/")+1:]]
            
        print(" ".join(cmd))
        wd = os.getcwd()
        print(os.getcwd())
        os.chdir(settings["path"]["trees"])
        print(os.getcwd())
        

        child= subprocess.check_output(cmd , text=True)
        os.chdir(wd)
        print(child)
        print(os.getcwd())
    else: print("concatenated file not found : ", concatenatedAlignmentFile)

if __name__ == '__main__':
    concatenatedAlignmentFile = concatAllRBH()
    generateTree(concatenatedAlignmentFile)