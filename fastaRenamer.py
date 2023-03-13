# imports
import os, sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


################################################################################
# functions


# function fastaToDict
# this function read a fasta and return a dictionnary where the keys are the sequence ID
# and the value the sequence
# this function also rename the sequence id according to the species name (indicated in the file) and an increasing index
# @param
# @fastaFile, path to the fasta to read
# @return, a dictionnary like {seqID:seq}
def fastaToDict(fastaFile):
    speciesName = fastaFile[fastaFile.rfind("/")+1: fastaFile.find("_")]
    return {speciesName+ "_" + str(i) : record.seq for i, record in enumerate(SeqIO.parse(fastaFile, "fasta"))}
# function fastaToTSV,
# this function transform a fasta file into a tsv file (seqID, sequence)
# @param,
# @fastaFile, path to the fasta to convert
# @fastaDict, dictionnary {seqID:sequence}, if provided, the function will use it instead of re reading the fasta
def fastaToTSV(fastaFile, fastaDict= {}):
    speciesName = fastaFile[fastaFile.rfind("/"): fastaFile.find("_")]
    dict = {}
    if len(fastaDict) == 0: dict = fastaToDict(fastaFile)
    # check if we have to read the fasta
    else : dict = fastaDict
    
    dict = {"seqId": dict.keys(), "sequence": dict.values()}
    df = pd.DataFrame.from_dict(dict)
    if not os.path.isdir("data/") : os.mkdir("data/")
    if not os.path.isdir("data/tbl_files") : os.mkdir("data/tbl_files")
    df.to_csv("data/tbl_files/" + speciesName + ".tab", sep= "\t", index=False)
    print("tbl convertion done for ", fastaFile)



# function renameSeqInFasta,
# this function write the the fasta file with the reanmed sequenceID
# @param,
# @fastaFile, path to the fasta to 
def renameSeqInFasta(fastaFile, fastaDict={}):
    speciesName = fastaFile[fastaFile.rfind("/"): fastaFile.find("_")]
    dict = {}
    if len(fastaDict) == 0 : dict = fastaToDict(fastaFile)
    else : dict = fastaDict
    
    # dict = {"seqId": dict.keys(), "sequence": dict.values()}
    df = pd.DataFrame(zip(dict.keys(), dict.values()))
    if not os.path.isdir("data/") : os.mkdir("data/")
    if not os.path.isdir("data/fasta_files") : os.mkdir("data/fasta_files")
    df.to_csv("data/fasta_files/" + speciesName + ".tab", sep= "\t", index=False)
    SeqIO.convert("data/fasta_files/" + speciesName + ".tab", 'tab', "data/fasta_files/" + speciesName + ".fasta", 'fasta-2line')
    os.remove("data/fasta_files/" + speciesName + ".tab")
    print("sequence id renamed for ", speciesName)

# main function
def main():
    if len(sys.argv) == 1 : print("error: you need to provide a file")
    elif len(sys.argv) == 2:
        fastaDict = fastaToDict(sys.argv[1])
        renameSeqInFasta(sys.argv[1], fastaDict)
        fastaToTSV(sys.argv[1], fastaDict)
    elif len(sys.argv) == 3:
        if int(sys.argv[2]) == 1: renameSeqInFasta(sys.argv[1])
        elif int(sys.argv[2]) == 2 : fastaToTSV(sys.argv[1])
        

if __name__ == '__main__':
    main()