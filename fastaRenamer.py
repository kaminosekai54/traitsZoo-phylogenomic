# imports
import os, sys, time
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool
import settings

################################################################################
settings = settings.getSettings()
# seqIdDict = {}


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
    if "_" in fastaFile: speciesName = fastaFile[fastaFile.rfind("/")+1: fastaFile.find("_")]
    else: speciesName = fastaFile[fastaFile.rfind("/")+1: fastaFile.find(".")]
    # speciesName = fastaFile[fastaFile.rfind("/")+1: fastaFile.find("_")]
    dictFasta = {}
    seqIdDict = {}
    for i, record in enumerate(SeqIO.parse(fastaFile, "fasta")):
        dictFasta[speciesName+ "_" + str(i+1)] = record.seq
        seqIdDict[speciesName+ "_" + str(i+1)] = record.id
    # speciesName+ "_" + str(i+1) : record.seq for i, record in enumerate(SeqIO.parse(fastaFile, "fasta"))}
    return dictFasta, seqIdDict  


# function fastaToTSV,
# this function transform a fasta file into a tsv file (seqID, sequence)
# @param,
# @fastaFile, path to the fasta to convert
# @fastaDict, dictionnary {seqID:sequence}, if provided, the function will use it instead of re reading the fasta
def fastaToTSV(fastaFile, fastaDict= {}):
    if "_" in fastaFile: speciesName = fastaFile[fastaFile.rfind("/")+1: fastaFile.find("_")]
    else: speciesName = fastaFile[fastaFile.rfind("/")+1: fastaFile.find(".")]
    # speciesName = fastaFile[fastaFile.rfind("/"): fastaFile.find("_")]
    file= fastaFile[fastaFile.rfind("/"): fastaFile.rfind(".")]
    dict = {}
    if len(fastaDict) == 0: dict = fastaToDict(fastaFile)[0]
    else : dict = fastaDict
    
    dict = {"seqId": dict.keys(), "sequence": dict.values()}
    df = pd.DataFrame.from_dict(dict)
    df.to_csv(settings["path"]["tblFiles"] + speciesName + "_pep.tsv", sep= "\t", index=False)
    # print("tbl convertion done for ", fastaFile)

# function renameSeqInFasta,
# this function write the fasta file with the renamed sequenceID
# @param,
# @fastaFile, path to the fasta to 
def renameSeqInFasta(fastaFile, fastaDict={}):
    if "_" in fastaFile: speciesName = fastaFile[fastaFile.rfind("/")+1: fastaFile.find("_")]
    else: speciesName = fastaFile[fastaFile.rfind("/")+1: fastaFile.find(".")]
    # speciesName = fastaFile[fastaFile.rfind("/"): fastaFile.find("_")]
    file= fastaFile[fastaFile.rfind("/"): fastaFile.rfind(".")]
    dict = {}
    if len(fastaDict) == 0 : dict = fastaToDict(fastaFile)[0]
    else : dict = fastaDict
    
    df = pd.DataFrame(zip(dict.keys(), dict.values()))
    if not os.path.isdir(settings["path"]["data"]) : os.mkdir(settings["path"]["data"])
    if not os.path.isdir(settings["path"]["renamedFasta"]) : os.mkdir(settings["path"]["renamedFasta"])
    df.to_csv(settings["path"]["renamedFasta"] + speciesName + ".tab", sep= "\t", index=False, header=False)
    SeqIO.convert(settings["path"]["renamedFasta"]+ speciesName + ".tab", 'tab', settings["path"]["renamedFasta"]  + speciesName+ "_pep.fasta", 'fasta-2line')
    os.remove(settings["path"]["renamedFasta"] + speciesName + ".tab")
    # print("sequence id renamed for ", speciesName)



# function processFile,
# This function creat a tsv and renamed fasta for the given file,
# needed to be able to paralelise the code
# @param
# @fastaFile, the file to process
def processFile(fastaFile):
    fastaDict, seqIdDict = fastaToDict(fastaFile)
    renameSeqInFasta(fastaFile, fastaDict)
    fastaToTSV(fastaFile, fastaDict)
    return seqIdDict
        
# function convertAll,
# this function apply parallely the processFiles function on all the file found in the rawData folder
def convertAll():
    t1 = time.time()
    listDict = []
    with Pool() as p:
        listDict= p.map(processFile, [settings["path"]["rawData"] +file for file in os.listdir(settings["path"]["rawData"])])
    t2= time.time()
    seqIdDict = {key: value for d in listDict for key, value in d.items()}

    df = pd.DataFrame.from_dict({"newId":seqIdDict .keys(), "oldId": seqIdDict .values()})
    df.to_csv(settings["path"]["tblFiles"] + "renamedIdList.tsv", sep= "\t", index=False)
    print("all file renamed in " + str((t2-t1)/60) + " sec")
    

# main function
def main():
    if len(sys.argv) == 1 :
        convertAll() 
    elif len(sys.argv) == 2:
        fastaDict = fastaToDict(sys.argv[1])
        renameSeqInFasta(sys.argv[1], fastaDict)
        fastaToTSV(sys.argv[1], fastaDict)
    elif len(sys.argv) == 3:
        if int(sys.argv[2]) == 1: renameSeqInFasta(sys.argv[1])
        elif int(sys.argv[2]) == 2 : fastaToTSV(sys.argv[1])
    return

if __name__ == '__main__':
    main()