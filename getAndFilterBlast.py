# import  
import sys, os, subprocess, time
import pandas as pd
from multiprocessing import Pool
import settings


################################################################################
settings= settings.getSettings()

################################################################################
# functions


# function indexFile,
# this function creat a bash script and run it to creat a diamond db file
# @param,
# @refFile, reference file to index
def indexFile(refFile):
    # if not os.path.isdir(settings["path"]["data"]) : os.mkdir(settings["path"]["data"])
    # if not os.path.isdir(settings["path"]["indexedFiles"]) : os.mkdir(settings["path"]["indexedFiles"])
    if os.path.isfile(refFile):
        t1 = time.time()
        file= refFile[refFile.rfind("/")+1: refFile.rfind(".")]
        # cmd = ["./diamond", "makedb", "--in", refFile, "-d", settings["path"]["indexedFiles"] + file]
        cmd = ["diamond", "makedb", "--in", refFile, "-d", settings["path"]["indexedFiles"] + file]
        # print(cmd)

        child= subprocess.check_output(cmd, text=True)
        t2 = time.time()
        # print("indexation done for " + refFile+ " in " + str((t2-t1)/60) + " sec")
    else: print("no file found")

# function indexAllFile,
# this function call the function indexFile on all file to index
# in a parallele way
def indexAllFile():
    if not os.path.isdir(settings["path"]["data"]) : os.mkdir(settings["path"]["data"])
    if not os.path.isdir(settings["path"]["indexedFiles"]) : os.mkdir(settings["path"]["indexedFiles"])
    t1 = time.time()
    with Pool() as p:
        p.map(indexFile, [settings["path"]["renamedFasta"] +file for file in os.listdir(settings["path"]["renamedFasta"])])
    t2 = time.time()
    print("all indexation done in " + str((t2-t1)/60) + " sec")

# function blastRefToOther,
# this function performe a "blast" search using diamond
# it will use the ref proteom as a querry file
# @ param
# @dbFile, database file to search match in
# @querryFile, the file to use as querry, by default, it's the renamed fasta reference proteom
def blastRefToOther(dbFile, querryFile= settings["path"]["renamedFasta"]+ settings["referenceProteom"] + "_pep.fasta"):
    if os.path.isfile(querryFile) and os.path.isfile(dbFile):
        if querryFile[querryFile.rfind("/") +1 :querryFile.rfind(".")] in dbFile: return
        t1 = time.time()
        outputFile=settings["path"]["diamondMatchs"] + querryFile[querryFile.rfind("/")+1: querryFile.rfind(".")] + "-matchs-in-" + dbFile[dbFile.rfind("/")+1: dbFile.rfind(".")] + ".tsv"
        # cmd = ["./diamond", "blastp", "--very-sensitive", "-q", querryFile, "-d", dbFile, "-o", outputFile, "--max-target-seqs", "1", "--outfmt", "6", "qseqid","sseqid","qlen","slen","length","ppos","pident","evalue","bitscore","full_qseq","full_sseq"]
        cmd = ["diamond", "blastp", "--very-sensitive", "-q", querryFile, "-d", dbFile, "-o", outputFile, "--max-target-seqs", "1", "--outfmt", "6", "qseqid","sseqid","qlen","slen","length","ppos","pident","evalue","bitscore","full_qseq","full_sseq"]
        # print(cmd)

        child= subprocess.check_output(cmd , text=True)
        print(child)
        t2 = time.time()
        # print("blast search completedfor " + querryFile + " in ", dbFile + " in " + str((t2-t1)/60) + " sec")
    else: print("no file found")


# function blastOtherToRef,
# this function performe a "blast" search using diamond
# it will use the ref proteom as a db filefile
# @ param
# @dbFile, database file to search match in, reference proteom by default
# @querryFile, the file to use as querry, 
def blastOtherToRef(querryFile, dbFile= settings["path"]["indexedFiles"]+ settings["referenceProteom"] + "_pep.dmnd"):
    if os.path.isfile(querryFile) and os.path.isfile(dbFile):
        
        if querryFile[querryFile.rfind("/") +1 :querryFile.rfind(".")] in dbFile: return
        t1 = time.time()
        outputFile=settings["path"]["diamondMatchs"] + querryFile[querryFile.rfind("/")+1: querryFile.rfind(".")] + "-matchs-in-" + dbFile[dbFile.rfind("/")+1: dbFile.rfind(".")] + ".tsv"
        # cmd = ["./diamond", "blastp", "--very-sensitive", "-q", querryFile, "-d", dbFile, "-o", outputFile, "--max-target-seqs", "1", "--outfmt", "6", "qseqid","sseqid","qlen","slen","length","ppos","pident","evalue","bitscore","full_qseq","full_sseq"]
        cmd = ["diamond", "blastp", "--very-sensitive", "-q", querryFile, "-d", dbFile, "-o", outputFile, "--max-target-seqs", "1", "--outfmt", "6", "qseqid","sseqid","qlen","slen","length","ppos","pident","evalue","bitscore","full_qseq","full_sseq"]
        # print(cmd)

        child= subprocess.check_output(cmd , text=True)
        t2 = time.time()
        # print("blast search completedfor " + querryFile + " in ", dbFile + " in " + str((t2-t1)/60) + " sec")
    else: print("one of the file is missing")

# function allBlastSearch,
# this function will perform the blast search from proteom to other file and the oposit
# in a parallele way
def allBlastSearch():
    if not os.path.isdir(settings["path"]["data"]) : os.mkdir(settings["path"]["data"])
    if not os.path.isdir(settings["path"]["diamondMatchs"]) : os.mkdir(settings["path"]["diamondMatchs"])
    t1 = time.time()
    with Pool() as p:
        p.map(blastRefToOther, [settings["path"]["indexedFiles"] + file for file in os.listdir(settings["path"]["indexedFiles"])])
        p.map(blastOtherToRef, [settings["path"]["renamedFasta"] + file for file in os.listdir(settings["path"]["renamedFasta"])])
    t2 = time.time()
    print("all blast search completed in " + str((t2-t1)/60) + " sec")



# main function
def main():
    if len(sys.argv) == 1 :
        indexAllFile()
        allBlastSearch()
    elif len(sys.argv) == 2:
        indexFile(sys.argv[1])
    elif len(sys.argv) == 3:
        blastRefToOther(sys.argv[2], sys.argv[1])
        blastOtherToRef(sys.argv[1], sys.argv[2])
    return

if __name__ == '__main__':
    main()