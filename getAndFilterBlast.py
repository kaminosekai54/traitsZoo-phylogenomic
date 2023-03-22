# import  
import sys, os, platform, subprocess, re, time
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
        print(refFile)
        print(file)
        cmd = "diamond makedb --in " + refFile + " -d " +settings["path"]["indexedFiles"] + file
        cmdTest = ["diamond", "makedb", "--in", refFile, "-d", settings["path"]["indexedFiles"] + file]
        print(cmd)
        # with open(settings["path"]["indexedFiles"] + file + ".sh", "w") as bashFile : bashFile.write("#!/bin/bash \n" + cmd)
        with open(settings["path"]["indexedFiles"] + file + ".sh", "w") as bashFile : bashFile.write("#!/bin/bash \n pwd")
        right = subprocess.run(["chmod", "777", settings["path"]["indexedFiles"] + file + ".sh"], text=True)
        print(right)
        bash= subprocess.check_call(settings["path"]["indexedFiles"] + file + ".sh")
        print(bash)

        # child= subprocess.check_output(r"diamond help", text=True)
        # child= subprocess.check_output(cmd, text=True)
        t2 = time.time()
        print("indexation done for " + refFile+ " in " + str((t2-t1)/60) + " sec")
    else: print("no file found")


def indexAllFile():
    if not os.path.isdir(settings["path"]["data"]) : os.mkdir(settings["path"]["data"])
    if not os.path.isdir(settings["path"]["indexedFiles"]) : os.mkdir(settings["path"]["indexedFiles"])
    t1 = time.time()
    with Pool() as p:
        p.map(indexFile, [settings["path"]["renamedFasta"] +file for file in os.listdir(settings["path"]["renamedFasta"])])
    t2 = time.time()
    print("all indexation done in " + str((t2-t1)/60) + " sec")

if __name__ == '__main__':
    indexAllFile()