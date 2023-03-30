# import
import sys, os, subprocess, time
from multiprocessing import Pool
import settings

################################################################################
settings= settings.getSettings()

################################################################################
# functions
# def aligneSequenceWithMafft,
# This function use mafft to align a file
# with the --auto parameter only
def aligneSequenceWithMafft(fasta, outputFolder=settings["path"]["rawAlignments"]):
    if os.path.isfile(fasta):
        outputFile= outputFolder + fasta[fasta.rfind("/")+1:-6] + "_mafft_align.fasta"

        cmd = ["mafft", "--auto", "--out", outputFile, fasta]
        child= subprocess.check_output(cmd , text=True)
        # print(child)
    else: print("one of the file is missing")


# function alignAll,
# this function process all the alignment,
# in a parallele way
def alignAll():
    if not os.path.isdir(settings["path"]["data"]) : os.mkdir(settings["path"]["data"])
    if not os.path.isdir(settings["path"]["alignments"]) : os.mkdir(settings["path"]["alignments"])
    if not os.path.isdir(settings["path"]["rawAlignments"]) : os.mkdir(settings["path"]["rawAlignments"])
    t1 = time.time()
    with Pool() as p:
        p.map(aligneSequenceWithMafft, [settings["path"]["filteredRBH"]+ file for file in os.listdir(settings["path"]["filteredRBH"]) if file.endswith(".fasta")])
    t2 = time.time()
    print("all alignment done in " + str((t2-t1)/60) + " sec")

def trimAlignment(alignmentFile, outputFolder=settings["path"]["trimmedAlignments"]):
    if os.path.isfile(alignmentFile):
        outputFile= outputFolder + alignmentFile[alignmentFile.rfind("/")+1:-6] + "_trimmed.fasta"

        cmd = [settings["pathToTrimal"] + "trimal", "-in", alignmentFile, "-out", outputFile, "-noallgaps"]
        child= subprocess.check_output(cmd , text=True)
        # print(child)
    else: print("one of the file is missing")

def trimAllAlignment():
    if not os.path.isdir(settings["path"]["data"]) : os.mkdir(settings["path"]["data"])
    if not os.path.isdir(settings["path"]["alignments"]) : os.mkdir(settings["path"]["alignments"])
    if not os.path.isdir(settings["path"]["trimmedAlignments"]) : os.mkdir(settings["path"]["trimmedAlignments"])
    t1 = time.time()
    with Pool() as p:
        p.map(trimAlignment, [settings["path"]["rawAlignments"]+ file for file in os.listdir(settings["path"]["rawAlignments"]) if file.endswith(".fasta")])
    t2 = time.time()
    print("all alignment done in " + str((t2-t1)/60) + " sec")

if __name__ == '__main__':
    alignAll()
    trimAllAlignment()