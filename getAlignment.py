# import
import sys, os, subprocess, time
from multiprocessing import Pool
import settings

################################################################################
settings= settings.getSettings()

################################################################################
# functions

def aligneSequenceWithMafft(fasta, outputLocation=settings["path"]["alignments"]):
    if os.path.isfile(fasta):
        outputFile = outputLocation + fasta[fasta.rfind("/")+1:-6]+ "_mafft_align.phy"
        tmpFile= outputLocation + fasta[fasta.rfind("/")+1:-6] + "_mafft_align.fasta"

        cmd = ["mafft", "--auto", "--out", tmpFile, fasta]
        child= subprocess.check_output(cmd , text=True)
        # print(child)
    else: print("one of the file is missing")



def alignAll():
    if not os.path.isdir(settings["path"]["data"]) : os.mkdir(settings["path"]["data"])
    if not os.path.isdir(settings["path"]["alignments"]) : os.mkdir(settings["path"]["alignments"])
    t1 = time.time()
    with Pool() as p:
        p.map(aligneSequenceWithMafft, [settings["path"]["filteredRBH"]+ file for file in os.listdir(settings["path"]["filteredRBH"]) if file.endswith(".fasta")])
    t2 = time.time()
    print("all alignment done in " + str((t2-t1)/60) + " sec")

if __name__ == '__main__':
    alignAll()