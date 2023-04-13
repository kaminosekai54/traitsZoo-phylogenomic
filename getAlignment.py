# import
import sys, os, subprocess, time
from multiprocessing import Pool
import settings
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
import numpy as np


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

        cmd = [settings["pathToTrimal"] + "trimal", "-in", alignmentFile, "-out", outputFile, "-gt", "0.9", "-cons", "25", "-w", "3", "-st", "0.001"]
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
    print("all alignment trimming done in " + str((t2-t1)/60) + " sec")


# function getDictFromFasta
# this function read a fasta and return a dict {seqId : seq}
# @param,
# @ fastaFile, file to read
def getDictFromFasta(fastaFile):
    return {record.id[:record.id.find("_")]: record.seq for record in SeqIO.parse(fastaFile, "fasta")}
# function evaluateTrimming,
# this function create metrics to evaluate the quality of an alignment file and it's trimming
# @param,
# @trimmedFile, path to the trimmed file
# @rawFile, path to the raw file
def evaluateTrimming(trimmedFile, rawFile):
    # print(trimmedFile)
    # print(rawFile)
    trimmingData = {"nbBlock_strict":0, "nbBlock_permisif":0}
    
    trimmedAln = AlignIO.read(trimmedFile, "fasta")
    rawAln = AlignIO.read(rawFile, "fasta")
    trimmingData |= {"rawLength": rawAln.get_alignment_length(), "trimmedLength": trimmedAln.get_alignment_length(), "nbAATrimmed": rawAln.get_alignment_length() - trimmedAln.get_alignment_length()}
    # print(trimmingData)
    trimmedDict = getDictFromFasta(trimmedFile)
    trimmedDf = pd.DataFrame.from_dict({"seqId": trimmedDict.keys(), "aligned_seq": trimmedDict.values()})
    log = ""
    
    log += " block information for a strick consencsus of 0.9 identical AA : \n"
    log += findBlockInMSA(trimmedDf, 0.9, 5, trimmingData , "strict")
    log += " block information for a more permisif consencsus of 0.5 identical AA : \n"
    log += findBlockInMSA(trimmedDf, 0.5, 5, trimmingData , "permisif")

    # Compute the gap percentage for each sequence
    # trimmedDf["gapCount"] = trimmedDf.iloc[:, 1].str.count("P")
    trimmedDf["gapCount"] = trimmedDf.aligned_seq.apply(lambda x : str(x).count("-"))
    trimmedDf["sequenceLength"] = trimmedDf.iloc[:, 1].str.len()
    trimmedDf["gapPercentage"] = trimmedDf["gapCount"] / trimmedDf["sequenceLength"] * 100
    trimmedDf["gapPercentage"] = trimmedDf["gapPercentage"] .round(2)
    # print(trimmedDf)

    # Find the IDs of sequences with gap percentage more than 2 standard deviations from the mean
    mean_gap_pct = trimmedDf["gapPercentage"].mean()
    std_gap_pct = trimmedDf["gapPercentage"].std()
    threshold = mean_gap_pct + 2 * std_gap_pct

    outlier_ids = trimmedDf.loc[trimmedDf["gapPercentage"] > threshold, "seqId"].tolist()
    trimmingData["meanGapePerc"] = round(mean_gap_pct, 2)
    trimmingData["gapePercOutliers"] = outlier_ids
    trimmingDataLog =""
    for k,v in trimmingData.items() : trimmingDataLog += f"{k} : {v} \n"


    log  = trimmingDataLog + "\n" + log

    with open(settings["path"]["alignmentsLogs"] + trimmedFile[trimmedFile.rfind("/") +1 :].replace(".fasta", ".log"), "w") as logFile: logFile.write(log)
    return (trimmingData, trimmingData["trimmedLength"] / trimmingData["rawLength"] * 100, trimmedFile)


# function findBlockInMSA,
# This function find block in a msa, compute their size and the length betwwen each block
# @param,
# @trimmedDf, the dataframe representing the msa
# @consensus_threshold , Percentage of identical amino acids needed for consensus
# @min_block_size, Minimum block size to be considered a consensus block 
# @trimmingData , dictionnary where to add the data, 
# @suffix, string to add to the keys
def findBlockInMSA(trimmedDf, consensus_threshold , min_block_size,  trimmingData , suffix):
    blockList = []
    currentBlock = []
    log = ""
    for i in range(len(trimmedDf.loc[0, 'aligned_seq'])):
        col = trimmedDf.loc[:, 'aligned_seq'].str[i]
        consensus_freq = col.value_counts(normalize=True).iloc[0]
        if consensus_freq >= consensus_threshold:
            currentBlock.append(i)
        else:
            if len(currentBlock) >= min_block_size:
                blockList.append(currentBlock)
                trimmingData["nbBlock_" + suffix] +=1

            currentBlock = []

    if len(currentBlock) >= min_block_size:
        blockList.append(currentBlock)
        trimmingData["nbBlock_"+ suffix] +=1

    for i, block in enumerate(blockList):
        colName = "block_" + suffix + str(i+1) + "_infos"
        trimmingData[colName] = {"length":len(block), "start":block[0], "end":block[-1]}
        log += str(trimmingData[colName]) + "\n"
        # print(trimmedDf.loc[:, 'aligned_seq'].str[block[0]:block[-1]+1])

        if i+1 < len(blockList):
            nextBlock = blockList[i+1]
            trimmingData["dist_block"+ suffix+ str(i+1)+"-block"+ str(i+2)] = nextBlock[0] - block[-1]
            log += f"Distance between block {i+1} and block {i+2} = {nextBlock[0] - block[-1]} with the following sequences : \n"
            for seq in trimmedDf.loc[:, 'aligned_seq'].str[block[-1]: nextBlock[0]].tolist():
                log += str(seq) + "\n"

    return log

def evalAllTrimming():
    if not os.path.isdir(settings["path"]["alignmentsLogs"]): os.makedirs(settings["path"]["alignmentsLogs"])
    couple = []
    trimmedFiles = os.listdir(settings["path"]["trimmedAlignments"])
    rawFiles = os.listdir(settings["path"]["rawAlignments"])
    for raw in rawFiles:
        prefix = raw[raw.rfind("/")+1:-6]
        for trimmed in trimmedFiles:
            if prefix in trimmed :
                couple.append((settings["path"]["trimmedAlignments"]+trimmed, settings["path"]["rawAlignments"]+raw))
    res = []

    with Pool() as p:
        res = p.starmap(evaluateTrimming, couple)
        # res.append(evaluateTrimming(trimmedFile, rawFile))

    keep0 = []
    keep25 = []
    keep50 = []
    keep75 = []
    keep100 = []
    speciesIssueCount = {species[:species .find("_")] :0 for species in os.listdir(settings["path"]["renamedFasta"])}
    finalSummary = {"index":[], "nbBlockStrickt":[], "nbBlockPermisif":[], "rawLength":[], "trimmedLength":[], "keeptSeqPerc":[], "meanGapePerc":[], "meanGapePercOutliers":[]}
    for data, perc, file in res:
        finalSummary["index"].append(file[file.rfind("/"):file.find("rbhs")])
        finalSummary["nbBlockStrickt"].append(data["nbBlock_strict"])
        finalSummary["nbBlockPermisif"].append(data["nbBlock_permisif"])
        finalSummary["rawLength"].append(data["rawLength"])
        finalSummary["trimmedLength"].append(data["trimmedLength"])
        finalSummary["meanGapePerc"].append(data["meanGapePerc"])
        finalSummary["meanGapePercOutliers"].append(data["gapePercOutliers"])
        finalSummary["keeptSeqPerc"].append(perc)
        for species in data["gapePercOutliers"]: speciesIssueCount[species] +=1


        if perc <= 25 : keep25.append(file)
        if perc <= 50 and perc > 25: keep50.append(file)
        if perc <= 75 and perc > 50: keep75.append(file)
        if perc == 0:
            keep0.append(file)
            print("protein completly removed ", file)
        if perc >= 100:
            print(f"protein completely keppt ", file) 
            keep100.append(file)
        

    print(f'nb protein with 25 % of length = {len(keep25 )}')
    print(f'nb protein with 50 % of length = {len(keep50)}')
    print(f'nb protein with 75 % of length = {len(keep75)}')
    print(f'nb protein with 100 % of length = {len(keep100)}')
    print(f'nb protein with 0 % of length = {len(keep0)}')
    finalSummaryDf = pd.DataFrame.from_dict(finalSummary)
    finalSummaryDf.to_csv(settings["path"]["alignmentsLogs"] +"finalTrimmingMetricSummary.csv", sep="\t", index=False)
    speciesIssueCountDf = pd.DataFrame.from_dict({"speci": speciesIssueCount.keys(), "nbTimeCauseIssue": speciesIssueCount.values()})
    speciesIssueCountDf= speciesIssueCountDf.sort_values(by='nbTimeCauseIssue', ascending=False)
    speciesIssueCountDf .to_csv(settings["path"]["alignmentsLogs"] +"speciesOutliersSummary.csv", sep="\t", index=False)

if __name__ == '__main__':
    # alignAll()
    # trimAllAlignment()
    evalAllTrimming()