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



# function getBRH,
# this function will found the BRH  between two diamond output file,
# the function also apply filter on identity > 30 %, and evalue < 1e-6
# @param
# @matchCouple, a tupple of diamond output file
# @return, return a pandas data frame containing the BRH fount in the couple of file
def getBRH(args):
    filter =False
    if len(args) == 3 : filter = args[2]
    df1 = pd.read_csv(args[0], sep="\t", names=["qseqid","sseqid","qlen","slen","length","ppos","pident","evalue","bitscore","full_qseq","full_sseq"])
    df2 = pd.read_csv(args[1], sep="\t", names=["qseqid","sseqid","qlen","slen","length","ppos","pident","evalue","bitscore","full_qseq","full_sseq"])
    df1 = df1.assign(brhCouple = df1.qseqid + "__" +df1.sseqid)
    df2 = df2.assign(brhCouple = df2.sseqid + "__" +df2.qseqid)
    dfBRH = df1[(df1.brhCouple.isin(df2.brhCouple.tolist()))]
    if filter :
        dfBRH = df1[(df1.brhCouple.isin(df2.brhCouple.tolist())) & (df1.pident >= 30) & (df1.evalue < 1e-5)]
        dfBRH = dfBRH.assign(minSeqLen = dfBRH[["qlen", "slen"]].min(axis=1))
        dfBRH = dfBRH[dfBRH.length >= dfBRH.minSeqLen *0.6 ]
    dfBRH = dfBRH[["qseqid","sseqid","full_qseq","full_sseq"]]
    # dfBRH=dfBRH.set_index(dfBRH.qseqid)
    
    for col in dfBRH.columns:
        if col not in  ["qseqid", "full_qseq"]:
            dfBRH[col] = dfBRH[col].apply(lambda x: x + ",")
    return dfBRH




# function mergeBrhDF,
# this function merge two brh dataframe outputed by getBRH function
# the subject id and sequence will be concatenated if a querry seq already exist, or be added if not
# @param,
# @df1, first data frame to merge
# @df2, second data frame to merge
# @colNotToMerge, a list of column name not to merge
# @return, a merged data frame
def mergeBrhDF(df1, df2, colNotToMerge):
    df = df1.merge(df2, how="outer", on=colNotToMerge)
    df = df.fillna("")
    for col in df1.columns:
        if col  not in colNotToMerge:
            name =  col +"_final"
            df.insert(len(df.columns), name, df[col+"_x"] + df[col + "_y"])
    colToReturn = colNotToMerge + [col for col in df.columns if col not in  colNotToMerge and "final" in col]
    df = df[colToReturn]
    renameDict = {col: col.replace("_final","") for col in colToReturn }
    df = df.rename(renameDict, axis=1)
    return df

#  function finalFilter,
# this function apply the final filter on the brh data frame
# this function will get only the gene that have a rbh in all other species (n-1 brh) 
# n being the number of species in our dataset
# @param
# @df, the dataframe on whish to apply the filter,
# @colNotToMerge, a list of column name not to merge
# @return the filtered data frame
def finalFilter(df, colNotToMerge, outputFolder = settings["path"]["filteredRBH"]):
    nbRBHRequiered = len(os.listdir(settings["path"]["renamedFasta"])) -1
    for col in df.columns:
        if col not in colNotToMerge:
            df[col] = df[col].apply(lambda x :[s.replace(",", "") for s in x.split(",") if s != ""])
    df = df.assign(nbRBH = df.sseqid.apply(lambda x : len(x)))
    df = df[df.nbRBH == nbRBHRequiered]

    if not os.path.isdir(settings["path"]["data"]) : os.mkdir(settings["path"]["data"])
    if not os.path.isdir(settings["path"]["RBH"]) : os.mkdir(settings["path"]["RBH"])
    if not os.path.isdir(outputFolder) : os.mkdir(outputFolder)
    df.to_csv(outputFolder+ "full_rbh_tab.tsv", sep="\t", index=False)
    return df[["qseqid","sseqid","full_qseq","full_sseq"]]

# function createRBHFile,
# this function creat a fasta file containing a sequence of the reference proteom and all its BRH.
# @param,
# @row, a dictionnary like object, basically a row of the final rbh data frame
def createRBHFile(args):
    row = args[0]
    outputFolder = args[1]
    t1= time.time()
    reccordList = [SeqRecord(Seq(row["full_qseq"]), id = row["qseqid"], description="")]
    for seqId, sSeq in zip(row["sseqid"], row["full_sseq"]):
        reccordList.append(SeqRecord(Seq(sSeq), id = seqId, description=""))

    if not os.path.isdir(settings["path"]["data"]) : os.mkdir(settings["path"]["data"])
    if not os.path.isdir(settings["path"]["RBH"]) : os.mkdir(settings["path"]["RBH"])
    if not os.path.isdir(outputFolder) : os.mkdir(outputFolder)
    fileName = row["qseqid"] + "_rbhs_pep.fasta" 
    SeqIO.write(reccordList, outputFolder+ fileName, "fasta-2line")
    t2 = time.time()
    return fileName

# function createAllRBHFile,
# this function apply the createRBHFile on each row of the data frame
# in a parallele way
# @param
# @df, the data frame on wish to call the function
def createAllRBHFile(df, outputFolder = settings["path"]["filteredRBH"]):
    t1 = time.time()
    fileList = []
    with Pool() as p:
        fileList = p.map(createRBHFile, [(v, outputFolder) for k,v in df.iterrows()])
    t2 = time.time()
    print(str(len(fileList)) + " rbh files created in " + str((t2-t1)/60) + " sec")


# function getAllBRH
# this function create a fasta file for each brh fount,
# apply parallele computation when possible
def getAllBRH():
    coupleList= []
    fileList = os.listdir(settings["path"]["diamondMatchs"])
    for file in fileList :
        s1 = file[:file.find("-matchs-in-")]
        s2 = file[file.find("-matchs-in-") + len("-matchs-in-"): file.rfind(".")]
        if s2 + "-matchs-in-"+ s1 +".tsv" in fileList:
            if settings["referenceProteom"] in s1:
                couple = (file, s2 + "-matchs-in-"+ s1 +".tsv")
            elif settings["referenceProteom"] in s2:
                couple = (s2 + "-matchs-in-"+ s1 +".tsv", file)
            if not couple in coupleList and not (couple[1], couple[0]) in coupleList: coupleList.append(couple)
    finalRawRBHDf= pd.DataFrame.from_dict({"qseqid":"",  "sseqid":[],  "full_qseq":[], "full_sseq":[]})
    finalFilteredRBHDf= pd.DataFrame.from_dict({"qseqid":"",  "sseqid":[],  "full_qseq":[], "full_sseq":[]})
    colNotToMerge = ["qseqid","full_qseq"]
    rawRBHDfList = []
    filteredRBHDfList = []
    with Pool() as p:
        rawRBHDfList = p.map(getBRH, [(settings["path"]["diamondMatchs"] + c1, settings["path"]["diamondMatchs"] +c2)for c1, c2 in coupleList])
        filteredRBHDfList = p.map(getBRH, [(settings["path"]["diamondMatchs"] + c1, settings["path"]["diamondMatchs"] +c2, True) for c1, c2 in coupleList])
    for rbhDf in rawRBHDfList:
        finalRawRBHDf= mergeBrhDF(finalRawRBHDf, rbhDf, colNotToMerge)

    for rbhDf in filteredRBHDfList:
        finalFilteredRBHDf = mergeBrhDF(finalFilteredRBHDf, rbhDf, colNotToMerge)
    finalRawRBHDf= finalFilter(finalRawRBHDf, colNotToMerge, outputFolder=settings["path"]["rawRBH"] )
    finalFilteredRBHDf = finalFilter(finalFilteredRBHDf, colNotToMerge)
    # print(finalFilteredRBHDf)
    createAllRBHFile(finalRawRBHDf, settings["path"]["rawRBH"] )
    createAllRBHFile(finalFilteredRBHDf, settings["path"]["filteredRBH"])

if __name__ == '__main__':
    getAllBRH()