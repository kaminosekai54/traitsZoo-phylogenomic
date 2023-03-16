import os
rootDir = "traitsZoo-phylogenomic/"
fileNames = [file.replace(".pep","") for file in os.listdir(rootDir+ "rawData") if file.endswith(".pep")]
print(fileNames)
rule convertToTBLAndFasta:
    input:
        expand("{rootDir}rawData/{name}.pep", name=fileNames, rootDir= rootDir)
    output:
        expand("{rootDir}data/fasta_files/{name}.fasta", name=fileNames, rootDir= rootDir),
        expand("{rootDir}data/tbl_files/{name}.tsv", name=fileNames, rootDir=rootDir)
    shell:
        "python3 traitsZoo-phylogenomic/fastaRenamer.py {input}"