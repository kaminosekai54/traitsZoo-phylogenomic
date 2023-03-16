import os
rootDir = "traitsZoo-phylogenomic/"
fileNames = [file.replace(".pep","") for file in os.listdir(rootDir+ "rawData") if file.endswith(".pep")]

rule all:
    input:
        expand("traitsZoo-phylogenomic/rawData/{name}.pep", name=fileNames),
        expand("traitsZoo-phylogenomic/data/fasta_files/{name}.fasta", name=fileNames),
        expand("traitsZoo-phylogenomic/data/tbl_files/{name}.tsv", name=fileNames)

rule convertToTBLAndFasta:
    input:
        "traitsZoo-phylogenomic/rawData/{name}.pep"
    output:
        "traitsZoo-phylogenomic/data/fasta_files/{name}.fasta",
        "traitsZoo-phylogenomic/data/tbl_files/{name}.tsv" 
    shell:
        "python3 traitsZoo-phylogenomic/fastaRenamer.py {input}"