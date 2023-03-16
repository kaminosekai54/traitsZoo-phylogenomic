import os
rootDir = "traitsZoo-phylogenomic/"


fileNames = [file.replace(".pep","") for file in os.listdir(rootDir+ "rawData") if file.endswith(".pep")]
rule convertToTBLAndFasta:
    input:
        expand("{rootDir}rawData/{name}.pep", name=fileNames)
    output:
        expand("{rootDir}data/fasta_files/{name}.fasta", name=fileNames),
        expand("{rootDir}data/tbl_files/{name}.tsv", name=fileNames)
    shell:
        "python3 {rootDir}fastaRenamer.py {input}"























