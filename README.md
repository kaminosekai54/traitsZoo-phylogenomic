# traitsZoo-phylogenomic
A python pipeline to reconstruct copepods phylogenomic tree from transcriptoms / proteoms data.
This pipeline use different script for each step in order to allow you to re do every step individually.

A setting file is also provided  for small customisation.

The work flow of the pipeline is the following :


#Environement
You need to have several package installed,
- pandas
-biopython
- matplotlib
- seaborn
- numpy

You can also use the ** environement.yml ** provided like that :
`conda env create -f environment.yml`
 
You might need to change the environement name in the yml file to fit your own path.

# Usage :
Put your proteom file (fasta format) in the rawData folder.
The file should start by the species name :
Tigriopus-californicus.pep.fa
The first "_" or "." will be used later to find the species name. So make sure to name your file correctly.

For a quick and easy use simply run :
- `python fastaRenamer.py`
- `python getAndFilterBlast.py`
- `python getAndFilterRBH.py`
- `python getAlignment.py`
- `python getTree.py`

Those five commands will perform the all process.
Make sure you have choosen your "reference proteom" in the setting file.

⚠️ **Warning:**⚠️ 
Make sure to have [diamond](https://github.com/bbuchfink/diamond), [mafft](https://mafft.cbrc.jp/alignment/software/), and [RAxML](https://github.com/stamatak/standard-RAxML) and [TrimAL](https://github.com/inab/trimal) correctly installed.

# scripts :

** fastaRenamer.py **
This script is in charge of the data formating.
For each file, two new files will be created: a fasta and a tabular file. The name of the file corresponds to the species name followed by "_pep" (to indicate that it's amino acid and not nucleotide), and the file extension:
* Oithona-nana_female_ind17_Trinity.fasta.transdecoder.pep - > Oithona-nana_pep.fasta *

Inside the file, the sequences ID are renamed with the name of the species followed by an index:
* TRINITY_DN10000_c0_g1::TRINITY_DN10000_c0_g1_i1::g.33892::m.33892 -> Oithona-nana_1 *

When executing the file without argument, it will execute the renaming and file generation on all the files in the raw data folder. 
You can give it a file path as a command line argument to perform the formating on this specific file.
Additionnaly, you can give it another argument (1, 2) to choose to generate the fasta or tbl file only.

To summarise :

To format all the files in the "rawData" folder :
`python fastaRenamer.py`

To format a single file :
`python fastaRenamer.py my_file.fasta`

To create only the fasta file for a single file :
`python fastaRenamer.py my_file.fasta 1`
To create only the tbl file for a single file :
`python fastaRenamer.py my_file.fasta 2`


** getAndFilterBlast.py **
This script will perform the blast search.
It will index all the files in the "data/fastaFiles" folder, and then perform all blast search of our proteom against all other, and reciprocally.
You can also give it command line argument to perform only the indexing or blast search for an individual file.

To run all indexing and blast search :
`python getAndFilterBlast.py`

To index a fasta file :
`python getAndFilterBlast.py file_to_index.fasta`

To index and blast file1 vs file2 and file2 vs file1:
`python getAndFilterBlast.py file1.fasta file2.fasta`

** getAndFilterRBH.py **


This script will filter the diamond output file and compute the RBH.
It will first create couple for the diamond searches (ref-porteom1, proteom1-ref) and then apply filter on it (pid >= 30 %, evalue <1E-6, coverage > 60 % of the shortest sequences)
A final filter will be applied to keep only proteins with a RBH in all other species.


To perform the RBH identification and filtering :
`python getAndFilterRBH.py`

** getAlignments.py **
This script will align all RBH files using * MAFFT * Trim the aligned file using * TrimAL * and finally it will create bunch of log files and plot the alignments length distribution before and after triming.

To perform the alignment and triming : 
- `python getAlignment.py`

** getTree.py **
This script will concatenate all the trimed alignment files and then call * RAxML * to construct the phylogenomic tree.

to generate the phylogenomic tree :
`python getTree.py`