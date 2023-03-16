#!/bin/bash
#SBATCH --job-name=pipPhylo
#SBATCH --chdir=/scratch/aculpin
#SBATCH --output=pipPhylo_log.o
#SBATCH --error=pipPhylo_log.e
#SBATCH --mail-user=alexis.culpin@cri-paris.org
#SBATCH --mail-type=all

#langages
. /softs/local/env/envperl-5.26.2.sh
. /softs/local/env/envopenjdk-17.0.3.sh
. /softs/local/env/envjava-1.8.0.sh
. /softs/local/env/envpython-3.7.1.sh

#tools
. /softs/local/env/envblast-2.9.0.sh
. /softs/local/env/envdiamond-2.0.2.sh
. /softs/local/env/envmuscle-3.8.1551.sh
. /softs/local/env/envmafft-7.475.sh
. /softs/local/env/envraxml_8.2.12.sh
./local/env/envconda.sh
. /local/env/envsnakemake-6.0.5.sh

#install Gblocks
#wget http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_Linux64_0.91b.tar.Z
#timeout et git hub n'existe plus.
#Si probleme persistant, install/utilisation de trimAI http://trimal.cgenomics.org/
#wget https://github.com/scapella/trimal/archive/v1.4.1.tar.gz


conda activate ~/pipeline
 # python3 traitsZoo-phylogenomic/fastaRenamer.py  rawData/Oithona-nana_female_ind17_Trinity.fasta.transdecoder.pep

snakemake -n -r -s traitsZoo-phylogenomic/rules.smk












