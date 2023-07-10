#!/bin/bash
#SBATCH --job-name=pipPhylo
#SBATCH --chdir=/scratch/aculpin/traitsZoo-phylogenomic/
#SBATCH --output=pipeline_log.o
#SBATCH --error=pipeline_log.e
#SBATCH --mail-user=alexis.culpin@cri-paris.org
#SBATCH --mail-type=none

#langages
. /softs/local/env/envperl-5.26.2.sh
. /softs/local/env/envpython-3.9.5.sh

#tools
. /softs/local/env/envblast-2.9.0.sh
. /softs/local/env/envdiamond-2.0.2.sh
. /softs/local/env/envmafft-7.475.sh
. /softs/local/env/envraxml-8.2.12.sh
. /local/env/envconda.sh


conda activate ~/pipeline


python3 fastaRenamer.py
python3 getAndFilterBlast.py
python3 getAndFilterRBH.py
python3 getAlignment.py
python3 getTree.py

