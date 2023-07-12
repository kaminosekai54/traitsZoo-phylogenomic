#!/bin/bash
#SBATCH --job-name=raxmlTest
#SBATCH --chdir=/scratch/aculpin/traitsZoo-phylogenomic/
#SBATCH --output=raxmlTest_log.o
#SBATCH --error=raxmlTest_log.e
#SBATCH --mail-user=alexis.culpin@cri-paris.org
#SBATCH --mail-type=none

#langages
. /softs/local/env/envperl-5.26.2.sh
. /softs/local/env/envpython-3.9.5.sh

#tools
. /softs/local/env/envdiamond-2.0.2.sh
. /softs/local/env/envmafft-7.475.sh

conda activate ~/pipeline
. /softs/local/env/envraxml-8.2.12.sh

python3 fastaRenamer.py
python3 getAndFilterBlast.py
python3 getAndFilterRBH.py
python3 getAlignment.py
python3 getTree.py
