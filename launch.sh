#!/bin/bash
#SBATCH --job-name=pipePhylo
#SBATCH --chdir=./
#SBATCH --output=pipePhylo_log.o
#SBATCH --error=pipePhylo_log.e
#SBATCH --mail-user=your_email@mon_domaine.org
#SBATCH --mail-type=none

#langages
. /softs/local/env/envperl-5.26.2.sh
. /softs/local/env/envpython-3.9.5.sh

#tools
conda activate pipelinePhylo

. /softs/local/env/envdiamond-2.0.2.sh
. /softs/local/env/envmafft-7.475.sh
. /softs/local/env/envraxml-8.2.12.sh

python3 fastaRenamer.py
python3 getAndFilterBlast.py
python3 getAndFilterRBH.py
python3 getAlignment.py
python3 getTree.py