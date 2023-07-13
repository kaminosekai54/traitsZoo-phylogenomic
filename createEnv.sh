#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=createEnv
#SBATCH --output=create_env.o
#SBATCH --error=create_env.e

#langages
. /softs/local/env/envpython-3.9.5.sh

#tools
. /local/env/envconda.sh


conda env create -f environment.yml

conda activate pipelinePhylo



