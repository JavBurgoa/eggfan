#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --output=/g/arendt/Javier/Python/geneannotator/tests/genan.log
#SBATCH --mem-per-cpu=300MB
#SBATCH --time=4:00:00

srun python /g/arendt/Javier/Python/geneannotator/geneannotator/orthogroup.py
