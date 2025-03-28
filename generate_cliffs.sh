#!/bin/bash
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -p amd_256
source /public1/soft/modules/module.sh 
module load miniforge
source activate trans
export PYTHONUNBUFFERED=1
# python generate_cliffs.py -s 0.85

python generate_cliffs.py -d ./data/grampa_e_coli_7_25.csv -c "blosum62 average" -f 2 -t 0.85
