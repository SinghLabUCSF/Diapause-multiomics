#!/bin/sh
#SBATCH --time=80:00:00
#SBATCH --mem=200G
#SBATCH --account=USER

CURRENT_DIR=''

cd $CURRENT_DIR/Subscripts/

module load python3

python3 Trans_Novo.py $ONE

echo $ONE

