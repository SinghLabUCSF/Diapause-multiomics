#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --account=USER

cd $ZERO

./Packages/lavToPsl $ONE $TWO

