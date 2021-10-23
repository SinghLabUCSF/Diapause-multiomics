#!/bin/sh
#SBATCH --time=8:00:00
#SBATCH --mem=20G
#SBATCH --account=USER

cd $ZERO

./Packages/axtChain -psl -minScore=2000 -linearGap=loose $THREE $ONE $TWO $FOUR -minScore=2000

