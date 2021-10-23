#!/bin/sh
#SBATCH --time=96:00:00
#SBATCH --mem=100G
#SBATCH --account=USER

cd $ZERO

./Packages/lastz_D $ONE $TWO --notransition --step=20 --gap=400,30 --hspthresh=3000 --ydrop=9400 --gappedthresh=3000 --scores=./Sub_Scripts/HoxD55.q > $THREE
