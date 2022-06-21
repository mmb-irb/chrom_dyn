#!/bin/bash
#$ -N Jue_Chrom_Float_MC
#$ -q cpu.q
#$ -cwd
export PATH="$HOME/Programs/SCHNArP/bin":'/orozco/homes/pluto/jwalther/Programs/x3dna-v2.1/bin':$PATH
SCHNA_AR_P=$HOME/Programs/SCHNArP
export SCHNA_AR_P
cd /home/jwalther/Programs/Floating_part-DNA/src
./MC_Chromatin
