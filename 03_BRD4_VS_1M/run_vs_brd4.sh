#!/bin/bash
#PBS -N BRD4_1M_Screen
#PBS -q serial
#PBS -V
#PBS -j oe
#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=1:mem=8gb

cd $PBS_O_WORKDIR

echo "Job started at: $(date)"

module load anaconda
source activate brd4_env

python -u hpc_mega_screen.py

echo "Job finished at: $(date)"