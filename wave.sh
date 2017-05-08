#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 5:00
#SBATCH -p node
#SBATCH -n 1

mpirun -np 1 ./wave 1 1 400 100
