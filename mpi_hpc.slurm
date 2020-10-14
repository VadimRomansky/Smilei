#!/bin/bash
#SBATCH --nodes=10
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH -p tornado
#SBATCH -J Smilei
#SBATCH -o pic-%j.out
#SBATCH -e pic-%j.err
#SBATCH -t 240:00:00
module load python/3.5.2 mpi/openmpi/3.0.0/gcc/7.2.0 library/hdf5/1.10.1/gcc72
mpirun -np 10 /home/ipntsr/romansky/Smilei/smilei sim2d.py
