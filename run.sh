#!/bin/bash

#SBATCH -n 24
#SBATCH -N 1
#SBATCH -p short
#SBATCH -t 0-7:00
#SBATCH -A angio
#SBATCH --export=ALL
#SBATCH -o out.txt
#SBATCH -e err.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andre.p.duarte@gmail.com

module load libs/mvapich2/2.1rc1-gcc-4.4.7
module load utils/taskfarmer
module load comp/gcc-4.9.2

mpirun -np 24 taskfarmer -f tasks.txt >sched.log
