#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=128
#SBATCH -p fpga
#SBATCH -t 8:00:00
#SBATCH -A pc2-mitarbeiter


module reset
module load toolchain/intel/2024a
export OMP_NUM_THREADS=128

./build/ptb_dev systems/$1.xyz -par .atompara -chrg $2 -bas .basis_vDZP -filter $3 -check -purify >  outputs/$1_$2_$3.out 2> outputs/$1_$2_$3.err
