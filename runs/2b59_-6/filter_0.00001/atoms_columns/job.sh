#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=128
#SBATCH -p fpga
#SBATCH -t 4:00:00
#SBATCH -A pc2-mitarbeiter


module reset
module load toolchain/intel/2024a
export OMP_NUM_THREADS=128
lscpu

numactl -b -m all ../../../../build/ptb_dev ../../../../systems/2b59.xyz -par ../../../../.atompara -chrg   -bas ../../../../.basis_vDZP -filter 0.00001 -check -purify > out 2> err
