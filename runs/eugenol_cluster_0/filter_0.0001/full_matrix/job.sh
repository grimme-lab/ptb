#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=64
#SBATCH -p fpga
#SBATCH -t 1-0
#SBATCH -A pc2-mitarbeiter

module reset
module load toolchain/intel/2024a
export OMP_NUM_THREADS=64
lscpu
export OMP_STACKSIZE=16G

numactl -i 0-3 --cpunodebind=0-3 ../../../../build/ptb_dev ../../../../systems/eugenol_cluster.xyz -par ../../../../.atompara -chrg   -bas ../../../../.basis_vDZP -filter 0.0001 -check -purify > out 2> err
#numactl -m 0 --cpunodebind=0 ../../../../build/ptb_dev ../../../../systems/eugenol_cluster.xyz -par ../../../../.atompara -chrg   -bas ../../../../.basis_vDZP -filter 0.0001 -check -purify > out 2> err
