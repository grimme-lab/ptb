set -e
source /opt/intel/oneapi/setvars.sh
cd source
#make clean
make
cd ..
./build/ptb_dev systems/$1.xyz -par .atompara -chrg $2 -bas .basis_vDZP -filter $3 -check -purify $4 $5 $6 $7 $8 $9 | tee $1_$2_$3.out 
