filter=0.00001
mkdir outputs
for sys in `ls systems/*.xyz | sed "s/\.xyz//g" | sed "s/systems\///g"`;
do
    charge=`cat systems/$sys.charge`
    sbatch job.sh $sys $charge $filter
done
