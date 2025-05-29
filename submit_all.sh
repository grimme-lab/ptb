wd=`pwd`
for filter in 0.0001 0.00001 0.000001;
do
  #for col in full_matrix kmeans_r_heuristic atoms_columns;
  for col in full_matrix kmeans_r_heuristic;
  do
    cd $wd
    for sys in `ls systems/*.xyz | sed "s/\.xyz//g" | sed "s/systems\///g"`;
  #  for sys in 2b59 h2o_512;
  #  for sys in ADACAT;
    do
        cd $wd
        charge=`cat systems/$sys.charge`
        dir="runs/${sys}_$charge/filter_$filter/$col"
        echo $dir
        rm -rf $dir
        mkdir -p $dir
        cp job.sh $dir
        cd $dir
        sed -i "s/_SYS_/$sys/g" job.sh
        sed -i "s/_CHRG_/$chrg/g" job.sh
        sed -i "s/_FILTER_/$filter/g" job.sh
        echo "mode = submatrix" > .PUR
        echo "check" >> .PUR
        echo "submatrix_mode = submatrix_sygv" >> .PUR
        echo "submatrix_columns = $col " >> .PUR
        sbatch job.sh
    done  
  done
done
