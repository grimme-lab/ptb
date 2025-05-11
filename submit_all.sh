filter=0.00001
wd=`pwd`
for col in full_matrix kmeans_r_heuristic atoms_columns single_columns;
do
  cd $wd
  for sys in `ls systems/*.xyz | sed "s/\.xyz//g" | sed "s/systems\///g"`;
  do
      cd $wd
      charge=`cat systems/$sys.charge`
      dir="runs/${sys}_$charge/filter_$filter/$col"
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
