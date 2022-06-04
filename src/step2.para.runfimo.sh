### bsub -q Z-ZQF -n 20 -eo MIIGV.rand.eo "./para.runfimo.sh MIIGV.rand.fa" 

PARA=/src/parallel  ### speficy a path to GNU parallel
FIMO=/meme_4.11.2-app/bin/fimo  ### specify a path to MEME suite tool "fimo"

dir=`pwd`
lib=$dir/DB_motif_collection.out   ### change the motif library in your run
file=$dir/$1  ### fasta file extracted from step1

rundir=$dir/${1%.fa}
mkdir -p $rundir

cp $lib $file $rundir
cd $rundir

split -d -a 3 -l 100 $file SPLIT

$PARA -j 20 "$FIMO --verbosity 1 --text --norc $lib {} >fimo.{}.out" ::: SPLIT*

cd ../
