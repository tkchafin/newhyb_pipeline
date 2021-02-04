#!/bin/bash

#set arguments
popmap="example.popmap"
phylip="example.phy"

#target populations
p1="PopA"
p2="PopB"
pH="HybA+HybB"

#data filtering 
missPop="0.5" #missing data allower per SNP per population (phylip2newhybrids.pl)
missGlobal="0.5" #missing data per SNP globally (phylip2newhybrids.pl)
topLoci=200 #how many top loci to retain 
fst_threshold="0.01" #Fst threshold to retain a locus (Weird & Cockerham Fst)
fst_type="Fst" #Options are gst, Gst, GGst, D, Fst

#other settings 
sim_type="GT" #options are GT = freqbasedsim.GTFreq; and AS = freqbasedsim_AlleleSample
num_sims=2 #number of independent simulations
num_reps=2 #of replicates per simulation 
burnin=500
mcmc=1000
sample_size_sims=100
cores=2 #how many cores you want to use for parallel newhybrid runs. Requires GNU parallel 
pop_prior="T" #set to T to have the simulations specify the pure individuals as such

#convergence checking
PofZ_threshold=0.1 #assignment threshold of Pures to F2 category
Prop_threshold=0.5 #proportion of Pures allowed to exceed F2 threshold

#paths
rpath=/usr/local/bin/Rscript #path to Rscript 
pipeline=/Users/tyler/programs/newhyb_pipeline #path to newhyb_pipeline directory
nh_path=/Users/tyler/programs/newhyb_pipeline/newhybrids/bin/OSX/newhybsng #path to newhybng binary

#note: you can also randomly sample a max of X SNPs using <-r #>
echo "Selecting populations from phylip file... (Logfile: out.phy2newhyb.log)"
$pipeline/phylip2newhybrids.pl -p $popmap -i $phylip -1 $p1 -2 $p2 -a $pH -n $missPop -N $missGlobal -P >phy2newhyb.log

echo "Creating new population map (out.popmap)..."
paste -d"\t" ind_order.txt cat_order.txt > out.popmap

#output phylip file will be called out.phy
#convert to structure format 
echo "Converting to structure format... (Logfile: out.phy2str.log)"
$pipeline/phylip2structure.pl -i out.phy -p out.popmap -o out.str > phy2str.log

#convert to genepop
#stole script from Steve 
echo "Creating genepop file using Steve's str2genepop.pl script.."
$pipeline/str2genepop.pl -m out.popmap -s out.str -o "out.gen"

#Get top X Fst loci that are above a threshold
echo "Getting top $topLoci above $fst_type threshold of $fst_threshold... (Full Fst log: out.diffCalc)"
$rpath --vanilla --no-save --no-restore $pipeline/get_top_fst.R $topLoci $fst_threshold $fst_type> out.getTopLoc.log 2>&1

#simulate hybrids 
echo "Simulating hybrid genotypes across $num_sims simulations, each with $sum_reps replicates..."
$rpath --vanilla --no-save --no-restore $pipeline/sim_hybrids.R $num_sims $num_reps $sim_type $sample_size_sims > out.freqBasedSim.log 2>&1

#prep simulations to run
mkdir -p "out.sims"
mv out*_S*R*_NH.txt out.sims/
for f in out.sims/out*_S*R*_NH.txt
do 
   new=`echo $f | sed "s/.txt/_individuals.txt/"`
   touch $new
   for grp in Pure_PopA Pure_PopB F1.out F2.out BC_PopA BC_PopB
   do
     for ((i=1; i<= $sample_size_sims; i++))
     do
       echo $grp"_"$i >> $new
     done
   done
done
if [ $pop_prior == "T" ]
then
  echo "Adding population priors to simulated datasets..."
  for f in `find "out.sims" -name "*_NH.txt"`
  do
    #head -5 $f > .temp.1
    #tail -n +6 $f | sed -e "s"
		cp $f temp
    for i in `seq 1 $sample_size_sims`
    do
       sed -ie "s/^$i /$i z0 /" temp
    done
    for i in `seq $((sample_size_sims + 1)) $((sample_size_sims + sample_size_sims))`
    do
       sed -ie "s/^$i /$i z1 /" temp
    done
    mv temp $f
    #exit 0
  done
fi

#run newhybrids on the simulations
echo "Running NewHybrids on simulated datasets with $burnin and $mcmc sweeps..."
#$rpath --vanilla --no-save --no-restore $pipeline/run_newhyb.R $nh_path "out.sims" $burnin $mcmc > out.runNHsims.log 2>&1
#doing this myself because parallelnewhybrid doesn't work.
sim_dir="out.sims.NH.Results"
mkdir -p $sim_dir
run_nh (){
	base_name=`echo $1 | sed "s|.*/||" | sed "s/_NH.txt//"`
	sim_num=`echo $base_name | awk 'BEGIN{FS="_"}{print $2}' | awk 'BEGIN{FS=""}{print $2}'`
	rep_num=`echo $base_name | awk 'BEGIN{FS="_"}{print $2}' | awk 'BEGIN{FS=""}{print $4}'`
	seed1=$((rep_num + $RANDOM / sim_num))
	seed2=$((rep_num * $RANDOM - sim_num))
	dir_name=$base_name"_Results"
	mkdir -p $3/$dir_name
	cd $3/$dir_name
	#echo $seed1 "  " $seed2
	$2 -d ../../$1 --burn-in $4 --num-sweeps $5 -s $seed1 $seed2 --no-gui > stdout.log 2>&1
	cd - 2>&1 >/dev/null
}
export -f run_nh
#find "out.sims" -name "*_NH.txt" | parallel -j 1 'run_nh {} /Users/tyler/programs/newhyb_pipeline/newhybrids/bin/OSX/newhybsng out.sims.NH.Results 10 20'

for d in `find "out.sims" -name "*_NH.txt"`; do
  run_nh $d /Users/tyler/programs/newhyb_pipeline/newhybrids/bin/OSX/newhybsng out.sims.NH.Results $burnin $mcmc
done

#check simulations for convergence
echo "Checking simulated dataset runs for signs of non-convergence..."
to_get=$((sample_size_sims + sample_size_sims))
num_files=0
num_failed=0
for f in `find $sim_dir -name "aa-PofZ.txt"`
do
  num_files=$((num_files + 1))
  count=`tail -n +2 $f | head -n $to_get | awk -v thresh=$PofZ_threshold 'BEGIN{FS=" "; count=0}{if ($6 > thresh) count+=1}END{print count}'`
  prop=`echo $count / $to_get | bc -l`
	echo $prop "  " $Prop_threshold
  if (( $(echo "$prop > $Prop_threshold" | bc -l) ))
  then
    num_failed=$((num_failed + 1))
  fi
done
if [ $num_failed -ge 1 ]
then
  echo "$num_failed of $num_files simulations failed to converge"
  exit 0
fi

#

