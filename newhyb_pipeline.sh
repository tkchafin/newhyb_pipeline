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
topLoci="500" #how many top loci to retain 
fst_threshold="0.01" #Fst threshold to retain a locus (Weird & Cockerham Fst)
fst_type="Fst" #Options are gst, Gst, GGst, D, Fst

#other settings 
sim_type="AS" #options are GT = freqbasedsim.GTFreq; and AS = freqbasedsim_AlleleSample
num_sims="3" #number of independent simulations
num_reps="3" #of replicates per simulation 

#paths
rpath=/usr/local/bin/Rscript #path to Rscript 
pipeline=/Users/tyler/programs/newhyb_pipeline #path to newhyb_pipeline directory

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
$rpath --vanilla --no-save --no-restore $pipeline/sim_hybrids.R $num_sims $num_reps $sim_type > out.freqBasedSim.log 2>&1




