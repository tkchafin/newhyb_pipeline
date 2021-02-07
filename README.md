# newhyb_pipeline

--work in progress--

## Installation 
You will need the hybriddetective, genepopedit, and parallelnewhybrid R packages. In R, install them like so:

devtools::install_github("rstanley/genepopedit")
devtools::install_github("bwringe/parallelnewhybrid")
devtools::install_github("bwringe/hybriddetective")
install.packages("diveRsity")

## Running the pipeline
You should be able to run the entire pipeline using ./newhyb_pipeline.sh. The options are currently hard-coded... You will need to specify various parameters and file names at the top of the newhyb_pipeline.sh script using your favorite text editor. 

The starting files that are required are a phylip-formatted sequence alignment (typically of concatenated SNPs), and a tab-delimited popmap:

```
Ind1	Pop1
Ind2	Pop2
Ind3	Pop1
...
...
```

You will specify the population IDs (identical to those in the popmap file) at the top of the newhyb_pipeline script:
```
p1="PopA"
p2="PopB"
pH="HybA+HybB"
```

You can combine population IDs for any of the P1 (=reference1), P2 (=reference2), or hybrid populations using "+". IMPORTANT: No spaces!!

The script will first extract your chosen samples to make a new phylip file (out.phy). It will then convert to create a genepop file using a script by Steve Mussmann (from github.com/stevemussmann/file_converters), which I've placed in this directory for convenience. An intermediate structure-formatted file will also be created (out.str)



