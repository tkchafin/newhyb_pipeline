#devtools::install_github("rystanley/genepopedit")
#devtools::install_github("bwringe/hybriddetective")
#devtools::install_github("bwringe/parallelnewhybrid")
#modified from script by Bradley Martin
suppressMessages(library("parallelnewhybrid"))
suppressMessages(library("hybriddetective"))
suppressMessages(library("genepopedit"))
suppressMessages(library("diveRsity"))
# suppressMessages(library("genepop"))
##*** SET WORKING DIRECTORY HERE ***###

# Store working directory as object.
my_path <- getwd()

###############################################################################

args = commandArgs(trailingOnly=TRUE)

# Name of GENEPOP input file
GenePopData <- "out.gen"
#print(args)
top_n_loci <- as.numeric(args[1])
fst.threshold <- as.numeric(args[2])
fst.type <- args[3]

# Filter genepop to only the pure populations 
subset_genepop(GenePopData, subs=NULL, spop=c("PureA", "PureB"), path="out.pure.gen")

#Calculate FST
#writeLines("Calculating locus-specific Fst\n\n")

FST.dat <- diveRsity::diffCalc(infile = "out.pure.gen", outfile = "out.diffCalc", fst = T, bs_locus = F)
FST.dat <- FST.dat$std_stats[1:length(FST.dat$std_stats$loci)-1,]
FSTs <- FST.dat[,fst.type]

FST.df <- data.frame(genepopedit::genepop_detective("out.pure.gen",variable="Loci"), FSTs)
names(FST.df)[1] <- "loci"

#sort dataframe
FST.df <- FST.df[base::order(FST.df$FSTs, decreasing = TRUE),]

#get Fst's above threshold 
FST.Filter.Vec <- as.character(FST.df[which(FST.df$FSTs >=fst.threshold), 1])

# subset_genepop("out.pure.gen.txt", subs=FST.Filter.Vec, keep=TRUE, path="out.temp.gen.txt")

# ld<- test_LD("out.pure.gen.txt", outputFile = "out.ld", settingsFile = "",
# 		dememorization = 10000, batches = 100,  iterations = 500, verbose = F)
#print(top_n_loci)
#print(FST.Filter.Vec[1:top_n_loci])
#get only the top X number of them
if (length(FST.Filter.Vec) > top_n_loci){
	FST.Filter.Vec <- FST.Filter.Vec[1:top_n_loci]
}
#print(FST.Filter.Vec)

#filter genepop data to only the top loci 
subset_genepop(GenePopData, subs=FST.Filter.Vec, keep=TRUE, spop=c("PureA", "PureB"), path="out.pure.topLoci.gen")
subset_genepop(GenePopData, subs=FST.Filter.Vec, keep=TRUE, path="out.all.topLoci.gen")

