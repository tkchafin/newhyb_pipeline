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
#print(args)
sims <- as.numeric(args[1])
reps <- as.numeric(args[2])
sim.type <- as.character(args[3])
ref.file <- "out.pure.gen"

#print(sim.type)
if (sim.type == "GT"){
	freqbasedsim_GTFreq(ref.file, NumSims = sims, NumReps=reps, pop.groups=c("PureA", "PureB"))
}else if (sim.type == "AS"){
	freqbasedsim_AlleleSample(ref.file, NumSims = sims, NumReps=reps, pop.groups=c("PureA", "PureB"))
}else{
	print("ERROR: Unsupported sim.type")
}
