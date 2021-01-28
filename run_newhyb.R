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
my.NH <- as.character(args[1])
my.Files <- as.character(args[2])
burn <- as.numeric(args[3])
sweeps <- as.numeric(args[4])


get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
  os <- sysinf['sysname']
  if (os == 'Darwin')
    os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

#print(sim.type)
os<-get_os()
if (os == "osx"){
	parallelnh_OSX(my.Files, where.NH = my.NH, burnin=burn, sweeps=sweeps)
}else if (sim.type == "linux"){
	parallelnh_LINUX(my.Files, where.NH = my.NH, burnin=burn, sweeps=sweeps)
}else{
	print("ERROR: Unsupported OS")
}
