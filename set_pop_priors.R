
# Store working directory as object.
my_path <- getwd()

###############################################################################

args = commandArgs(trailingOnly=TRUE)
#print(args)
my.file <- as.character(args[1])
z0 <- as.numeric(args[2])
z1 <- as.character(args[3])


fpath<-file("out.sims/out.pure.topLoci.gen_S2R1_NH.txt", "r")
opath<-file("test.txt", "w")
stuff<-""
while ( TRUE ) {
	line = readLines(con, n = 1)
	if ( length(line) == 0 ) {
		break
	}
	print(line)
}

close(con)
}
