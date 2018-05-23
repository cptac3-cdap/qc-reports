
source("parsemzid.R")

args <- commandArgs(TRUE)
mzid.file=args[1]

mzid = parsemzid(mzid.file)
summarize("PSMs",mzid$psms)
summarize("Alignments",mzid$alignments)
