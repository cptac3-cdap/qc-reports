
source("parsemzid.R")

args <- commandArgs(TRUE)
mzid.file=args[1]

results = parsemzid(mzid.file)
summarize("PSMs",results$psms)
summarize("Alignments",results$alignments)
