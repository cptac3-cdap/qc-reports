
library(rmarkdown)

# figure out the directory containing this script
# when run by Rscript

cmdArgs <- commandArgs(trailingOnly=FALSE)
needle <- "--file="
match <- grep(needle, cmdArgs)
this.dir <- dirname(sub(needle, "", cmdArgs[match]))

args <- commandArgs(TRUE)
QCmetrics_file <- args[1]
Base_filename <- sub(".qcmetrics.tsv","",QCmetrics_file)
Summary_file <- paste0(Base_filename,".summary.tsv",sep="")
Peptide_file <- paste0(Base_filename,".peptides.tsv",sep="")
Mayu_file <- paste0(Base_filename,".mayu.tsv",sep="")
expdes_file <- paste0(Base_filename,".sample.csv",sep="")

types  = c("",".phosphosite",".acetylsite",".ubiquitylsite",".glycosite",".phosphopeptide",".acetylpeptide","ubiquitylpeptide",".glycopeptide",".peptide")
rawopt = c("","-raw")
labels = c("tmt10","tmt11","tmt6","tmt16","tmt18","itraq","spectral_counts","precursor_area")

for (type in types) {
  for (raw in rawopt) {
    for (label in labels) {
      Quant_file <- paste0(Base_filename,type,".",label,raw,".tsv",sep="")
      # print(Quant_file)
      if (file.exists(Quant_file)) {
        break;
      }
    }
    if (file.exists(Quant_file)) {
      break;
    }
  } 
  if (file.exists(Quant_file)) {
    break;
  }
}
stopifnot(file.exists(Quant_file))

report_file <- args[2]

extraparams = list(dochecks=FALSE,noexpdesign=FALSE)
if (length(args) > 2) {
  for (i in 3:length(args)) {
      extraparams[args[i]] = TRUE
  } 
}

work.dir <- dirname(report_file)
for (rmd in Sys.glob(file.path(this.dir,"*.Rmd"))) {
    dummy <- file.copy(rmd,work.dir)
}

if (endsWith(report_file,'.pdf') == TRUE) {
    outformat <- 'pdf_document'
} else if (endsWith(report_file, '.html') == TRUE) {
    outformat <- 'html_document'
}

renderparams <- list(qcmetricsfile = QCmetrics_file,
                     summaryfile = Summary_file,
	             peptidefile = Peptide_file,
                     quantfile = Quant_file,
                     mayufile = Mayu_file,
		     expdesfile = expdes_file,
	             format = outformat,
		     dochecks = extraparams$dochecks,
                     noexpdesign = extraparams$noexpdesign) 

rmarkdown::render(file.path(work.dir,'qcreport.Rmd'),
                  output_format = outformat, 
                  params = renderparams,
                  output_file = report_file,
                  quiet = TRUE)
