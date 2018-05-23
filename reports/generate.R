
library(rmarkdown)

# figure out the directory containing this script
# when run by Rscript

cmdArgs <- commandArgs(trailingOnly=FALSE)
needle <- "--file="
match <- grep(needle, cmdArgs)
this.dir <- dirname(sub(needle, "", cmdArgs[match]))

args <- commandArgs(TRUE)
QCmetrics_file <- args[1]
report_file <- args[2]

if (endsWith(report_file,'.pdf')) {
    outformat <- 'pdf_document'
} else if (endsWith(report_file, '.html')) {
    outformat <- 'html_document'
}
rmarkdown::render(file.path(this.dir,'qcreport.Rmd'),
                  output_format = outformat,
                  params = list(qcmetricsfile = QCmetrics_file,format = outformat), 
                  output_file = report_file)
