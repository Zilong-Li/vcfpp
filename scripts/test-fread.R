library(Rcpp)
library(data.table)

args <- commandArgs(trailingOnly = FALSE)
vcffile <- args[1]
system.time(vcf <- fread(vcffile, header=TRUE, sep="\t", data.table=FALSE))
