library(Rcpp)
library(vcfR)

args <- commandArgs(trailingOnly = FALSE)
vcffile <- args[1]
system.time(vcf <- read.vcfR(vcffile))
