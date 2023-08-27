library(Rcpp)
library(vcfR)

args <- commandArgs(trailingOnly = TRUE)
vcffile <- args[1]
system.time(vcf <- read.vcfR(vcffile))
