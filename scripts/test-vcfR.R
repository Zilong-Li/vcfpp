library(Rcpp)
library(stringr)
library(vcfR)

run <- 1
args <- commandArgs(trailingOnly = TRUE)
vcffile <- args[1]
run <- as.integer(args[2])

system.time(vcf <- read.vcfR(vcffile))

if(run == 2) {
  gt <- extract.gt(vcf[is.biallelic(vcf),], element = 'GT')
  hets <- apply(gt, 2, function(a) {
    o <- sapply(str_split(a, fixed("|")), as.numeric)
    sum(colSums(o)==1)
  })
}
