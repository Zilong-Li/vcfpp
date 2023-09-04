library(stringr)
library(vcfR)

run <- 1
args <- commandArgs(trailingOnly = TRUE)
vcffile <- args[1]
run <- as.integer(args[2])

system.time(vcf <- read.vcfR(vcffile))

if(run == 2) {
  gt <- extract.gt(vcf[is.biallelic(vcf),], element = 'GT', as.numeric = TRUE)
  system.time(hets <- apply(gt, 2, function(g) sum(g==1)))
}
