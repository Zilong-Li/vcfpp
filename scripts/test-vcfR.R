library(stringr)
library(vcfR)

run <- 1
args <- commandArgs(trailingOnly = TRUE)
vcffile <- args[1]
run <- as.integer(args[2])

system.time(vcf <- read.vcfR(vcffile))

calc_hets1 <- function(gt) {
  hets <- apply(gt, 2, function(a) {
    o <- sapply(str_split(a, fixed("|")), as.numeric)
    sum(colSums(o)==1)
  })
  hets
}

if(run == 2) {
  gt <- extract.gt(vcf[is.biallelic(vcf),], element = 'GT')
  system.time(hets<-calc_hets1(gt))
}
