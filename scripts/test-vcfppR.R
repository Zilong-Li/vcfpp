## devtools::install_github("Zilong-Li/vcfppR")

library(vcfppR)

run <- 1
args <- commandArgs(trailingOnly = TRUE)
vcffile <- args[1]
run <- as.integer(args[2])

system.time(vcf <- tableGT(vcffile, "chr21"))

if(run == 2) {
  res <- sapply(vcf[["gt"]], function(a) {
    n=length(a)
    abs(a[seq(1,n,2)]-a[seq(2,n,2)])
  })
  hets<-rowSums(res)
}
