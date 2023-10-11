library(vcfR)

run <- 1
args <- commandArgs(trailingOnly = TRUE)
vcffile <- args[1]

print(system.time(vcf <- read.vcfR(vcffile)))

gt <- extract.gt(vcf[is.biallelic(vcf),], element = 'GT', as.numeric = TRUE)
print(system.time(hets <- colSums(gt==1, na.rm = TRUE)))

