## devtools::install_github("Zilong-Li/vcfppR")

library(vcfppR)

args <- commandArgs(trailingOnly = TRUE)
vcffile <- args[1]
run <- as.integer(args[2])

if(run == 1) {
  print(paste("run", run))
  print(system.time(res1 <- heterozygosity(vcffile)))
  q(save="no")
}

if(run == 2) {
  print(paste("run", run))
  print(system.time(vcf <- vcftable(vcffile, vartype = "snps")))
  print(system.time(res2 <- colSums(vcf[["gt"]]==1, na.rm = TRUE)))
  q(save="no")
}

