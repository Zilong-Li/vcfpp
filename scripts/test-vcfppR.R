## devtools::install_github("Zilong-Li/vcfppR")

library(vcfppR)

args <- commandArgs(trailingOnly = TRUE)
vcffile <- args[1]
run <- as.integer(args[2])

if(run == 1) {
  print(paste("run", run, "-- wrapper of C++ code"))
  print(system.time(res1 <- heterozygosity(vcffile)))
  q(save="no")
}

if(run == 2) {
  print(paste("run", run, "-- with vcftable function in vcfppR"))
  print(system.time(vcf <- vcftable(vcffile, vartype = "snps")))
  print(system.time(res2 <- colSums(vcf[["gt"]]==1, na.rm = TRUE)))
  q(save="no")
}

if(run == 3) {
  print(paste("run", run, "-- with vcfreader API in vcfppR"))
  vcf <- vcfreader$new(vcffile)
  res <- rep(0L, vcf$nsamples())
  while(vcf$variant()) {
    if(vcf$isSNP()) {
      gt <- vcf$genotypes(TRUE) == 1
      gt[is.na(gt)] <- FALSE
      res <- res + gt
    }
  }
  q(save="no")
}
