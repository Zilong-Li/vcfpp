library(Rcpp)
Sys.setenv("PKG_LIBS"="-I/home/rlk420/mambaforge/envs/R/include -lhts")
sourceCpp("test-vcfpp-r.cpp", verbose=TRUE, rebuild=TRUE)

run <- 1
args <- commandArgs(trailingOnly = TRUE)
vcffile <- args[1]
run <- as.integer(args[2])

if (run==1) {
  system.time(gt <- genotypes(vcffile))
  system.time(info <- varinfos(vcffile, length(gt)))
}

if(run == 2) {
  system.time(hets <- hetrate(vcffile))
}
