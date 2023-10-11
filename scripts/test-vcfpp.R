library(Rcpp)
Sys.setenv("PKG_LIBS"="-I/home/rlk420/mambaforge/envs/R/include -lhts")
sourceCpp("test-vcfpp-r.cpp", verbose=TRUE, rebuild=TRUE)

args <- commandArgs(trailingOnly = TRUE)
vcffile <- args[1]


system.time(hets <- heterozygosity(vcffile))
