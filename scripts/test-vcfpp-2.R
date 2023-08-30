library(Rcpp)
Sys.setenv("PKG_LIBS"="-I/home/rlk420/mambaforge/envs/R/include -lhts")
system.time(sourceCpp("test-vcfpp-r.cpp", verbose=TRUE, rebuild=TRUE))

args <- commandArgs(trailingOnly = TRUE)
vcffile <- args[1]

system.time(vcf <- readtable(vcffile, "chr21"))

gt <- do.call(rbind, vcf[["gt"]])
n <- ncol(gt)

hets <- colSums(abs(gt[,seq(1,n,2)] - gt[,seq(2,n, 2)]))
head(hets)

