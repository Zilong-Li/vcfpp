library(Rcpp)
Sys.setenv("PKG_LIBS"="-I~/mambaforge/envs/R/include -L~/mambaforge/envs/R/include -lhts")
sourceCpp("vcfpp-read.cpp", verbose=TRUE, rebuild=TRUE)

vcfgt <- tableGT("vcf.gz", "chr20")
vcfgl <- tableGL("vcf.gz", "chr20")
