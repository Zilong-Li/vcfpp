library(Rcpp)
Sys.setenv("PKG_LIBS"="-I~/mambaforge/envs/R/include -L~/mambaforge/envs/R/include -lhts")
sourceCpp("Rcpp_example.cpp", verbose=TRUE, rebuild=TRUE)

## vcffile <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
hets <- heterozygosity(vcffile)


