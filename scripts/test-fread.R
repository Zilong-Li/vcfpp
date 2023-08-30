library(stringr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
vcffile <- args[1]

cmd=paste("bcftools","view","-H", vcffile)
system.time(vcf <- fread(cmd=cmd, header=FALSE, sep="\t", data.table=FALSE))
object.size(vcf)

calc_hets1 <- function(gt) {
  hets <- apply(gt, 2, function(a) {
    o <- sapply(str_split(a, fixed("|")), as.numeric)
    sum(colSums(o)==1)
  })
  hets
}

calc_hets2 <- function(gt) {
  out <- apply(gt, 2, str_split, fixed("|"))
  hets <- sapply(out, function(d) {
    o <- sapply(d, as.integer)
    sum(abs(o[1,] - o[2,]))
  })
  hets
}

gt <- vcf[,-c(1:9)]

system.time(hets<-calc_hets1(gt))
