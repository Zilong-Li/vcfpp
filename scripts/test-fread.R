library(data.table)
args <- commandArgs(trailingOnly = TRUE)
vcffile <- args[1]
cmd=paste("bcftools","view","-H", vcffile)
system.time(vcf <- fread(cmd=cmd, header=FALSE, sep="\t", data.table=FALSE))
