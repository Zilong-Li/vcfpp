library(vcfppR)

## Listing 4
vcffile <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.recalibrated_variants.vcf.gz"
vcf <- vcftable(vcffile, region="chr21:1-10000000", samples="NA12878,HG00118,HG00119", format="DP", vartype="snps", pass = TRUE, info = FALSE)
boxplot(vcf$dp, names=vcf$samples, ylab="Read Depth (DP)")
