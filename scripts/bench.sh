#!/bin/bash

# benchmarking performance of vcfpp-r against vcfR and fread

# use gnu time to record the RAM and TIME
gtime="/usr/bin/time"

# download vcf file from 1000 genome project

wget  -N -r --no-parent https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz

# big vcf file with very long INFO field
vcffile="1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

$gtime -vvv Rscript test-vcfR.R $vcffile &> test-vcfR.llog1 &
$gtime -vvv Rscript test-vcfpp.R $vcffile &> test-vcfpp.llog1 &
$gtime -vvv Rscript test-fread.R $vcffile &> test-fread.llog1 &
wait

# small vcf file with no INFO field
bcftools annotate -x INFO --threads 4 -Oz -o $vcffile2 $vcffile && bcftools index $vcffile2
vcffile2="small."$vcffile

$gtime -vvv Rscript test-vcfR.R $vcffile2 &> test-vcfR.llog2 &
$gtime -vvv Rscript test-vcfpp.R $vcffile2 &> test-vcfpp.llog2 &
$gtime -vvv Rscript test-fread.R $vcffile2 &> test-fread.llog2 & 
wait
