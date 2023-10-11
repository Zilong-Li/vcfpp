#!/bin/bash

# benchmarking performance of vcfppR against cyvcf2, vcfR and fread

# use gnu time to record the RAM and TIME
gtime="/usr/bin/time"

# download vcf file from 1000 genome project
vcffile="1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
if [ ! -f $vcffile ];then
   wget -N -r --no-parent --no-directories https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
   wget -N -r --no-parent --no-directories https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi
fi

   
# first compile a C++ script
x86_64-conda-linux-gnu-c++ test-vcfpp.cpp -o test-vcfpp -std=c++11 -O3 -Wall -I/home/rlk420/mambaforge/envs/R/include -lhts


# to be fair, make sure cyvcf2 and vcfpp are compiled against same version of htslib.
# do the benchmarking in same conda env

$gtime -vvv Rscript test-fread.R $vcffile &> test-fread.llog.1 &
$gtime -vvv Rscript test-vcfR.R $vcffile 2 &> test-vcfR.llog.2 &
$gtime -vvv Rscript test-vcfppR.R $vcffile 1 &> test-vcfppR.llog.1 &
$gtime -vvv Rscript test-vcfppR.R $vcffile 2 &> test-vcfppR.llog.2 &
$gtime -vvv python test-cyvcf2.py $vcffile &> test-cyvcf2.llog &
$gtime -vvv ./test-vcfpp -i $vcffile &> test-vcfpp.llog &
wait

echo "jobs done. god day"
exit

