import sys
import cyvcf2
import numpy as np

vcffile = sys.argv[1]
# vcffile = "1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
vcf = cyvcf2.VCF(vcffile)
sample_counts = np.zeros(len(vcf.samples), dtype=float)

for var in vcf:
    sample_counts[(var.gt_types == vcf.HET)] += 1
print(sample_counts)
