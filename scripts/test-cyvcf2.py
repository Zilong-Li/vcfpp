import cyvcf2
import numpy as np

vcf = cyvcf2.VCF("1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
sample_counts = np.zeros(len(vcf.samples), dtype=float)

for var in vcf:
    sample_counts[(var.gt_types == vcf.HET)] += 1
print(zip(vcf.samples, sample_counts))
