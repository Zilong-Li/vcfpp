#include "../vcfpp.h"
#include "catch.hh"

using namespace vcfpp;
using namespace std;

TEST_CASE("Calculate the heterozygosity rate", "[bcf-reader]")
{
    BcfReader vcf("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/"
                  "working/20220422_3202_phased_SNV_INDEL_SV/"
                  "1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
                  "-", "chr22:19000000-20000000");
    BcfRecord var(vcf.header); // construct a variant record
    vector<char> gt; // genotype can be of bool, char or int type
    vector<int> hetsum(vcf.nsamples);
    while(vcf.getNextVariant(var))
    {
        if(!var.isSNP()) continue; // skip other type of variants
        var.getGenotypes(gt);
        for(int i = 0; i < gt.size() / var.nploidy; i++)
            hetsum[i] += abs(gt[var.nploidy * i + 0] - gt[var.nploidy * i + 1]);
    }
    for (auto i : hetsum) { cout << i << endl; }
}
