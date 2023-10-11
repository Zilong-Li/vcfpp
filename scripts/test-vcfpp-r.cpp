#include <Rcpp.h>
#include <vcfpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
IntegerVector heterozygosity(std::string vcffile, std::string region = "", std::string samples = "-")
{
    vcfpp::BcfReader vcf(vcffile, region, samples);
    vcfpp::BcfRecord var(vcf.header);
    vector<int> gt;
    vector<int> hetsum(vcf.nsamples, 0);  // store the het counts
    while (vcf.getNextVariant(var)) {
        var.getGenotypes(gt);
        if (!var.isSNP()) continue; // analyze SNPs only
        assert(var.ploidy() == 2);  // make sure it is diploidy
        for (int i = 0; i < gt.size() / 2; i++)
            hetsum[i] += abs(gt[2 * i + 0] - gt[2 * i + 1]) == 1;
    }
    return wrap(hetsum);
}
