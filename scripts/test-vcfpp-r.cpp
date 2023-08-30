#include <Rcpp.h>
#include <vcfpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
int getRegionIndex(const string & vcffile, const std::string & region)
{
    vcfpp::BcfReader vcf(vcffile);
    return vcf.getRegionIndex(region);
}

// [[Rcpp::export]]
List readtable(const string & vcffile, const std::string & region)
{
    vcfpp::BcfReader vcf(vcffile);
    vcfpp::BcfRecord var(vcf.header);
    int nsnps = vcf.getRegionIndex(region);
    CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps), info(nsnps);
    IntegerVector pos(nsnps);
    NumericVector qual(nsnps);
    vector<vector<bool>> GT(nsnps);
    vector<bool> gt;
    for(int i = 0; i < nsnps; i++)
    {
        vcf.getNextVariant(var);
        var.getGenotypes(gt);
        GT[i] = gt;
        pos(i) = var.POS();
        qual(i) = var.QUAL();
        chr(i) = var.CHROM();
        id(i) = var.ID();
        ref(i) = var.REF();
        alt(i) = var.ALT();
        filter(i) = var.FILTER();
        info(i) = var.INFO();
    }
    return List::create(Named("chr") = chr, Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                        Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                        Named("info") = info, Named("gt") = GT);
}

// [[Rcpp::export]]
IntegerVector hetrate(const string & vcffile)
{
    vcfpp::BcfReader vcf(vcffile);
    vcfpp::BcfRecord var(vcf.header);
    vector<char> gt;
    vector<int> hetsum(vcf.nsamples, 0); // store the het counts
    while(vcf.getNextVariant(var))
    {
        var.getGenotypes(gt);
        // analyze SNP variant with no genotype missingness
        if(!var.isSNP() || !var.isNoneMissing()) continue;
        assert(var.ploidy() == 2); // make sure it is diploidy
        for(int i = 0; i < gt.size() / 2; i++) hetsum[i] += abs(gt[2 * i + 0] - gt[2 * i + 1]);
    }
    return wrap(hetsum);
}
