#include <Rcpp.h>
#include <vcfpp.h>
using namespace Rcpp;
using namespace std;

List genotypes_bool(const string & vcffile)
{
    vcfpp::BcfReader vcf(vcffile);
    vcfpp::BcfRecord var(vcf.header);
    vector<vector<bool>> G;
    vector<bool> gt;
    while(vcf.getNextVariant(var))
    {
        var.getGenotypes(gt);
        G.push_back(gt);
    }
    return wrap(G);
}

List genotypes_int(const string & vcffile)
{
    vcfpp::BcfReader vcf(vcffile);
    vcfpp::BcfRecord var(vcf.header);
    vector<vector<int>> G;
    vector<int> gt;
    while(vcf.getNextVariant(var))
    {
        var.getGenotypes(gt);
        G.push_back(gt);
    }
    return wrap(G);
}

// [[Rcpp::export]]
List genotypes(const string & vcffile, bool use_bool = true)
{
    if(use_bool)
        return genotypes_bool(vcffile);
    else
        return genotypes_int(vcffile);
}

// [[Rcpp::export]]
List varinfos(const string & vcffile, int nsnps)
{
    vcfpp::BcfReader vcf(vcffile);
    vcfpp::BcfRecord var(vcf.header);
    CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps);
    IntegerVector pos(nsnps);
    NumericVector qual(nsnps);
    int i = 0;
    while(vcf.getNextVariant(var))
    {
        chr(i) = var.CHROM();
        pos(i) = var.POS();
        id(i) = var.ID();
        ref(i) = var.REF();
        alt(i) = var.ALT();
        qual(i) = var.QUAL();
        filter(i) = var.FILTER();
        i++;
    }
    return List::create(Named("chr") = chr, Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                        Named("alt") = alt, Named("qual") = qual, Named("filter") = filter);
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
