#include <Rcpp.h>
#include <vcfpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List table.GT(const string & vcffile, const std::string & region, std::string samples="-")
{
    vcfpp::BcfReader vcf(vcffile, samples);
    vcfpp::BcfRecord var(vcf.header);
    int nsnps = vcf.getRegionIndex(region);
    CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps), info(nsnps);
    IntegerVector pos(nsnps);
    NumericVector qual(nsnps);
    vector<vector<int>> GT(nsnps);
    vector<int> gt;
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
List table.GL(const string & vcffile, const std::string & region, std::string samples="-")
{
    vcfpp::BcfReader vcf(vcffile, samples);
    vcfpp::BcfRecord var(vcf.header);
    int nsnps = vcf.getRegionIndex(region);
    CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps), info(nsnps);
    IntegerVector pos(nsnps);
    NumericVector qual(nsnps);
    vector<vector<int>> GL(nsnps);
    vector<int> gl;
    for(int i = 0; i < nsnps; i++)
    {
        vcf.getNextVariant(var);
        var.getFORMAT("GL",gl);
        GL[i] = gl;
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
                        Named("info") = info, Named("gl") = GL);
}

