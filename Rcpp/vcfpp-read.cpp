#include <Rcpp.h>
#include "vcfpp.h"

using namespace Rcpp;
using namespace std;

//' parse GT of a VCF file into tables in R
//'
//' @param vcffile path to the VCF file with index
//' @param region  region to extract 
//' @param samples samples to extract
//' @return A list of genotypes and other fixed fields in VCF
//' @export
// [[Rcpp::export]]
List tableGT(std::string vcffile, std::string region, std::string samples="-")
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
        pos(i) = var.POS();
        qual(i) = var.QUAL();
        chr(i) = var.CHROM();
        id(i) = var.ID();
        ref(i) = var.REF();
        alt(i) = var.ALT();
        filter(i) = var.FILTER();
        info(i) = var.INFO();
        var.getGenotypes(gt);
        GT[i] = gt;
    }
    return List::create(Named("chr") = chr, Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                        Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                        Named("info") = info, Named("gt") = GT);
}

//' parse GL of a VCF file into tables in R
//'
//' @param vcffile path to the VCF file with index
//' @param region  region to extract 
//' @param samples samples to extract
//' @return A list of genotypes and other fixed fields in VCF
//' @export
// [[Rcpp::export]]
List tableGL(std::string vcffile, std::string region, std::string samples = "-")
{
    vcfpp::BcfReader vcf(vcffile, samples);
    vcfpp::BcfRecord var(vcf.header);
    int nsnps = vcf.getRegionIndex(region);
    CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps), info(nsnps);
    IntegerVector pos(nsnps);
    NumericVector qual(nsnps);
    vector<vector<float>> GL(nsnps); 
    vector<float> gl;
    for(int i = 0; i < nsnps; i++)
    {
        vcf.getNextVariant(var);
        pos(i) = var.POS();
        qual(i) = var.QUAL();
        chr(i) = var.CHROM();
        id(i) = var.ID();
        ref(i) = var.REF();
        alt(i) = var.ALT();
        filter(i) = var.FILTER();
        info(i) = var.INFO();
        var.getFORMAT("GL",gl);
        GL[i] = gl;
    }
    return List::create(Named("chr") = chr, Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                        Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                        Named("info") = info, Named("gl") = GL);
}


//' parse PL of a VCF file into tables in R
//'
//' @param vcffile path to the VCF file with index
//' @param region  region to extract 
//' @param samples samples to extract
//' @return A list of genotypes and other fixed fields in VCF
//' @export
// [[Rcpp::export]]
List tablePL(std::string vcffile, std::string region, std::string samples = "-")
{
    vcfpp::BcfReader vcf(vcffile, samples);
    vcfpp::BcfRecord var(vcf.header);
    int nsnps = vcf.getRegionIndex(region);
    CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps), info(nsnps);
    IntegerVector pos(nsnps);
    NumericVector qual(nsnps);
    vector<vector<int>> PL(nsnps); 
    vector<int> pl;
    for(int i = 0; i < nsnps; i++)
    {
        vcf.getNextVariant(var);
        pos(i) = var.POS();
        qual(i) = var.QUAL();
        chr(i) = var.CHROM();
        id(i) = var.ID();
        ref(i) = var.REF();
        alt(i) = var.ALT();
        filter(i) = var.FILTER();
        info(i) = var.INFO();
        var.getFORMAT("PL",pl);
        PL[i] = pl;
    }
    return List::create(Named("chr") = chr, Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                        Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                        Named("info") = info, Named("pl") = PL);
}
