#include <Rcpp.h>
#include <vcfpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List genotypes(const string & vcffile)
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
