#include "../vcfpp.h"

using namespace std;
using namespace vcfpp;

int main(int argc, char ** argv)
{
    std::string invcf = "-", outvcf = "-", samples = "-", region = "";
    BcfReader vcf(invcf, region, samples);
    const int nsamples = vcf.nsamples;
    // BcfWriter bw(outvcf, vcf.header);
    BcfWriter bw(outvcf, vcf.header);
    bw.header.setSamples(samples);
    bw.header.addFORMAT("DS", "1", "Float", "Diploid Genotype Dosage"); // add DS tag into the header
    BcfRecord var(vcf.header);
    vector<char> gts;
    vector<float> ds(nsamples);
    while(vcf.getNextVariant(var))
    {
        var.getGenotypes(gts);
        for(int i = 0; i < nsamples; i++) ds[i] = gts[i * 2] + gts[i * 2 + 1];
        var.setFORMAT("DS", ds);
        bw.writeRecord(var);
    }
}
