#include "../vcfpp.h"

using namespace std;
using namespace vcfpp;

int main(int argc, char ** argv)
{
    std::string invcf = "-", outvcf = "-", samples = "-", region = "";
    BcfReader vcf(invcf, region, samples);
    BcfRecord var(vcf.header);
    vector<char> gts;
    int nsamples = vcf.nsamples;
    float sum{0};
    while(vcf.getNextVariant(var))
    {
        var.getGenotypes(gts);
        for(int i = 0; i < nsamples; i++) sum += gts[i * 2] + gts[i * 2 + 1];
        cout << sum << endl;
    }
}
