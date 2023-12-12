#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../vcfpp.h"

namespace py = pybind11;

using namespace vcfpp;
using namespace std;

vector<int> heterozygosity(std::string vcffile, std::string region = "", std::string samples = "")
{
    BcfReader vcf(vcffile, region, samples);
    BcfRecord var(vcf.header); // construct a variant record
    vector<int> gt; // genotype can be bool, char or int type
    vector<int> hetsum(vcf.nsamples, 0);
    while(vcf.getNextVariant(var))
    {
        var.getGenotypes(gt);
        // analyze SNP variant with no genotype missingness
        if(!var.isSNP()) continue; // analyze SNPs only
        assert(var.ploidy() == 2); // make sure it is diploidy
        for(size_t i = 0; i < gt.size() / 2; i++) hetsum[i] += abs(gt[2 * i + 0] - gt[2 * i + 1]) == 1;
    }
    return hetsum;
}

PYBIND11_MODULE(py_example, m)
{
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("heterozygosity", &heterozygosity, "A function that calculates the heterozygosity for each sample in the VCF");
}
