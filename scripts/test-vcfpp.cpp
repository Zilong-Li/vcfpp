// -*- compile-command: "x86_64-conda-linux-gnu-c++ test-vcfpp.cpp -o test-vcfpp -std=c++11 -O3 -Wall -I../ -I/home/rlk420/mambaforge/envs/R/include -lhts" -*-

#include "vcfpp.h"
using namespace std;
using namespace vcfpp;

int main(int argc, char * argv[])
{
    // ========= helper message and parameters parsing ============================
    std::vector<std::string> args(argv + 1, argv + argc);
    if(argc <= 1 || args[0] == "-h" || args[0] == "-help")
    {
        std::cout << "Author: Zilong-Li (zilong.dk@gmail.com)\n"
                  << "Description:\n"
                  << "     calculate heterzygous genotypes per sample \n\n"
                  << "Usage example:\n"
                  << "     " + (std::string)argv[0] + " -i in.bcf \n"
                  << "     " + (std::string)argv[0] + " -i in.bcf -o out.bcf -s ^S1,S2 -r chr1:1-1000 \n"
                  << "     bcftools view in.bcf | " + (std::string)argv[0] + " -i - -o out.bcf \n"
                  << "\nOptions:\n"
                  << "     -i    input vcf/bcf file\n"
                  << "     -s    list of samples to be included or excluded\n"
                  << "     -r    specific region to be included\n"
                  << std::endl;
        return 1;
    }
    string vcffile, samples = "-", region = "";
    for (int i = 0; i < args.size(); i++)
    {
        if (args[i] == "-i")
            vcffile = args[++i];
        if (args[i] == "-s")
            samples = args[++i];
        if (args[i] == "-r")
            region = args[++i];
    }
    // ========= core calculation part ===========================================
    BcfReader vcf(vcffile, region, samples);
    BcfRecord var(vcf.header); // construct a variant record
    vector<int> gt; // genotype can be bool, char or int type
    vector<int> hetsum(vcf.nsamples, 0);
    while(vcf.getNextVariant(var))
    {
        var.getGenotypes(gt);
        // analyze SNP variant with no genotype missingness
        if (!var.isSNP()) continue; // analyze SNPs only
        assert(var.ploidy() == 2);  // make sure it is diploidy
        for (int i = 0; i < gt.size() / 2; i++)
            hetsum[i] += abs(gt[2 * i + 0] - gt[2 * i + 1]) == 1;
    }
    // for(auto i : hetsum) cout << i << endl;
    return 0;
}
