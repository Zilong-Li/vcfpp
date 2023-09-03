// -*- compile-command: "g++ flipGT.cpp -std=c++11 -g -O3 -Wall -lhts -lz -lm -lbz2 -llzma -lcurl" -*-
#include "vcfpp.h"

using namespace std;
using namespace vcfpp;

int main(int argc, char* argv[])
{
    // ========= helper message and parameters parsing ============================
    std::vector<std::string> args(argv + 1, argv + argc);
    if(argc <= 1 || args[0] == "-h" || args[0] == "-help" || args[0] == "--help")
    {
        std::cout << "Author: Zilong-Li (zilong.dk@gmail.com)\n"
                  << "Description:\n"
                  << "     swap REF and ALT and flip corresponding genotypes for only biallelics \n\n"
                  << "Usage example:\n"
                  << "     " + (std::string)argv[0] + " -i in.bcf \n"
                  << "     " + (std::string)argv[0] + " -i in.bcf -o out.bcf -s ^S1,S2 -r chr1:1-1000 \n"
                  << "     bcftools view in.bcf | " + (std::string)argv[0] + " -i - -o out.bcf \n"
                  << "\nOptions:\n"
                  << "     -i    input vcf/bcf file\n"
                  << "     -o    ouput vcf/bcf file [stdout]\n"
                  << "     -s    list of samples to be included or excluded\n"
                  << "     -r    specific region to be included\n"
                  << std::endl;
        return 1;
    }
    std::string invcf, outvcf = "-", samples = "-", region = "";
    for (int i = 0; i < args.size(); i++)
    {
        if (args[i] == "-i")
            invcf = args[++i];
        if (args[i] == "-o")
            outvcf = args[++i];
        if (args[i] == "-s")
            samples = args[++i];
        if (args[i] == "-r")
            region = args[++i];
    }
    // ========= core calculation part ===========================================
    BcfReader vcf(invcf, region, samples);
    BcfRecord var(vcf.header);
    BcfWriter bw(outvcf, vcf.header);
    std::vector<char> gts;
    while (vcf.getNextVariant(var))
    {
        var.getGenotypes(gts);
        var.swap_REF_ALT();
        for (auto& g : gts)
        {
            g = !g;
        }
        var.setGenotypes(gts);
        bw.writeRecord(var);
    }

    return 0;
}
