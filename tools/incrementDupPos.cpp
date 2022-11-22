// -*- compile-command: "g++ incrementDupPos.cpp -std=c++11 -g -O3 -Wall -lhts -lz -lm -lbz2 -llzma -lcurl" -*-
#include "../vcfpp.h"

using namespace std;
using namespace vcfpp;


int main(int argc, char* argv[])
{
    std::vector<std::string> args(argv + 1, argv + argc);
    if (argc <= 1 || args[0] == "-h" || args[0] == "-help")
    {
        std::cout << "Author: Zilong-Li (zilong.dk@gmail.com)\n"
                  << "Description:\n"
                  << "     increment POS for duplicated biallelics sites from sorted vcf\n\n"
                  << "Usage example:\n"
                  << "     " + (std::string)argv[0] + " -i in.bcf -o out.bcf \n"
                  << "     " + (std::string)argv[0] + " -i in.bcf -o out.bcf -s ^S1,S2 -r chr1:1-1000 \n"
                  << "     bcftools view in.bcf | " + (std::string)argv[0] + " -i - -o out.bcf \n"
                  << "\nOptions:\n"
                  << "     -i    input vcf/bcf file\n"
                  << "     -o    ouput vcf/bcf file\n"
                  << "     -s    list of samples to be included or excluded\n"
                  << "     -r    specific region to be included\n"
                  << std::endl;
        return 1;
    }
    std::string invcf, outvcf, samples = "-", region = "";
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
    BcfReader vcf(invcf, samples, region);
    BcfRecord var(vcf.header);
    BcfWriter bw(outvcf);
    bw.initalHeader(vcf.header);
    bw.writeHeader();
    int64_t pos{-1}, count{0};
    if (vcf.getNextVariant(var))
    {
        pos = var.POS();
        bw.writeRecord(var);
    }
    while (vcf.getNextVariant(var))
    {
        if (var.POS() == pos)
        {
            count++;
            var.setPOS(pos + count);
        }
        else
        {
            pos = var.POS();
            count = 0;
        }
        bw.writeRecord(var);
    }

    return 0;
}
