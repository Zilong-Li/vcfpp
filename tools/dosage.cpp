// -*- compile-command: "g++ dosage.cpp -std=c++11 -g -O3 -Wall -lhts -lz -lm -lbz2 -llzma -lcurl" -*-
#include "../vcfpp.h"

using namespace std;
using namespace vcfpp;

int main(int argc, char* argv[])
{

    std::vector<std::string> args(argv + 1, argv + argc);
    if (argc <= 1 || args[0] == "-h" || args[0] == "-help")
    {
        std::cout << "Description:\n"
                  << "      create DS tag for diploid samples given GP tag from input vcf file\n\n"
                  << "Usage example:\n"
                  << "     " + (std::string)argv[0] + " -i in.bcf -o out.bcf \n"
                  << "     " + (std::string)argv[0] + " -i in.bcf -o out.bcf -s ^S1,S2 -r chr1:1-1000 \n"
                  << "      bcftools view in.bcf | " + (std::string)argv[0] + " -i - -o out.bcf \n"
                  << "\nOptions:\n"
                  << "     -i    vcf/bcf file of reference panel\n"
                  << "     -o    vcf/bcf file of query sequence\n"
                  << "     -s    list of samples to be included or excluded\n"
                  << "     -r    specific region be included\n"
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
    bw.header.addFORMAT("DS", "1", "Float", "Diploid Genotype Dosage"); // add DS tag into the header
    bw.writeHeader();                                                   // output header
    int nsamples = vcf.nsamples;
    std::vector<float> gp, ds(nsamples);
    while (vcf.getNextVariant(var))
    {
        var.getFORMAT("GP", gp); // get GP values
        for (int i = 0; i < nsamples; i++)
        {
            ds[i] = gp[i * 3 + 1] + gp[i * 3 + 2] * 2;
        }
        var.setFORMAT("DS", ds);
        bw.writeRecord(var);
    }
    // bw.close();
    return 0;
}
