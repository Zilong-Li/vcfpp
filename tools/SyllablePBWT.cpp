// -*- compile-command: "g++ SyllablePBWT.cpp -o SyllablePBWT -std=c++11 -O3 -Wall -lhts -lz -lm -lbz2 -llzma -lcurl" -*-
#include "SyllablePBWT.h"

int main(int argc, char * argv[])
{

    // ========= helper message and parameters parsing ============================
    const std::vector<std::string> args(argv + 1, argv + argc);
    if(argc <= 1 || args[0] == "-h" || args[0] == "-help")
    {
        std::cout << "Author: Zilong-Li (zilong.dk@gmail.com)\n"
                  << "Usage example:\n"
                  << "     " + (std::string)argv[0]
                         + " -panel panel.vcf.gz -save index.bin -region chr22 -samples '^HG00103'\n"
                  << "     " + (std::string)argv[0]
                         + " -query panel.vcf.gz -load index.bin -region chr22 -samples HG00103 -L 1000\n"
                  << "\nOptions:\n"
                  << "     -panel      vcf/bcf file of reference panel\n"
                  << "     -query      vcf/bcf file to query\n"
                  << "     -save       save mspbwt indicies as binary file\n"
                  << "     -load       load mspbwt indicies from binary file\n"
                  << "     -region     chromosome to be included, bcftools-like format\n"
                  << "     -samples    sample id to be included, bcftools-like format\n"
                  << "     -L          min number of matching sites\n"
                  << std::endl;
        return 1;
    }
    std::string vcfpanel, vcfquery, binfile, outfile,samples = "-", region = "";
    int L{0};
    for(size_t i = 0; i < args.size(); i++)
    {
        if(args[i] == "-panel") vcfpanel = args[++i];
        if(args[i] == "-query") vcfquery = args[++i];
        if(args[i] == "-out") outfile = args[++i];
        if(args[i] == "-save") binfile = args[++i];
        if(args[i] == "-load") binfile = args[++i];
        if(args[i] == "-region") region = args[++i];
        if(args[i] == "-samples") samples = args[++i];
        if(args[i] == "-L") L = stoi(args[++i]);
    }

    // ========= core calculation part ===========================================

    SyllablePBWT<unsigned long long> syp;
    if(!vcfpanel.empty())
    {
        syp.build(vcfpanel, region, samples);
        if(!binfile.empty()) syp.save(binfile);
    }
    else
    {
        if(!binfile.empty())
            syp.load(binfile);
        else
            throw std::invalid_argument("binary file with syllable pbwt indicies does not exist!\n");
        if(!vcfquery.empty())
            syp.query(vcfquery, region, samples, L);
        else
            throw std::invalid_argument("please provide -query file!\n");
    }

    return 0;
}
