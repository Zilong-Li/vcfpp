// -*- compile-command: "g++ vcf2beagle.cpp -std=c++11 -g -O3 -Wall -lhts -lz" -*-
#include "vcfpp.h"
#include <zlib.h>
#include <cmath>

using namespace std;
using namespace vcfpp;

int main(int argc, char * argv[])
{
    // ========= helper message and parameters parsing ============================
    std::vector<std::string> args(argv + 1, argv + argc);
    if(argc <= 1 || args[0] == "-h" || args[0] == "-help" || args[0] == "--help")
    {
        std::cout << "Author: Zilong-Li (zilong.dk@gmail.com)\n"
                  << "Description:\n"
                  << "     convert vcf with PL or GL tag to beagle file\n\n"
                  << "Usage example:\n"
                  << "     " + (std::string)argv[0] + " -i in.bcf -o beagle.gz -t PL\n"
                  << "     " + (std::string)argv[0] + " -i in.bcf -o beagle.gz -t GL -s ^S1,S2 -r chr1:1-1000 \n"
                  << "\nOptions:\n"
                  << "     -i    input vcf/bcf file\n"
                  << "     -o    ouput vcf/bcf file [stdout]\n"
                  << "     -s    list of samples to be included or excluded\n"
                  << "     -r    specific region to be included\n"
                  << "     -t    tag of genotype likelihoods in VCF. [PL or GL]\n"
                  << std::endl;
        return 1;
    }
    std::string invcf, outfile = "begale.gz", samples = "-", region = "", tag = "PL";
    for(size_t i = 0; i < args.size(); i++)
    {
        if(args[i] == "-i") invcf = args[++i];
        if(args[i] == "-o") outfile = args[++i];
        if(args[i] == "-s") samples = args[++i];
        if(args[i] == "-r") region = args[++i];
        if(args[i] == "-t") tag = args[++i];
    }
    // ========= core calculation part ===========================================
    BcfReader vcf(invcf, region, samples);
    BcfRecord var(vcf.header);
    int nsamples = vcf.nsamples;
    vector<int> pl;
    vector<float> gl;
    std::string hdr{"marker\tallele1\tallele2"}, line;
    for(auto s : vcf.header.getSamples()) hdr += "\t" + s + "\t" + s + "\t" + s;
    hdr += "\n";
    gzFile gzfp = gzopen(outfile.c_str(), "wb");
    gzwrite(gzfp, hdr.c_str(), hdr.size());
    vector<double> out;
    while(vcf.getNextVariant(var))
    {
        if(!var.isSNP()) continue;
        if(tag == "PL")
        {
            var.getFORMAT("PL", pl); // try get PL values
            out.resize(pl.size());
            for(size_t i = 0; i < pl.size(); i++) out[i] = std::pow(10, -pl[i] / 10);
        }
        else if(tag == "GL")
        {
            var.getFORMAT("GL", gl); // try get GL values
            out.resize(gl.size());
            for(size_t i = 0; i < gl.size(); i++) out[i] = std::pow(10, gl[i]);
        }
        else
        {
            throw runtime_error("please specify -t PL or GL\n");
        }
        // normalize it
        line = var.CHROM() + "_" + to_string(var.POS()) + "\t" + var.REF() + "\t" + var.ALT();
        for(int i = 0; i < nsamples; i++)
        {
            double s = out[i] + out[i + 1] + out[i + 2];
            out[i] /= s;
            out[i + 1] /= s;
            out[i + 2] /= s;
            line += "\t" + to_string(out[i]) + "\t" + to_string(out[i + 1]) + "\t" + to_string(out[i + 2]);
        }
        line += "\n";
        gzwrite(gzfp, line.c_str(), line.size());
    }
    gzclose(gzfp);

    return 0;
}
