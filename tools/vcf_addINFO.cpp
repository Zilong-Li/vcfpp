// -*- compile-command: "g++ vcf_addINFO.cpp -std=c++11 -g -O3 -Wall -I.. -lhts" -*-
#include "vcfpp.h"
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
                  << "     calculate INFO score per site given GP from input vcf file\n\n"
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
    for(size_t i = 0; i < args.size(); i++)
    {
        if(args[i] == "-i") invcf = args[++i];
        if(args[i] == "-o") outvcf = args[++i];
        if(args[i] == "-s") samples = args[++i];
        if(args[i] == "-r") region = args[++i];
    }
    // ========= core calculation part ===========================================
    BcfReader vcf(invcf, region, samples);
    BcfWriter bw(outvcf, vcf.header);
    bw.header.addSample(samples);
    bw.header.addINFO("INFO", "1", "Float", "INFO score given genotype probability");
    bw.header.addINFO("EAF", "1", "Float", "Estimated Allele Frequency");
    int N = vcf.nsamples, i{0};
    vector<float> gps;
    float info, eaf, eij, fij, a0, a1, thetaHat;
    BcfRecord var(bw.header);
    while(vcf.getNextVariant(var))
    {
        var.getFORMAT("GP", gps); // try get GP values
        for(eij = 0.0, fij = 0.0, i = 0; i < N; i++)
        {
            a0 = gps[i * 3 + 1] + gps[i * 3 + 2] * 2;
            a1 = gps[i * 3 + 1] + gps[i * 3 + 2] * 4;
            eij += a0;
            fij += a1 - a0 * a0;
        }
        eaf = eij / 2 / N;
        info = 1.0 - fij / (eij * (1 - eaf));
        thetaHat = std::round(1e2 * eaf) / 1e2;
        if(thetaHat == 0 || thetaHat == 1)
            info = 1.0;
        else if(info < 0)
            info = 0.0;
        else
            info = std::round(1e3 * info) / 1e3;
        eaf = std::round(1e6 * eaf) / 1e6;
        var.setINFO("INFO", info);
        var.setINFO("EAF", eaf);
        bw.writeRecord(var);
    }

    return 0;
}
