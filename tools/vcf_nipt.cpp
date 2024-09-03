// -*- compile-command: "g++ vcf_nipt.cpp -std=c++11 -g -O3 -Wall -I.. -lhts" -*-
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
                  << "     process the VCF outputted by QUILT2 NIPT\n\n"
                  << "Usage example:\n"
                  << "     " + (std::string)argv[0] + " -i in.bcf \n"
                  << "     " + (std::string)argv[0] + " -i in.bcf -o out.bcf -s ^S1,S2 -r chr1:1-1000 \n"
                  << "     bcftools view in.bcf | " + (std::string)argv[0] + " -i - -o out.bcf \n"
                  << "\nOptions:\n"
                  << "     -i    input vcf/bcf file\n"
                  << "     -o    ouput prefix. \n"
                  << "     -s    list of samples to be included or excluded\n"
                  << "     -r    specific region to be included\n"
                  << std::endl;
        return 1;
    }
    std::string invcf, outvcf = "nipt", samples = "-", region = "";
    for(size_t i = 0; i < args.size(); i++)
    {
        if(args[i] == "-i") invcf = args[++i];
        if(args[i] == "-o") outvcf = args[++i];
        if(args[i] == "-s") samples = args[++i];
        if(args[i] == "-r") region = args[++i];
    }
    // ========= core calculation part ===========================================
    BcfReader vcf(invcf, region, samples);
    std::string outmat = outvcf + ".mat.vcf.gz";
    std::string outfet = outvcf + ".fet.vcf.gz";
    BcfWriter mw(outmat, "VCF4.2"), fw(outfet, "VCF4.2");
    {
        mw.header.addFORMAT("GT", "1", "String", "Phased genotypes");
        mw.header.addFORMAT("GP", "3", "Float", "Genotype probabilities");
        mw.header.addFORMAT("DS", "1", "Float", "Genotype dosages");
        mw.header.addFILTER("PASS", "All filters passed");
        for(auto & s : vcf.header.getSamples()) mw.header.addSample(s);
        for(auto & s : vcf.header.getSeqnames()) mw.header.addContig(s);
    }
    {
        fw.header.addFORMAT("GT", "1", "String", "Phased genotypes");
        fw.header.addFORMAT("GP", "3", "Float", "Genotype probabilities");
        fw.header.addFORMAT("DS", "1", "Float", "Genotype dosages");
        fw.header.addFILTER("PASS", "All filters passed");
        for(auto & s : vcf.header.getSamples()) fw.header.addSample(s);
        for(auto & s : vcf.header.getSeqnames()) fw.header.addContig(s);
    }
    BcfRecord war(mw.header);
    int N = vcf.nsamples, i{0};
    vector<int> gts, mgt(N * 2), fgt(N * 2);
    vector<float> mgp, fgp, mds, fds;
    const vector<char> phased(N, 1);
    war.setPhasing(phased);
    war.setQUAL('.');
    war.setFILTER("PASS");
    BcfRecord var(vcf.header);
    while(vcf.getNextVariant(var))
    {
        var.getGenotypes(gts);
        var.getFORMAT("MGP", mgp);
        var.getFORMAT("FGP", fgp);
        var.getFORMAT("MDS", mds);
        var.getFORMAT("FDS", fds);
        for(i = 0; i < N; i++)
        {
            mgt[i * 2 + 0] = gts[i * 3 + 0];
            mgt[i * 2 + 1] = gts[i * 3 + 1];
            fgt[i * 2 + 0] = gts[i * 3 + 0];
            fgt[i * 2 + 1] = gts[i * 3 + 2];
        }
        war.setCHR(var.CHROM().c_str());
        war.setPOS(var.POS());
        war.setRefAlt(var.REF() + "," + var.ALT());
        war.setID(var.ID().c_str());
        war.setGenotypes(mgt);
        war.setFORMAT("GP", mgp);
        war.setFORMAT("DS", mds);
        mw.writeRecord(war);
        war.setGenotypes(fgt);
        war.setFORMAT("GP", fgp);
        war.setFORMAT("DS", fds);
        fw.writeRecord(war);
    }
    mw.close();
    fw.close();

    return 0;
}
