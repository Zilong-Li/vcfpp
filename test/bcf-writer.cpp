#include "catch.hh"
#include "../vcfpp.h"

using namespace vcfpp;
using namespace std;

TEST_CASE("make a vcf from scratch", "[bcf-writer]")
{
    BcfWriter bw("new.bcf.gz", "VCF4.3");
    bw.header.addContig("21");
    bw.header.addFORMAT("GT", "1", "String", "Genotype");
    bw.header.addFORMAT("PL", "G", "Integer", "List of Phred-scaled genotype likelihoods");
    bw.header.addFORMAT("DP", "1", "Integer", "Number of high-quality bases");
    bw.header.addINFO("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
    for (auto & s : {"S001", "S002", "S003", "S004"}) {
        bw.header.addSample(s);
    }
    REQUIRE(bw.header.nSamples() == 4);
    bw.writeHeader();
}
