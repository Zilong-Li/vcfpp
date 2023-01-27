#include "../vcfpp.h"
#include "catch.hh"

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
    for (auto& s : {"S001", "S002", "S003", "S004"})
    {
        bw.header.addSample(s);
    }
    REQUIRE(bw.header.nSamples() == 4);
    bw.writeHeader();
}

TEST_CASE("Write BCF with custome header and variants", "[bcf-writer]")
{
    BcfWriter bw("test.bcf");
    bw.header.addFORMAT("GT", "1", "String", "Genotype");
    bw.header.addINFO("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
    bw.header.addContig("chr20");
    bw.header.addSample("NA12878");
    bw.writeLine("chr20\t2006060\trs146931526\tG\tC\t100\tPASS\tAF=0.000998403\tGT\t1|0");
}

TEST_CASE("Write VCF with custome header and variants", "[bcf-writer]")
{
    BcfWriter bw("test.vcf", "VCF4.3");
    bw.header.addFORMAT("GT", "1", "String", "Genotype");
    bw.header.addINFO("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
    bw.header.addContig("chr20");
    bw.header.addSample("NA12878");
    bw.writeLine("chr20\t2006060\trs146931526\tG\tC\t100\tPASS\tAF=0.000998403\tGT\t1|0");
}


TEST_CASE("Write VCF by copying header from another VCF", "[bcf-writer]")
{
    BcfReader br("../htslib/test/test-vcf-hdr-in.vcf");
    BcfRecord v(br.header);
    BcfWriter bw("test.vcf", br.header);
    bw.writeHeader();
    br.getNextVariant(v);
    bw.writeRecord(v);
}
