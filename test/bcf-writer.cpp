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
    for(auto & s : {"S001", "S002", "S003", "S004"})
    {
        bw.header.addSample(s);
    }
    REQUIRE(bw.header.nSamples() == 4);
    bw.writeHeader();
}

TEST_CASE("Write BCF with custome header and variants", "[bcf-writer]")
{
    BcfWriter bw("test.write.vcf.gz");
    bw.header.addFORMAT("GT", "1", "String", "Genotype");
    bw.header.addINFO("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
    bw.header.addContig("chr20");
    bw.header.addSample("NA12878");
    bw.writeLine("chr20\t2006060\trs146931526\tG\tC\t100\tPASS\tAF=0.000998403\tGT\t1|0");
    bw.close();
}

TEST_CASE("can use header from another bcf.gz", "[bcf-writer]")
{
    BcfReader br("test-vcf-read.bcf.gz");
    BcfRecord v(br.header);
    BcfWriter bw("test.write.bcf.gz", br.header);
    bw.writeHeader();
    br.getNextVariant(v);
    bw.writeRecord(v);
}

TEST_CASE("init header twice", "[bcf-writer]")
{
    BcfReader br("test-vcf-read.vcf");
    BcfRecord v(br.header);
    BcfWriter bw("test.write.vcf.gz", "VCF4.3");
    bw.initalHeader(br.header);
    bw.writeHeader();
    br.getNextVariant(v);
    bw.writeRecord(v);
}

TEST_CASE("Write VCF by copying header of another VCF and modifying it", "[bcf-writer]")
{
    BcfReader br("test-vcf-read.bcf");
    BcfWriter bw;
    bw.open("test-vcf-write.vcf");
    bw.copyHeader("test-vcf-read.bcf");
    bw.header.addINFO("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
    bw.header.updateSamples("X001");
    BcfRecord v(bw.header);
    br.getNextVariant(v);
    v.setINFO("AF", 0.2);
    bw.writeRecord(v);
    bw.close();
}
