#include "../vcfpp.h"
#include "catch.hh"

using namespace vcfpp;
using namespace std;

TEST_CASE("parse vcf.gz with no missing genotypes diploid -- vector<bool>", "[bcf-reader]")
{
    BcfReader br("test-vcf-read.vcf.gz");
    BcfRecord v(br.header);
    vector<bool> gt;
    int n = 0;
    while(br.getNextVariant(v))
    {
        v.getGenotypes(gt);
        for(auto g : gt) cout << g << endl;
        if(!v.allPhased()) n++;
    }
    REQUIRE(n == 10);
}

TEST_CASE("parse vcf.gz with no missing genotypes diploid -- vector<char>", "[bcf-reader]")
{
    BcfReader br("test-vcf-read.vcf.gz");
    BcfRecord v(br.header);
    vector<char> gt;
    int n = 0;
    while(br.getNextVariant(v))
    {
        v.getGenotypes(gt);
        for(auto g : gt) cout << g << endl;
        if(!v.allPhased()) n++;
    }
    REQUIRE(n == 10);
}

TEST_CASE("parse vcf.gz with missing genotypes diploid - vector<int>", "[bcf-reader]")
{
    BcfReader br("test-vcf-read.vcf.gz");
    BcfRecord v(br.header);
    vector<int> gt;
    int n = 0;
    while(br.getNextVariant(v))
    {
        v.getGenotypes(gt);
        for(auto g : gt) cout << g << endl;
        if(!v.allPhased()) n++;
    }
    REQUIRE(n == 10);
}

TEST_CASE("parse vcf with missing genotypes diploid - vector<int>", "[bcf-reader]")
{
    BcfReader br("test-vcf-read.vcf");
    BcfRecord v(br.header);
    vector<int> gt;
    int n = 0;
    while(br.getNextVariant(v))
    {
        v.getGenotypes(gt);
        for(auto g : gt) cout << g << endl;
        if(!v.isNoneMissing()) n++;
    }
    REQUIRE(n == 1);
}

TEST_CASE("parse PL in vcf - vector<int>", "[bcf-reader]")
{
    string vcffile{"test-GL.vcf.gz"};
    BcfWriter bw(vcffile, "VCF4.2");
    bw.header.addContig("chr20"); 
    bw.header.addINFO("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
    bw.header.addFORMAT("GT", "1", "String", "Genotype");
    bw.header.addFORMAT("PL", "G", "Integer", "List of Phred-scaled genotype likelihoods");
    for(auto & s : {"id01", "id02"}) bw.header.addSample(s); // add 3 samples
    bw.writeLine("chr20\t2006060\trs146931526\tG\tA\t100\tPASS\tAF=0.02\tGT:PL\t0/0:0,9,75\t0/0:0,3,30");
    bw.close();
    BcfReader br(vcffile);
    BcfRecord var(br.header);
    vector<int> pl;
    REQUIRE(br.getNextVariant(var)==true);
    var.getFORMAT("PL",pl);
    for(auto g : pl) cout << g << endl;
}

TEST_CASE("parse GL in vcf - vector<float>", "[bcf-reader]")
{
    string vcffile{"test-GL.vcf.gz"};
    BcfWriter bw(vcffile, "VCF4.2");
    bw.header.addContig("chr20"); 
    bw.header.addINFO("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
    bw.header.addFORMAT("GT", "1", "String", "Genotype");
    bw.header.addFORMAT("GL", "G", "Float", "List of log scale genotype likelihoods");
    for(auto & s : {"id01", "id02"}) bw.header.addSample(s); // add 3 samples
    bw.writeLine("chr20\t2006060\trs146931526\tG\tA\t100\tPASS\tAF=0.02\tGT:GL\t0/1:-323.03,-99.29,-802.53\t1/1:-133.03,-299.29,-902.53");
    bw.close();
    BcfReader br(vcffile);
    BcfRecord var(br.header);
    vector<float> gl;
    REQUIRE(br.getNextVariant(var)==true);
    var.getFORMAT("GL",gl);
    for(auto g : gl) cout << g << endl;
}

TEST_CASE("parse vcf with multialleles - vector<int>", "[bcf-reader]")
{
    string vcffile{"test-multialleles.vcf.gz"};
    BcfWriter bw(vcffile, "VCF4.2");
    bw.header.addContig("chr20"); 
    bw.header.addINFO("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
    bw.header.addFORMAT("GT", "1", "String", "Genotype");
    for(auto & s : {"id01", "id02", "id03"}) bw.header.addSample(s); // add 3 samples
    string l1 = "chr20\t2006060\trs146931526\tG\tA,T\t100\tPASS\tAF=0.02\tGT\t1|2\t1|1\t0|2";
    bw.writeLine(l1);
    bw.close();
    BcfReader br("test-multialleles.vcf.gz");
    BcfRecord var(br.header);
    vector<int> gt;
    REQUIRE(br.getNextVariant(var)==true);
    auto l2 = var.asString();
    REQUIRE(l2 == l1 + "\n");
    var.getGenotypes(gt);
    for(auto g : gt) cout << g << endl;
}
