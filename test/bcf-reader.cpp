#include "../vcfpp.h"
#include "catch.hh"
#include <exception>

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
    REQUIRE(br.getNextVariant(var) == true);
    var.getFORMAT("PL", pl);
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
    bw.writeLine("chr20\t2006060\trs146931526\tG\tA\t100\tPASS\tAF=0.02\tGT:GL\t0/"
                 "1:-323.03,-99.29,-802.53\t1/1:-133.03,-299.29,-902.53");
    bw.close();
    BcfReader br(vcffile);
    BcfRecord var(br.header);
    vector<float> gl;
    REQUIRE(br.getNextVariant(var) == true);
    var.getFORMAT("GL", gl);
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
    REQUIRE(br.getNextVariant(var) == true);
    auto l2 = var.asString();
    REQUIRE(l2 == l1 + "\n");
    var.getGenotypes(gt);
    for(auto g : gt) cout << g << endl;
}

TEST_CASE("parse EV in vcf - vector<string>", "[bcf-reader]")
{
    string vcffile{"test-GL.vcf.gz"};
    BcfWriter bw(vcffile, "VCF4.2");
    bw.header.addContig("chr20");
    bw.header.addINFO("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
    bw.header.addFORMAT("GT", "1", "String", "Genotype");
    bw.header.addFORMAT("EV", "G", "String", "Classes of evidence supporting final genotype");
    for(auto & s : {"id01", "id02"}) bw.header.addSample(s); // add 3 samples
    bw.writeLine("chr20\t2006060\trs146931526\tG\tA\t100\tPASS\tAF=0.02\tGT:EV\t0/1:RD\t1/1:SR,PE");
    bw.close();
    BcfReader br(vcffile);
    BcfRecord var(br.header);
    vector<string> ev;
    REQUIRE(br.getNextVariant(var) == true);
    var.getFORMAT("EV", ev);
    REQUIRE(ev[0] == "RD");
    REQUIRE(ev[1] == "SR,PE");
}

TEST_CASE("throw error when file is not valid", "[bcf-reader]")
{
    BcfReader br;
    CHECK_THROWS(br.open("no-test-GL.vcf.gz"));
    BcfReader bw;
    CHECK_THROWS(bw.open("ff://no-access.vcf.gz"));
}

TEST_CASE("throw error if the region is not valid", "[bcf-reader]")
{
    BcfReader * br = nullptr;
    try
    {
        br = new BcfReader("test-region.vcf.gz", "XXXX");
        delete br;
    }
    catch(exception & e)
    {
        cout << e.what();
    }
}

TEST_CASE("can check the status of a region", "[bcf-reader]")
{
    int status;
    BcfReader br("test-region.vcf.gz");
    status = br.getStatus("chr21:5030089-5030090");
    REQUIRE(status == 0); // valid but empty
    status = br.getStatus("chr21:5030089-");
    REQUIRE(status == 1); // valid and not empty
    status = br.getStatus("chr22:5030089");
    REQUIRE(status == -2); // not valid or not found in the VCF
    BcfReader br2("test-no-index.vcf.gz");
    status = br2.getStatus("chr21");
    REQUIRE(status == -1); // no index file found
}

TEST_CASE("can count the number of variants in a valid region", "[bcf-reader]")
{
    int nVariants = -1;
    BcfReader br("test-region.vcf.gz");
    nVariants = br.getVariantsCount("chr21:5030089-5030090");
    REQUIRE(nVariants == 0);
    nVariants = br.getVariantsCount("chr21:5030089-");
    REQUIRE(nVariants == 13);
}

TEST_CASE("can work with chr:pos form", "[bcf-reader]")
{
    int nVariants = -1;
    BcfReader br("test-region.vcf.gz");
    nVariants = br.getVariantsCount("chr21:5030088");
    REQUIRE(nVariants == 1);
    nVariants = br.getVariantsCount("chr21:5030087");
    REQUIRE(nVariants == 0);
    nVariants = br.getVariantsCount("chr21:5030089");
    REQUIRE(nVariants == 0);
}

TEST_CASE("subset both region and samples for bcf", "[bcf-reader]")
{
    BcfReader br("test.bcf", "1:10000-13000", "I1,I30");
    BcfRecord var(br.header);
    vector<float> ds;
    int n = 1;
    while(br.getNextVariant(var))
    {
        var.getFORMAT("DS", ds);
        if(n == 1) REQUIRE(ds == vector<float>{2, 2});
        if(n == 2) REQUIRE(ds == vector<float>{0.195, 0.000});
        if(n == 3) REQUIRE(ds == vector<float>{1.000, 0.000});
        if(n == 4) REQUIRE(ds == vector<float>{2.000, 0.036});
        n++;
    }
}
