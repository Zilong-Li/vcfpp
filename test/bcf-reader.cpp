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

TEST_CASE("parse vcf with multialleles - vector<int>", "[bcf-reader]")
{
    BcfWriter bw("test-multialleles.vcf", "VCF4.3");
    bw.header.addContig("chr20"); 
    bw.header.addFORMAT("GT", "2", "String", "Genotype");
    bw.header.addINFO("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
    for(auto & s : {"id01", "id02", "id03"}) bw.header.addSample(s); // add 3 samples
    bw.writeLine("chr20\t2006060\trs146931526\tG\tA,T\t100\tPASS\tAF=0.02\tGT\t1|2\t1|1\t0|2");
    bw.close();
    BcfReader br("test-multialleles.vcf");
    BcfRecord var(br.header);
    vector<int> gt;
    while(br.getNextVariant(var))
    {
        var.getGenotypes(gt);
        for(auto g : gt) cout << g << endl;
    }
}
