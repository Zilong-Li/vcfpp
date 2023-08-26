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
