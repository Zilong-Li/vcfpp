#include "../vcfpp.h"
#include "catch.hh"

using namespace vcfpp;
using namespace std;

TEST_CASE("set genotypes without missing ", "[bcf-variant]")
{
    BcfReader br("test-vcf-read.vcf.gz");
    BcfRecord v(br.header);
    vector<int> gt;
    vector<int> gi{1, 1};
    while(br.getNextVariant(v))
    {
        v.getGenotypes(gt);
        v.setGenotypes(gi);
        v.getGenotypes(gt);
        for(size_t i = 0; i < gt.size(); i++) REQUIRE(gt[i] == gi[i]);
    }
}

TEST_CASE("set genotypes with missing ", "[bcf-variant]")
{
    BcfReader br("test-vcf-read.vcf.gz");
    BcfRecord v(br.header);
    vector<int> gt;
    vector<int> gi{1, -9}; // -9 as missing
    vector<int> g2{-2147483648, -9}; // also, bcf_int32_missing (R_NA) as missing
    // br.getNextVariant(v);
    while(br.getNextVariant(v))
    {
        v.getGenotypes(gt);
        v.setGenotypes(gi);
        v.getGenotypes(gt);
        for(size_t i = 0; i < gt.size(); i++) REQUIRE(gt[i] == gi[i]);
        v.setGenotypes(g2);
        v.getGenotypes(gt);
        for(size_t i = 0; i < gt.size(); i++) REQUIRE(gt[i] == -9);
    }
}
