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
    if(br.getNextVariant(v))
    {
        v.getGenotypes(gt);
        for(auto g : gt) cerr << g << endl;
        v.setGenotypes(gi);
        v.getGenotypes(gt);
        for(auto g : gt) cerr << g << endl;
    }
}

TEST_CASE("set genotypes with missing ", "[bcf-variant]")
{
    BcfReader br("test-vcf-read.vcf.gz");
    BcfRecord v(br.header);
    vector<int> gt;
    vector<int> gi{1, -9};
    vector<int> g2{-2147483648, -9};
    br.getNextVariant(v);
    if(br.getNextVariant(v))
    {
        v.getGenotypes(gt);
        for(auto g : gt) cerr << g << endl;
        v.setGenotypes(gi);
        v.getGenotypes(gt);
        for(auto g : gt) cerr << g << endl;
        v.setGenotypes(g2);
        v.getGenotypes(gt);
        for(auto g : gt) cerr << g << endl;
    }
}
