#include "catch.hh"
#include "vcfpp.h"

using namespace vcfpp;
using namespace std;

TEST_CASE("Parsing BCF with subset samples in target region", "[bcf-reader]")
{
    BcfReader br("test-vcf-read.vcf.gz");
    BcfRecord v(br.header);
    vector<bool> gt;
    int n = 0;
    while (br.getNextVariant(v))
    {
        v.getGenotypes(gt);
        if(!v.allPhased()) n++;
    }
    REQUIRE(n == 10);
}
