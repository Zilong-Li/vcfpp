#include "../vcfpp.h"
#include "catch.hh"
#include <exception>

using namespace vcfpp;
using namespace std;

TEST_CASE("can get type of INFO field from header", "[bcf-header]")
{
    BcfReader br("test-region.vcf.gz");

    REQUIRE(br.header.getInfoType("AA") == -1); // no existing
    REQUIRE(br.header.getInfoType("AN") == 1); // int
    REQUIRE(br.header.getInfoType("AF") == 2); // float
    REQUIRE(br.header.getInfoType("VariantType") == 3); // string
    REQUIRE(br.header.getInfoType("POSITIVE_TRAIN_SITE") == 0); //  doesn't support flag type, give it error code
}

TEST_CASE("can get type of FORMAT field from header", "[bcf-header]")
{
    BcfReader br("test-region.vcf.gz");

    REQUIRE(br.header.getFormatType("AA") == 0); // no existing
    REQUIRE(br.header.getFormatType("AD") == 1); // int
    REQUIRE(br.header.getFormatType("AB") == 2); // float
    REQUIRE(br.header.getFormatType("GT") == 3); // string
    REQUIRE(br.header.getFormatType("PGT") == 3); // string
}
