#include "catch.hh"
#include "vcfpp.h"
#include <iostream>

using namespace vcfpp;
using namespace std;

TEST_CASE("Parsing VCF with specific tag", "[bcf-reader]")
{
    BcfReader br("test/test-vcf-read.vcf");
    BcfRecord v(br.header);
    vector<float> ad;
    vector<int> gq;
    vector<char> gatk;
    vector<bool> gt;
    int n = 0;
    while (br.GetNextVariant(v))
    {
        n++;
        v.GetGenotypes(gt);
        for (auto i : gt) {
            cout << i << ":";
        }
        cout << endl;

        v.GetFormat(ad, "AD");
        for (auto i : ad) {
            cout << i << ":";
        }
        cout << endl;

        v.GetFormat(gq, "GQ");
        for (auto i : gq) {
            cout << i << ":";
        }
        cout << endl;

        v.GetFormat(gatk, "GATK");
        for (int i = 0; i < v.header->nsamples; i++) {
            cout << string(gatk.begin()+ i*v.shape1, gatk.begin() + i * v.shape1 + v.shape1 - 1) << ":";
        }
        cout << endl;

    }
    std::shared_ptr<BcfHeader> h = std::make_shared<BcfHeader>(br.header);
    h->RemoveContig("1");
    h->RemoveContig("16");
    h->RemoveFilter("LowQual");
    h->RemoveInfo("TRAILING");
    cout << h->AsString() << endl;

    BcfReader br2("test/index.vcf");
    BcfRecord v2(br2.header);
    br2.GetNextVariant(v2);
    v2.RemoveInfo<int>("DP");
    br2.header.AddInfo("Str", "1", "String", "this is a test for adding string in INFO");
    v2.SetInfo("Str", string{"S1S2"});
    v2.SetInfo("Str", string{"str"});
    v2.SetInfo("DP", vector<int>{1,2});
    v2.SetInfo("DP", 2);

    cout << br2.header.AsString() << endl;
    cout << v2.GetInfo<float>("DP") << endl;
    cout << v2.AsString() << endl;

    REQUIRE(n == 10);
}

TEST_CASE("Parsing VCF", "[bcf-reader]")
{
    BcfReader br("test/chr20.2000001.2100000.vcf.gz");
    BcfRecord v(br.header);
    vector<bool> gt;
    int n = 0;
    while (br.GetNextVariant(v))
    {
        v.GetGenotypes(gt);
        REQUIRE(gt.size() == 2504 * 2);
        if(v.isAllPhased) n++;
    }
    REQUIRE(n == 3786);
}

TEST_CASE("Parsing VCF with subset samples", "[bcf-reader]")
{
    BcfReader br("test/chr20.2000001.2100000.vcf.gz", "HG00107,HG00108,HG00109,HG00110,HG00111,HG00112,HG00113,HG00114,HG00115,HG00116");
    BcfRecord v(br.header);
    vector<bool> gt;
    int n = 0;
    while (br.GetNextVariant(v))
    {
        v.GetGenotypes(gt);
        REQUIRE(gt.size() == 20);
        if(v.isAllPhased) n++;
    }
    REQUIRE(n == 3786);
}

TEST_CASE("Parsing VCF in target region", "[bcf-reader]")
{
    BcfReader br("test/chr20.2000001.2100000.vcf.gz", "-", "chr20:2006060");
    BcfRecord v(br.header);
    vector<bool> gt;
    int n = 0;
    while (br.GetNextVariant(v))
    {
        v.GetGenotypes(gt);
        REQUIRE(gt.size() == 2504 * 2);
        if(v.isAllPhased) n++;
    }
    REQUIRE(n == 3243);
}

TEST_CASE("Parsing VCF with subset samples in target region", "[bcf-reader]")
{
    BcfReader br("test/chr20.2000001.2100000.vcf.gz", "HG00107,HG00108,HG00109,HG00110,HG00111,HG00112,HG00113,HG00114,HG00115,HG00116", "chr20:2006060-");
    BcfRecord v(br.header);
    vector<bool> gt;
    int n = 0;
    while (br.GetNextVariant(v))
    {
        v.GetGenotypes(gt);
        REQUIRE(gt.size() == 20);
        if(v.isAllPhased) n++;
    }
    REQUIRE(n == 3243);
}

TEST_CASE("Parsing BCF", "[bcf-reader]")
{
    BcfReader br("test/chr20.2000001.2100000.bcf.gz");
    BcfRecord v(br.header);
    vector<int> gt;
    int n = 0;
    while (br.GetNextVariant(v))
    {
        v.GetGenotypes(gt);
        REQUIRE(gt.size() == 2504 * 2);
        if(v.isAllPhased) n++;
    }
    REQUIRE(n == 3786);
}

TEST_CASE("Parsing BCF in target region", "[bcf-reader]")
{
    BcfReader br("test/chr20.2000001.2100000.vcf.gz", "-", "chr20");
    BcfRecord v(br.header);
    vector<int> gt;
    int n = 0;
    while (br.GetNextVariant(v))
    {
        v.GetGenotypes(gt);
        REQUIRE(gt.size() == 2504 * 2);
        if(v.isAllPhased) n++;
    }
    REQUIRE(n == 3786);
}

TEST_CASE("Parsing BCF with subset samples", "[bcf-reader]")
{
    BcfReader br("test/chr20.2000001.2100000.bcf.gz", "HG00107,HG00108,HG00109,HG00110,HG00111,HG00112,HG00113,HG00114,HG00115,HG00116");
    BcfRecord v(br.header);
    vector<bool> gt;
    int n = 0;
    while (br.GetNextVariant(v))
    {
        v.GetGenotypes(gt);
        REQUIRE(gt.size() == 20);
        if(v.isAllPhased) n++;
    }
    REQUIRE(n == 3786);
}

TEST_CASE("Parsing BCF with subset samples in target region", "[bcf-reader]")
{
    BcfReader br("test/chr20.2000001.2100000.bcf.gz", "HG00107,HG00108,HG00109,HG00110,HG00111,HG00112,HG00113,HG00114,HG00115,HG00116", "chr20:2006060-2010000");
    BcfRecord v(br.header);
    vector<bool> gt;
    int n = 0;
    while (br.GetNextVariant(v))
    {
        v.GetGenotypes(gt);
        REQUIRE(gt.size() == 20);
        if(v.isAllPhased) n++;
    }
    REQUIRE(n == 129);
}
