#include "catch.hh"
#include "vcfpp.h"

using namespace vcfpp;
using namespace std;

TEST_CASE("make a vcf from scratch", "[bcf-writer]")
{
    BcfWriter bw("new.vcf.gz");
    bw.initalHeader("VCF4.3");
    bw.header.addContig("21");
    bw.header.addFORMAT("GT", "1", "String", "Genotype");
    bw.header.addFORMAT("PL", "G", "Integer", "List of Phred-scaled genotype likelihoods");
    bw.header.addFORMAT("DP", "1", "Integer", "Number of high-quality bases");
    bw.header.addINFO("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
    for (auto & s : {"S001", "S002", "S003", "S004"}) {
        bw.header.addSample(s);
    }
    REQUIRE(bw.header.nSamples() == 4);
    bw.writeHeader();

    BcfReader br("new.vcf");
}

// TODO missing INFO operations for variant;
void test_bcf_setvalues()
{
    BcfReader br2("test/index.vcf");
    BcfRecord v2(br2.header);
    br2.getNextVariant(v2);
    vector<float> b; v2.getINFO("QS", b);
    cout << b[0] << "\tqs\t" << b[1] << endl;
    cout << v2.asString() << endl;
    v2.removeINFO("DP");
    cout << v2.asString() << endl;
    br2.header.addINFO("Str", "1", "String", "this is a test for adding string in INFO");
    v2.setINFO("DP", 2); // int
    v2.setINFO("MQ0F", 2.2);        // double
    // v2.setInfo("MQ0F", (float)2.2); // float
    int a;
    v2.getINFO("DP", a);
    cout << v2.asString() << endl;
    v2.setINFO("Str", string{"S1S2"});
    string s; v2.getINFO("Str", s);
    cout << s << endl;
    vector<int> pl;
    v2.getFormat("PL", pl);
    for (auto i : pl) { cout << i << ",";}
    cout << endl;
}
void test_cyvcf2_paper_case()
{
    BcfReader br("ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz");
    BcfRecord v(br.header);
    vector<char> gt;
    int n = 0;
    while (br.getNextVariant(v))
    {
        if (v.isIndel()) continue;
        if (v.isSNP())
        {
            cout << n << endl;
            v.getGenotypes(gt);
            Eigen::Map<Eigen::MatrixX<char>> m(gt.data(), 2, gt.size()/2);
            // cout << m.rowwise().sum() << endl;
        }
        if (v.isIndel()) cout << "indel line " << n + 1 << endl;
        if (v.isDeletion()) cout << "deletion line " << n + 1 << endl;
        n++;
        if (n > 20) break;
    }
}
