#include "catch.hh"
#include "vcfpp.h"
#include <eigen3/Eigen/Dense>

using namespace vcfpp;
using namespace std;
using namespace Eigen;

TEST_CASE("Calculate the heterozygosity rate", "[bcf-eigen]")
{
    BcfReader vcf("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz", "-", "chr22:19000000-20000000");
    BcfRecord v(vcf.header);
    vector<int> gt;
    ArrayXi hetsum = ArrayXi::Zero(vcf.nsamples);
    while (vcf.getNextVariant(v)) {
        if (!v.isSNP()) continue; // skip other type of variants
        v.getGenotypes(gt);
        Map<const ArrayXXi> m(gt.data(), v.nploidy , gt.size() / v.nploidy); // read only
        hetsum += (m.row(0) - m.row(1)).abs(); // for diploid
    }
    cout << hetsum.mean() << endl;
}

TEST_CASE("Infer population structure with PCA", "[bcf-eigen]")
{
    BcfReader vcf("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz", "HG00096,HG00097,HG00099,HG00100,HG00101,HG00102,HG00103,HG00105,HG00106,HG00107", "chr22:19000000-20000000");
    BcfRecord v(vcf.header);
    int nsnps = 0, nsamples = vcf.nsamples;
    vector<bool> gt;
    vector<float> mat, gts(nsamples);
    while (vcf.getNextVariant(v)) {
        if (!v.isSNP()) continue; // skip other type of variants
        v.getGenotypes(gt);
        for (int i = 0; i < nsamples; i++) {
            gts[i] = gt[2 * i] + gt[2 * i + 1];
        }
        mat.insert(mat.end(), gts.begin(), gts.end());
        nsnps++;
    }
    Eigen::Map<Eigen::MatrixXf> M(mat.data(), nsamples, nsnps);
    Eigen::ArrayXf hwe = M.colwise().mean().array();
    hwe = (hwe * (1 - hwe/2)).sqrt() ;
    hwe = (hwe > 1e-9).select(hwe, 1); // in case denominator is smaller than 1e-9
    M = (-M).rowwise() + M.colwise().mean(); // centering by subtracting the mean
    M = (M.array().rowwise() / hwe.transpose()).matrix(); // standardize the matrix
    Eigen::JacobiSVD<Eigen::MatrixXf> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
    cout << svd.matrixU().leftCols(10) << endl; // save or print out top 10 PCs
    cout << svd.singularValues().array().square() / nsnps << endl; // print out eigenvalues
}
