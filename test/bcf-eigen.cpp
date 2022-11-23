#include "catch.hh"
#include "../vcfpp.h"
#include <eigen3/Eigen/Dense>

using namespace vcfpp;
using namespace std;
using namespace Eigen;

TEST_CASE("Build PBWT from VCF", "[bcf-eigen]")
{
    BcfReader vcf("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
                  "1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
                  "-", "chr22:16050075-16050655");
    BcfRecord v(vcf.header);
    vector<bool> gt;
    vector<vector<bool>> x;
    int nsnps = 0, nsamples = vcf.nsamples;
    while (vcf.getNextVariant(v))
    {
        if (!v.isSNP())
            continue; // skip other type of variants
        v.getGenotypes(gt);
        x.push_back(gt);
        nsnps++;
    }
    int N = nsnps, M = nsamples * 2;
    vector<vector<int>> a(N, vector<int>(M));     // this is the index matrix to be retured
    vector<int> a0(M), a1(M), d0(M), d1(M), d(M); // note: to output divergence matrix, make d as two dimensional N x M.
    int i = 0, j = 0, u_ = 0, v_ = 0, p = 0, q = 0, d_, a_;
    for (j = 0; j < N; j++)
    {
        for (i = 0; i < M; i++)
        {
            d_ = j > 0 ? d[i] : 0;
            a_ = j > 0 ? a[j - 1][i] : i;
            p = max(p, d_);
            q = max(q, d_);
            if (x[j][a_])
            {
                a1[v_] = a_;
                d1[v_] = q;
                v_++;
                q = 0;
            }
            else
            {
                a0[u_] = a_;
                d0[u_] = p;
                u_++;
                p = 0;
            }
        }
        for (i = 0; i < M; i++)
        {
            if (i < u_)
            {
                a[j][i] = a0[i];
                d[i] = d0[i];
            }
            else
            {
                a[j][i] = a1[i - u_];
                d[i] = d1[i - u_];
            }
        }
    }
    cout << "done\n" << endl;
}

TEST_CASE("Calculate the heterozygosity rate", "[bcf-eigen]")
{
    BcfReader vcf("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
                  "-", "chr22:19000000-20000000");
    BcfRecord v(vcf.header);
    vector<int> gt;
    ArrayXi hetsum = ArrayXi::Zero(vcf.nsamples);
    while (vcf.getNextVariant(v))
    {
        if (!v.isSNP())
            continue; // skip other type of variants
        v.getGenotypes(gt);
        Map<const ArrayXXi> m(gt.data(), v.nploidy, gt.size() / v.nploidy); // read only
        hetsum += (m.row(0) - m.row(1)).abs();                              // for diploid
    }
    cout << hetsum.mean() << endl;
}

TEST_CASE("Infer population structure with PCA", "[bcf-eigen]")
{
    BcfReader vcf("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
                  "HG00096,HG00097,HG00099,HG00100,HG00101,HG00102,HG00103,HG00105,HG00106,HG00107", "chr22:19000000-20000000");
    BcfRecord v(vcf.header);
    int nsnps = 0, nsamples = vcf.nsamples;
    vector<bool> gt;
    vector<float> mat, gts(nsamples);
    while (vcf.getNextVariant(v))
    {
        if (!v.isSNP())
            continue; // skip other type of variants
        v.getGenotypes(gt);
        for (int i = 0; i < nsamples; i++)
        {
            gts[i] = gt[2 * i] + gt[2 * i + 1];
        }
        mat.insert(mat.end(), gts.begin(), gts.end());
        nsnps++;
    }
    Eigen::Map<Eigen::MatrixXf> M(mat.data(), nsamples, nsnps);
    Eigen::ArrayXf hwe = M.colwise().mean().array();
    hwe = (hwe * (1 - hwe / 2)).sqrt();
    hwe = (hwe > 1e-9).select(hwe, 1);                    // in case denominator is smaller than 1e-9
    M = (-M).rowwise() + M.colwise().mean();              // centering by subtracting the mean
    M = (M.array().rowwise() / hwe.transpose()).matrix(); // standardize the matrix
    Eigen::JacobiSVD<Eigen::MatrixXf> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
    cout << svd.matrixU().leftCols(10) << endl;                    // save or print out top 10 PCs
    cout << svd.singularValues().array().square() / nsnps << endl; // print out eigenvalues
}
