// -*- compile-command: "g++ mspbwt.cpp -o mspbwt -std=c++17 -g -O3 -Wall -lhts -lz -lm -lbz2 -llzma -lcurl && ./mspbwt -i t.bcf -r chr22 > t" -*-
#include "../vcfpp.h"

#include <algorithm>
#include <bitset>
#include <map>
#include <numeric>
#include <set>


using namespace std;
using namespace vcfpp;

// const int B{sizeof(T) * 8};
static const int B = 32;

void mspbwt_encode(const std::string& vcffile, const std::string& samples, const std::string& region);
std::map<int, std::vector<int>> build_W(const vector<int>& x, int cap = 0);

int main(int argc, char* argv[])
{
    // ========= helper message and parameters parsing ============================
    std::vector<std::string> args(argv + 1, argv + argc);
    if (argc <= 1 || args[0] == "-h" || args[0] == "-help")
    {
        std::cout << "Author: Zilong-Li (zilong.dk@gmail.com)\n"
                  << "Description:\n"
                  << "     create DS tag for diploid samples given GP or GT tag from input vcf file\n\n"
                  << "Usage example:\n"
                  << "     " + (std::string)argv[0] + " -i in.bcf \n"
                  << "     " + (std::string)argv[0] + " -i in.bcf -o out.bcf -s ^S1,S2 -r chr1:1-1000 \n"
                  << "     bcftools view in.bcf | " + (std::string)argv[0] + " -i - -o out.bcf \n"
                  << "\nOptions:\n"
                  << "     -i    input vcf/bcf file\n"
                  << "     -o    ouput vcf/bcf file [stdout]\n"
                  << "     -s    list of samples to be included or excluded\n"
                  << "     -r    specific region to be included\n"
                  << std::endl;
        return 1;
    }
    std::string invcf, outvcf = "-", samples = "-", region = "";
    for (int i = 0; i < args.size(); i++)
    {
        if (args[i] == "-i")
            invcf = args[++i];
        if (args[i] == "-o")
            outvcf = args[++i];
        if (args[i] == "-s")
            samples = args[++i];
        if (args[i] == "-r")
            region = args[++i];
    }

    // ========= core calculation part ===========================================

    mspbwt_encode(invcf, samples, region);

    return 0;
}

void mspbwt_encode(const std::string& vcffile, const std::string& samples, const std::string& region)
{
    BcfReader vcf(vcffile, samples, region);
    BcfRecord var(vcf.header);
    vector<char> gt;
    uint64_t N{0}, M{0}, G{0}, k{0}, m{0}, i{0}; // N haplotypes, M SNPs, G Grids, k Current grid, m Current SNP
    N = vcf.nsamples * 2;
    while (vcf.getNextVariant(var))
    {
        if (!var.isSNP())
            continue; // skip other type of variants
        M++;
    }
    G = (M + B - 1) / B;
    cerr << N << "\t" << M << "\t" << G << endl;
    vector<vector<int>> X; // Grids x Haps
    X.resize(G, vector<int>(N));
    vcf.setRegion(region);
    while (vcf.getNextVariant(var))
    {
        if (!var.isSNP())
            continue; // skip other type of variants
        var.getGenotypes(gt);
        // update current grids
        for (i = 0; i < N; i++)
            X[k][i] = (X[k][i] << 1) | (gt[i] != 0);
        m++;
        if (m % B == 0)
            k++; // update next grid
    }
    cerr << m << "\t" << k << "\t" << B << endl;
    if (G == k + 1)
    {
        int pad = G * B - M;
        for (i = 0; i < N; i++)
            X[k][i] <<= pad;
    }
    else if (G == k)
    {
        cerr << "no need padding\n";
    }
    else
    {
        throw std::runtime_error("something wrong\n");
    }
    vector<vector<int>> A; // (Grids+1) x Haps
    A.resize(G + 1, vector<int>(N));
    // vector<std::map<int, std::vector<int>>> W(G);
    vector<int> a0(N), y0(N);
    for (k = 0; k < G; k++)
    {
        if (k == 0)
        {
            std::iota(a0.begin(), a0.end(), 0);
            A[k] = a0;
        }
        for (i = 0; i < N; i++)
            y0[i] = X[k][a0[i]];
        auto Wg = build_W(y0); // TODO: don't cap symbols when building indices
        // W[k] = Wg;          // save current W
        for (i = 0; i < N; i++)
        {
            A[k + 1][Wg[y0[i]][i] - 1] = a0[i];
        }
        a0 = A[k + 1]; // next run

        // for (const auto& [s, occ] : Wg)
        // {
        //     cout << s << ": ";
        //     for (auto& oi : occ)
        //     {
        //         cout << oi << "\t";
        //     }
        //     cout << endl;
        // }
        // for (k = 0; k < N; k++)
        //     cout << X[k][i] << "\t";
        // cout << endl;
    }

    for (auto & ai : A) {
        for (auto & v : ai) {
            cout << v << "\t";
        }
        cout << endl;
    }

}


std::map<int, std::vector<int>> build_W(const vector<int>& x, int cap)
{
    std::set<int> s(x.begin(), x.end()); // convert to set which unique sorted
    std::map<int, int> C;
    int c{0}, i{0}, n{0};
    for (auto& si : s)
    {
        c = 0;
        for (auto& xi : x)
        {
            if (xi < si)
                c++;
        }
        C[si] = c;
        if (++n == cap && cap > 0)
            break;
    }
    std::map<int, std::vector<int>> W;
    n = 0;
    for (auto& si : s)
    {
        W[si] = vector<int>(x.size());
        c = 0;
        for (i = 0; i < x.size(); i++)
        {
            if (x[i] == si)
                c++;
            W[si][i] = c + C[si];
        }
        if (++n == cap && cap > 0)
            break;
    }
    return W;
    // for (auto & i : s)
    //     cout << bitset<32>(i) << "\t";
    // cout << endl;
}
