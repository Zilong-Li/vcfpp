// -*- compile-command: "g++ mspbwt.cpp -std=c++17 -g -O3 -Wall -lhts -lz -lm -lbz2 -llzma -lcurl && ./a.out -p panel.vcf.gz -k 1 -r chr22 > t" -*-
#include "../vcfpp.h"

#include <algorithm>
#include <bitset>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <unordered_map>


using namespace std;
using namespace vcfpp;

template <typename T>
T reverseBits(T n, size_t B = sizeof(T) * 8)
{
    assert(B <= std::numeric_limits<T>::digits);
    T rv = 0;
    for (size_t i = 0; i < B; ++i, n >>= 1)
        rv = (rv << 1) | (n & 0x01);
    return rv;
}

using IntGridVec = vector<uint32_t>;
using IntSet = set<uint32_t, less<uint32_t>>;
using IntMap = unordered_map<uint32_t, uint32_t>;
using IntVecMap = unordered_map<uint32_t, std::vector<uint32_t>>;
using SymbolIdxMap = map<uint32_t, uint32_t, less<uint32_t>>;
using WgSymbolMap = map<uint32_t, SymbolIdxMap, less<uint32_t>>;

vector<uint32_t> encodeZgrid(const vector<bool>& z, int G);
vector<bool> randhapz(uint64_t M);
IntMap build_C(const IntGridVec& x, const IntSet& s);
WgSymbolMap save_W(const IntGridVec& x, const IntSet& s);
IntVecMap build_W(const IntGridVec& x, const IntSet& s, const IntMap& C);
void mspbwt(const std::string& vcfpanel, const std::string& samples, const std::string& region, int q, int ki);

void coutZYk(const vector<IntGridVec>& X, const vector<vector<int>>& A, const IntGridVec& zg, const vector<int>& az, int k, bool bit = 1);

int main(int argc, char* argv[])
{
    // ========= helper message and parameters parsing ============================
    std::vector<std::string> args(argv + 1, argv + argc);
    if (argc <= 1 || args[0] == "-h" || args[0] == "-help")
    {
        std::cout << "Author: Zilong-Li (zilong.dk@gmail.com)\n"
                  << "Usage example:\n"
                  << "     " + (std::string)argv[0] + " -p panel.vcf.gz -q 10 -r chr22 \n"
                  << "\nOptions:\n"
                  << "     -p    vcf/bcf file of reference panel\n"
                  << "     -r    chromosome to be included\n"
                  << "     -q    first n numer of samples as query [1]\n"
                  << "     -k    Y[k] to show. [1]\n"
                  << std::endl;
        return 1;
    }
    std::string vcfpanel, vcfquery, outvcf = "-", samples = "-", region = "";
    int q{1}, k{1};
    for (int i = 0; i < args.size(); i++)
    {
        if (args[i] == "-p")
            vcfpanel = args[++i];
        if (args[i] == "-o")
            outvcf = args[++i];
        if (args[i] == "-r")
            region = args[++i];
        if (args[i] == "-q")
            q = stoi(args[++i]);
        if (args[i] == "-k")
            k = stoi(args[++i]);
    }

    // ========= core calculation part ===========================================

    mspbwt(vcfpanel, samples, region, q, k);

    return 0;
}

vector<uint32_t> encodeZgrid(const vector<bool>& z, int G)
{
    vector<uint32_t> zg(G);
    const int B = 32;
    size_t m{0}, k{0}, M{z.size()};
    for (m = 0; m < M; m++)
    {
        zg[k] = (zg[k] << 1) | (z[m] != 0);
        if ((m + 1) % B == 0)
        {
            zg[k] = reverseBits(zg[k]);
            k++;
        }
    }
    if (G == k + 1)
    {
        zg[k] <<= G * B - M; // padding 0s
        zg[k] = reverseBits(zg[k]);
    }
    else if (G == k)
    {
        cerr << "no need padding\n";
    }
    else
    {
        throw std::runtime_error("something wrong\n");
    }
    return zg;
}

void coutZYk(const vector<IntGridVec>& X, const vector<vector<int>>& A, const IntGridVec& zg, const vector<int>& az, int k, bool bit)
{
    const size_t B = sizeof(X[0][0]) * 8;
    for (size_t i = 0; i < X[0].size(); i++)
    {
        for (int j = 0; j <= k + 1; j++)
        {
            auto rb = reverseBits(X[j][A[k + 1][i]]);
            if (bit)
                cout << std::bitset<B>(rb) << " ";
            else
                cout << rb << " ";
        }
        cout << endl;

        if (i == az[k])
        {
            // print out original Z bits
            cout << "========= zg is inserting here ========  k-1=" << k << ", az[k-1]=" << az[k] << endl;
            for (int j = 0; j <= k + 1; j++)
            {
                auto rb = reverseBits(zg[j]);
                if (bit)
                    cout << std::bitset<B>(rb) << " ";
                else
                    cout << rb << " ";
            }
            cout << endl;
            cout << "========= zg is inserting here ========  k-1=" << k << ", az[k]=" << az[k + 1] << endl;
        }
    }
}


vector<bool> randhapz(uint64_t M)
{
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dist(1, 6);
    vector<bool> z(M);
    for (uint64_t i = 0; i < M; i++)
    {
        z[i] = dist(gen) > 3;
    }
    return z;
}


IntMap build_C(const IntGridVec& x, const IntSet& s)
{
    uint32_t c{0}, n{0}, cap{0};
    IntMap C;
    for (const auto& si : s)
    {
        c = 0;
        for (const auto& xi : x)
        {
            if (xi < si)
                c++;
        }
        C[si] = c;
        if (++n == cap && cap > 0)
            break;
    }

    return C;
}

WgSymbolMap save_W(const IntGridVec& x, const IntSet& s)
{
    size_t i{0}, k{0};
    WgSymbolMap W;
    for (const auto& si : s)
    {
        k = 0;
        for (i = 0; i < x.size(); i++)
        {
            if (x[i] == si)
                W[si][i] = k++; // k: how many occurrences before i
        }
    }

    return W;
}

IntVecMap build_W(const IntGridVec& x, const IntSet& s, const IntMap& C)
{
    uint32_t c{0}, i{0};
    IntVecMap W;
    for (const auto& si : s)
    {
        W[si] = vector<uint32_t>(x.size());
        c = 0;
        for (i = 0; i < x.size(); i++)
        {
            if (x[i] == si)
                c++;
            W[si][i] = c + C.at(si);
        }
    }

    return W;
}

void mspbwt(const std::string& vcfpanel, const std::string& samples, const std::string& region, int q, int ki)
{
    const int B = 32;
    BcfReader vcf(vcfpanel, samples, region);
    BcfRecord var(vcf.header);
    // N haplotypes, M SNPs, G Grids, k Current grid, m Current SNP
    uint64_t Nq{0}, Np{0}, M{0}, G{0}, k{0}, m{0}, i{0}, j{0}, w{0};
    Nq = q * 2;
    Np = vcf.nsamples * 2 - Nq;
    M = vcf.get_region_records(region);
    G = (M + B - 1) / B;
    cerr << "Haps(Np):" << Np << "\tHaps(Nq):" << Nq << "\tSNPs(M):" << M << "\tGrids(G):" << G << "\tInt(B):" << B << endl;
    vector<IntGridVec> X, Z; // Grids x Haps
    X.resize(G, IntGridVec(Np));
    Z.resize(Nq, IntGridVec(G));
    vcf.setRegion(region); // seek back to region
    vector<bool> gt;
    while (vcf.getNextVariant(var))
    {
        var.getGenotypes(gt);
        // update current grids
        for (i = 0; i < gt.size(); i++)
        {
            if (i < Nq)
            {
                j = i;
                Z[j][k] = (Z[j][k] << 1) | (gt[i] != 0);
            }
            else
            {
                j = i - Nq;
                X[k][j] = (X[k][j] << 1) | (gt[i] != 0);
            }
        }
        m++;
        if (m % B == 0)
        {
            for (i = 0; i < gt.size(); i++)
            {
                if (i < Nq)
                {
                    j = i;
                    Z[j][k] = reverseBits(Z[j][k]); // reverset bits
                }
                else
                {
                    j = i - Nq;
                    X[k][j] = reverseBits(X[k][j]); // reverset bits
                }
            }
            k++; // update next grid
        }
    }
    if (G == k + 1)
    {
        for (i = 0; i < gt.size(); i++)
        {
            if (i < Nq)
            {
                j = i;
                Z[j][k] <<= G * B - M;
                Z[j][k] = reverseBits(Z[j][k]); // reverset bits
            }
            else
            {
                j = i - Nq;
                X[k][j] <<= G * B - M;
                X[k][j] = reverseBits(X[k][j]); // reverset bits
            }
        }
    }
    else if (G == k)
    {
        cerr << "no need padding\n";
    }
    else
    {
        throw std::runtime_error("something wrong\n");
    }


    vector<WgSymbolMap> W(G);
    vector<IntMap> C(G);
    IntGridVec y0(Np);
    vector<vector<int>> A; // (Grids+1) x Haps, we use int here so that R can take it
    A.resize(G + 1, vector<int>(Np));
    vector<int> a0(Np);
    IntSet symbols;
    for (k = 0; k < G; k++)
    {
        if (k == 0)
        {
            std::iota(a0.begin(), a0.end(), 0);
            A[k] = a0;
        }
        for (i = 0; i < Np; i++)
        {
            y0[i] = X[k][a0[i]];
            symbols.insert(y0[i]);
        }
        C[k] = build_C(y0, symbols);
        // here save current W, which differs from build whole W table
        W[k] = save_W(y0, symbols);
        // without Wg=build_W (S x N) is slow but memory efficient
        for (i = 0; i < Np; i++)
        {
            w = C[k][y0[i]];
            if (i >= W[k][y0[i]].begin()->first)
            {
                for (const auto& [sidx, rank] : W[k][y0[i]])
                {
                    if (sidx == i)
                    {
                        w += rank + 1;
                        break;
                    }
                }
            }
            A[k + 1][w - 1] = a0[i];
        }
        // next run
        a0 = A[k + 1];
        symbols.clear();
    }

    // finially insert zg back to A at each grid
    auto zg = Z[1];
    vector<int> az(G); // use int for index to be compatibable to R
    k = 0, j = 0;
    // binary search for the closest symbol to zg[k] in W[k] keys if not exists
    auto wz = W[k].upper_bound(zg[k]) == W[k].begin() ? W[k].begin() : prev(W[k].upper_bound(zg[k]));
    // for k=0, start with random occurrence of the symbol. and put zg at that position
    az[k] = wz->second.rbegin()->first; // W[k] : {symbol: {index: rank}};
    // print out
    for (const auto& si : W[k])
        cerr << bitset<B>(reverseBits(si.first)) << ", " << si.first << "\n";
    cerr << "\nk:" << k << "\nclosest symbol to zg[" << bitset<B>(reverseBits(zg[k])) << "]:\t" << bitset<B>(reverseBits(wz->first))
         << "\norder in a[k+1]:" << az[k] << endl;
    cerr << "\nk:" << k << "\nclosest symbol to zg[" << (zg[k]) << "]:\t" << (wz->first) << "\norder in a[k]:" << az[k] << endl;
    // for k > 0
    for (k = 1; k < G; k++)
    {
        auto wz = W[k].upper_bound(zg[k]) == W[k].begin() ? *W[k].begin() : *prev(W[k].upper_bound(zg[k]));
        if (az[k - 1] < C[k][wz.first])
        {
            az[k] = C[k][wz.first];
        }
        else
        {
            auto wz = W[k].upper_bound(zg[k]) == W[k].begin() ? W[k].begin() : prev(W[k].upper_bound(zg[k]));
            if (az[k - 1] < C[k][wz->first])
            {
                j++;
                az[k] = C[k][wz->first];
            }
            else
            {
                auto wi = wz->second.upper_bound(az[k - 1]) == wz->second.begin() ? *wz->second.begin() : *prev(wz->second.upper_bound(az[k - 1]));
                az[k] = wi.second + C[k][wz->first];
            }
        }
    }
    for (const auto& ai : az)
        cerr << ai << "\n";
    coutZYk(X, A, zg, az, ki);
}
