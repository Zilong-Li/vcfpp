// -*- compile-command: "g++ mspbwt-example.cpp -std=c++17 -g -O3 -Wall -lhts -lz -lm -lbz2 -llzma -lcurl && ./a.out -p panel.vcf.gz -q 10 -r chr22 > t" -*-
#include "../vcfpp.h"

#include <algorithm>
#include <bitset>
#include <chrono>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <unordered_map>


using namespace std;
using namespace vcfpp;

class Timer
{
protected:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_timing_clock, prev_timing_clock;

public:
    Timer();
    ~Timer();
    void clock();
    unsigned int reltime();
    unsigned int abstime();
};

Timer::Timer()
{
    start_timing_clock = std::chrono::high_resolution_clock::now();
}

Timer::~Timer()
{
}

void Timer::clock()
{
    prev_timing_clock = std::chrono::high_resolution_clock::now();
}

unsigned int Timer::reltime()
{
    return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - prev_timing_clock).count();
}

unsigned int Timer::abstime()
{
    return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_timing_clock).count();
}

template <typename T>
T reverseBits(T n, size_t B = sizeof(T) * 8)
{
    assert(B <= std::numeric_limits<T>::digits);
    T rv = 0;
    for (size_t i = 0; i < B; ++i, n >>= 1)
        rv = (rv << 1) | (n & 0x01);
    return rv;
}

using grid_t = uint32_t;
using IntGridVec = vector<grid_t>;
using IntSet = set<grid_t, less<grid_t>>;
using IntMap = unordered_map<grid_t, uint32_t>;                 // {symbol : index}
using IntVecMap = unordered_map<grid_t, std::vector<uint32_t>>; // {symbol : vec(index)}
using SymbolIdxMap = map<uint32_t, uint32_t, less<grid_t>>;     // {index: rank}
using WgSymbolMap = map<grid_t, SymbolIdxMap, less<grid_t>>;    // {symbol:{index:rank}}

void mspbwt(const std::string& vcfpanel, const std::string& samples, const std::string& region, int q, int fast);

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
                  << "     -q    first n numer of samples as query [1]\n"
                  << "     -r    chromosome to be included\n"
                  << "     -f    force fast mode for building indices more ram.[0]\n"
                  << std::endl;
        return 1;
    }
    std::string vcfpanel, vcfquery, outvcf = "-", samples = "-", region = "";
    int q{1}, fast{0};
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
        if (args[i] == "-f")
            fast = stoi(args[++i]);
    }

    // ========= core calculation part ===========================================

    mspbwt(vcfpanel, samples, region, q, fast);

    return 0;
}

IntGridVec encodeZgrid(const vector<bool>& z, int G)
{
    IntGridVec zg(G);
    const int B = sizeof(grid_t) * 8;
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


IntMap build_C(const IntGridVec& x, const IntSet& s)
{
    uint32_t c{0};
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
    uint32_t c{0}, i{0}, n{0}, cap{0};
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
        if (++n == cap && cap > 0)
            break;
    }

    return W;
}

IntMap build_Symbols(const IntSet& s)
{
    uint32_t n{0};
    IntMap symbol;
    for (const auto& si : s)
    {
        symbol[si] = n++;
    }
    return symbol;
}

void mspbwt(const std::string& vcfpanel, const std::string& samples, const std::string& region, int q, int fast)
{
    Timer tm;
    const int B = sizeof(grid_t) * 8;
    BcfReader vcf(vcfpanel, samples, region);
    BcfRecord var(vcf.header);
    uint64_t Nq{0}, Np{0}, M{0}, G{0}, k{0}, m{0}, i{0}, j{0}, w{0}; // N haplotypes, M SNPs, G Grids, k Current grid, m Current SNP
    Nq = q * 2;
    Np = vcf.nsamples * 2 - Nq;
    vector<bool> gt;
    while (vcf.getNextVariant(var))
    {
        var.getGenotypes(gt);
        if (!var.isNoneMissing() || !var.allPhased())
            continue;
        M++;
    }
    G = (M + B - 1) / B;
    cerr << "Haps(Np):" << Np << "\tHaps(Nq):" << Nq << "\tSNPs(M):" << M << "\tGrids(G):" << G << "\tInt(B):" << B << endl;
    vector<IntGridVec> X, Z; // Grids x Haps
    X.resize(G, IntGridVec(Np));
    Z.resize(Nq, IntGridVec(G));
    vcf.setRegion(region); // seek back to region
    while (vcf.getNextVariant(var))
    {
        // check if var has missing of phased values
        var.getGenotypes(gt);
        if (!var.isNoneMissing() || !var.allPhased())
            continue;
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
        if ((symbols.size() < 256) || fast)
        {
            auto Wg = build_W(y0, symbols, C[k]); // here Wg is S x N
            for (i = 0; i < Np; i++)
                A[k + 1][Wg[y0[i]][i] - 1] = a0[i];
        }
        else
        {
            // without Wg = build_W is slow but memory efficient
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
        }
        // next run
        a0 = A[k + 1];
        symbols.clear();
    }
    cerr << "elapsed time of buiding indices: " << tm.abstime() << " seconds." << endl;

    vector<int> az(G); // use int for index to be compatibable to R
    for (const auto& zg : Z)
    {
        tm.clock();
        k = 0, j = 0;
        // binary search for the closest symbol to zg[k] in W[k] keys if not exists
        auto wz = W[k].upper_bound(zg[k]) == W[k].begin() ? W[k].begin() : prev(W[k].upper_bound(zg[k]));
        // for k=0, start with random occurrence of the symbol. and put zg at that position
        az[k] = wz->second.rbegin()->first; // W[k] : {symbol: {index: rank}};
        // for k > 0
        for (k = 1; k < G; k++)
        {
            auto wz = W[k].upper_bound(zg[k]) == W[k].begin() ? W[k].begin() : prev(W[k].upper_bound(zg[k]));
            if (az[k - 1] < C[k][wz->first])
            {
                j++;
                az[k] = C[k][wz->first];
            }
            else
            {
                auto wi = wz->second.upper_bound(az[k - 1]) == wz->second.begin() ? wz->second.begin() : prev(wz->second.upper_bound(az[k - 1]));
                az[k] = wi->second + C[k][wz->first];
            }
        }
        cerr << "elapsed time of query hap z: " << tm.reltime() << " milliseconds. " << j << "/" << G << " grids skipped searching\n";
    }
    cerr << "elapsed time of whole program: " << tm.abstime() << " seconds." << endl;
}
