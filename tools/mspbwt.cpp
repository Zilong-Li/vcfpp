// -*- compile-command: "g++ mspbwt.cpp -o mspbwt -std=c++17 -g -O3 -Wall -lhts -lz -lm -lbz2 -llzma -lcurl && ./mspbwt -i t.bcf -r chr22 > t" -*-
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
T reverseBits(T n, size_t b = sizeof(T) * 8)
{
    assert(b <= std::numeric_limits<T>::digits);
    T rv = 0;
    for (size_t i = 0; i < b; ++i, n >>= 1) {
        rv = (rv << 1) | (n & 0x01);
    }
    return rv;
}

using IntGridVec = vector<uint32_t>;
using IntSet = set<uint32_t, less<uint32_t>>;
using IntMap = unordered_map<uint32_t, uint32_t>;
using IntVecMap = unordered_map<uint32_t, std::vector<uint32_t>>;
using WgSymbolMap = map<uint32_t, map<uint32_t, uint32_t, less<uint32_t>>, less<uint32_t>>;

vector<uint32_t> encode_z2grid(const vector<bool>& z, int G);
vector<bool> randhapz(uint64_t M);
IntMap build_C(const IntGridVec& x, const IntSet& s);
WgSymbolMap save_W(const IntGridVec& x, const IntSet& s);
IntVecMap build_W(const IntGridVec& x, const IntSet& s, const IntMap& C);
void mspbwt(const std::string& vcffile, const std::string& samples, const std::string& region, int ki);

void coutZYk(const vector<IntGridVec>& X, const vector<vector<int>>& A, const IntGridVec& zg, const vector<int>& ta, int k, bool bit = 1);
void coutWgSymbolMap(const WgSymbolMap& wg);

int main(int argc, char* argv[])
{
    // ========= helper message and parameters parsing ============================
    std::vector<std::string> args(argv + 1, argv + argc);
    if (argc <= 1 || args[0] == "-h" || args[0] == "-help")
    {
        std::cout << "Author: Zilong-Li (zilong.dk@gmail.com)\n"
                  << "Usage example:\n"
                  << "     " + (std::string)argv[0] + " -i in.bcf \n"
                  << "     " + (std::string)argv[0] + " -i in.bcf -o out.bcf -s ^S1,S2 -r chr1:1-1000 \n"
                  << "     bcftools view in.bcf | " + (std::string)argv[0] + " -i - -o out.bcf \n"
                  << "\nOptions:\n"
                  << "     -i    input vcf/bcf file\n"
                  << "     -o    ouput vcf/bcf file [stdout]\n"
                  << "     -s    list of samples to be included or excluded\n"
                  << "     -r    specific region to be included\n"
                  << "     -k    Y[k] show [1]\n"
                  << std::endl;
        return 1;
    }
    std::string invcf, outvcf = "-", samples = "-", region = "";
    int ki{1};
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
        if (args[i] == "-k")
            ki = stoi(args[++i]);
    }

    // ========= core calculation part ===========================================

    mspbwt(invcf, samples, region, ki);

    return 0;
}

vector<uint32_t> encode_z2grid(const vector<bool>& z, int G)
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

void coutZYk(const vector<IntGridVec>& X, const vector<vector<int>>& A, const IntGridVec& zg, const vector<int>& ta, int k, bool bit)
{
    const size_t B = sizeof(X[0][0]) * 4;
    for (size_t i = 0; i < X[0].size(); i++)
    {
        for (int j = 0; j < k + 1; j++)
        {
            auto rb = reverseBits(X[j][A[k][i]]);
            if (bit)
                cout << std::bitset<B>(rb) << " ";
            else
                cout << rb << " ";
        }
        cout << endl;

        if (i == ta[k])
        {
            // print out original Z bits
            cout << "========= zg is inserting here ========  k=" << k << ", ta[k]=" << ta[k] << endl;
            for (int j = 0; j < k + 1; j++)
            {
                auto rb = reverseBits(zg[j]);
                if (bit)
                    cout << std::bitset<B>(rb) << " ";
                else
                    cout << rb << " ";
            }
            cout << endl;
            cout << "========= zg is inserting here ========  k=" << k << ", ta[k]=" << ta[k] << endl;
        }
    }
}


void coutWgSymbolMap(const WgSymbolMap& wg)
{
    auto print_key_value = [](const auto& key, const auto& value) { std::cout << "Index:[" << key << "] rank:[" << value << "]\n"; };
    for (const auto& [s, v] : wg)
    {
        cout << "symbol:[" << s << "]\n";
        for (const auto& [idx, i] : v)
        {
            print_key_value(idx, i);
        }
    }
}

vector<bool> randhapz(uint64_t M)
{
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dist(0, 1);
    vector<bool> z(M);
    for (uint64_t i = 0; i < M; i++)
    {
        z[i] = dist(gen);
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

void mspbwt(const std::string& vcffile, const std::string& samples, const std::string& region, int ki)
{
    const int B = 32;
    BcfReader vcf(vcffile, samples, region);
    BcfRecord var(vcf.header);
    uint64_t N{0}, M{0}, G{0}, k{0}, m{0}, i{0}; // N haplotypes, M SNPs, G Grids, k Current grid, m Current SNP
    N = vcf.nsamples * 2;
    M = vcf.get_region_records(region);
    G = (M + B - 1) / B;
    cerr << "Haps(N):" << N << "\tSNPs(M):" << M << "\tGrids(G):" << G << "\tInt(B):" << B << endl;
    vector<IntGridVec> X; // Grids x Haps
    X.resize(G, IntGridVec(N));
    vcf.setRegion(region); // seek back to region
    vector<bool> gt;
    while (vcf.getNextVariant(var))
    {
        var.getGenotypes(gt);
        // update current grids
        for (i = 0; i < N; i++)
            X[k][i] = (X[k][i] << 1) | (gt[i] != 0);
        m++;
        if (m % B == 0)
        {
            X[k][i] = reverseBits(X[k][i]); // reverset bits
            k++; // update next grid
        }
    }
    if (G == k + 1)
    {
        int pad = G * B - M;
        for (i = 0; i < N; i++)
            X[k][i] <<= pad;
        X[k][i] = reverseBits(X[k][i]); // reverset bits
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
    IntGridVec y0(N);
    vector<vector<int>> A; // (Grids+1) x Haps, we use int here so that R can take it
    A.resize(G + 1, vector<int>(N));
    vector<int> a0(N);
    for (k = 0; k < G; k++)
    {
        if (k == 0)
        {
            std::iota(a0.begin(), a0.end(), 0);
            A[k] = a0;
        }
        for (i = 0; i < N; i++)
            y0[i] = X[k][a0[i]];
        IntSet s(y0.begin(), y0.end()); // convert to set which unique sorted
        C[k] = build_C(y0, s);
        auto Wg = build_W(y0, s, C[k]); // here Wg is S x N
        for (i = 0; i < N; i++)
            A[k + 1][Wg[y0[i]][i] - 1] = a0[i];
        // next run
        a0 = A[k + 1];
        // here save current W, which differs from the previous complete table
        W[k] = save_W(y0, s);
    }

    // finially insert zg back to A at each grid
    auto z = randhapz(M);
    auto zg = encode_z2grid(z, G);
    vector<int> ta(G+1); // use int for index to be compatibable to R
    k = 0;
    ta[k] = N;
    // binary search for the closest symbol to zg[k] in W[k] keys if not exists
    auto wz = W[k].lower_bound(zg[k]) != W[k].end() ? *W[k].lower_bound(zg[k]) : *W[k].crbegin();
    // for k=0, start with random occurrence of the symbol. and put zg at that position
    ta[k+1] = wz.second.crbegin()->first; // W[k] : {symbol: {index: rank}};
    // print out
    // for (const auto& si : W[k])
    //     cerr << bitset<B>(si.first) << ", " << si.first << "\n";
    cerr << "\nk:" << k << "\nclosest symbol to zg[" << bitset<B>(zg[k]) << "]:\t" << bitset<B>(wz.first) << "\norder in a[k+1]:" << ta[k+1] << endl;
    cerr << "\nk:" << k << "\nclosest symbol to zg[" << (zg[k]) << "]:\t" << (wz.first) << "\norder in a[k]:" << ta[k+1] << endl;
    // coutWgSymbolMap(W[k]);
    // for k > 0
    for (k = 1; k < G; k++)
    {
        auto wz = W[k].lower_bound(zg[k]) != W[k].end() ? *W[k].lower_bound(zg[k]) : *W[k].crbegin();
        auto wi = wz.second.lower_bound(ta[k]) != wz.second.end() ? *wz.second.lower_bound(ta[k - 1]) : *wz.second.crbegin();
        ta[k+1] = wi.second + C[k][wz.first];
    }
    for (const auto& ai : ta)
        cerr << ai << "\n";
    coutZYk(X, A, zg, ta, ki);
}
