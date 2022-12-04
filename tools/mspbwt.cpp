// -*- compile-command: "g++ mspbwt.cpp -o mspbwt -std=c++17 -g -O3 -Wall -lhts -lz -lm -lbz2 -llzma -lcurl && ./mspbwt -i t.bcf -r chr22 > t" -*-
#include "../vcfpp.h"

#include <algorithm>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <unordered_map>


using namespace std;
using namespace vcfpp;

using IntMap = unordered_map<int, int>;
using IntVecMap = unordered_map<int, std::vector<int>>;
using WgSymbolMap = map<int, map<int, int>>;

void coutWgSymbolMap(const WgSymbolMap& wg);
vector<int> encode_z2grid(const vector<bool>& z, int G);
vector<bool> randhapz(int M);
int find_closest_symbol2zg(int zg, const std::set<int>& symbols);
IntMap build_C(const vector<int>& x, const set<int>& s);
IntMap save_Symbols(const vector<int>& x);
WgSymbolMap save_W(const vector<int>& x, const set<int>& s);
IntVecMap build_W(const vector<int>& x, const set<int>& s, const IntMap& C);
void mspbwt(const std::string& vcffile, const std::string& samples, const std::string& region);

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

    mspbwt(invcf, samples, region);

    return 0;
}

vector<int> encode_z2grid(const vector<bool>& z, int G)
{
    vector<int> zg(G);
    const int B = 32;
    size_t m{0}, k{0}, M{z.size()};
    for (m = 0; m < M; m++)
    {
        zg[k] = (zg[k] << 1) | (z[m] != 0);
        if ((m + 1) % B == 0)
            k++;
    }
    if (G == k + 1)
    {
        zg[k] <<= G * B - M; // padding 0s
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

vector<bool> randhapz(int M)
{
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dist(0, 1);
    vector<bool> z(M);
    for (int i = 0; i < M; i++)
    {
        z[i] = dist(gen);
    }
    return z;
}

int find_closest_symbol2zg(int zg, const std::set<int>& symbols)
{
    // vec is ordered. so binary search
    // auto lower = std::lower_bound(symbols.begin(), symbols.end(), zg);
    auto lower = symbols.lower_bound(zg);
    if (lower != symbols.end())
        cerr << "no end\n";
    if (lower != symbols.end())
        return *lower;
    else
        return *symbols.rbegin();
}

IntMap save_Symbols(const vector<int>& x)
{
    set<int> s(x.begin(), x.end()); // convert to set which unique sorted
    IntMap symbols;                 // key: symbol; value: index in set
    int n{0};
    for (const auto& si : s)
        symbols[si] = n++;
    return symbols;
}

IntMap build_C(const vector<int>& x, const set<int>& s)
{
    IntMap C;
    int c{0}, n{0}, cap{0};
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

WgSymbolMap save_W(const vector<int>& x, const set<int>& s)
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

IntVecMap build_W(const vector<int>& x, const set<int>& s, const IntMap& C)
{
    int c{0}, i{0}, n{0}, cap{0};
    IntVecMap W;
    for (const auto& si : s)
    {
        W[si] = vector<int>(x.size());
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

void mspbwt(const std::string& vcffile, const std::string& samples, const std::string& region)
{
    const int B = 32;
    BcfReader vcf(vcffile, samples, region);
    BcfRecord var(vcf.header);
    uint64_t N{0}, M{0}, G{0}, k{0}, m{0}, i{0}; // N haplotypes, M SNPs, G Grids, k Current grid, m Current SNP
    N = vcf.nsamples * 2;
    M = vcf.get_region_records(region);
    G = (M + B - 1) / B;
    cerr << "Haps(N):" << N << "\tSNPs(M):" << M << "\tGrids(G):" << G << "\tInt(B):" << B << endl;
    vector<vector<int>> X; // Grids x Haps
    X.resize(G, vector<int>(N));
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
            k++; // update next grid
    }
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
    vector<WgSymbolMap> W(G);
    vector<IntMap> C(G);
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
        set<int> s(y0.begin(), y0.end()); // convert to set which unique sorted
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
    vector<int> ta(G);
    k = 0;
    // binary search for the closest symbol to zg[k] in W[k]->keys if not exists
    auto wz = W[k].lower_bound(zg[k]) != W[k].end() ? *W[k].lower_bound(zg[k]) : *W[k].crbegin();
    for (const auto& si : W[k])
        cerr << si.first << ",";
    cerr << "\nzg:" << zg[k] << "\t" << wz.first << endl;
    // coutWgSymbolMap(W[k]);
    // for k=0, start with random occurrence of the symbol. and put zg at that position
    ta[k] = wz.second.cbegin()->first;
    cerr << "k:" << k << "\tclosest symbol to zg:" << wz.first << "\torder in a[k]:" << ta[k] << endl;
    // for k > 0
    for (k = 1; k < G; k++)
    {
        auto wz = W[k].lower_bound(zg[k]) != W[k].end() ? *W[k].lower_bound(zg[k]) : *W[k].crbegin();
        auto wi = wz.second.lower_bound(ta[k - 1]) != wz.second.end() ? *wz.second.lower_bound(ta[k - 1]) : *wz.second.crbegin();
        ta[k] = wi.second + C[k][wz.first];
    }
}
