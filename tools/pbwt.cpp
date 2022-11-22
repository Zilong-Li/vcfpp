// -*- compile-command: "g++ pbwt.cpp -o pbwt -std=c++17 -g -O3 -Wall -lhts -lz -lm -lbz2 -llzma -lcurl && ./pbwt " -*-
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "../vcfpp.h"

using namespace std;
using namespace vcfpp;

typedef uint32_t uint;
typedef uint64_t uint2;

/// build and dump pbwt FM-index to disk;
void pbwt_index(const std::string& vcffile, const std::string& samples, const std::string& region)
{
    BcfReader vcf(vcffile, samples, region);
    BcfRecord var(vcf.header);
    std::vector<bool> gt;
    std::vector<std::vector<bool>> x;
    int nsnps = 0, nsamples = vcf.nsamples;
    while (vcf.getNextVariant(var))
    {
        if (!var.isSNP())
            continue; // skip other type of variants
        var.getGenotypes(gt);
        x.push_back(gt);
        nsnps++;
    }
    // build pbwt index
    int N = nsnps, M = nsamples * 2;
    std::ofstream ofpbwt(vcffile + ".pbwt", std::ios::binary);
    std::ofstream ofauxu(vcffile + ".auxu", std::ios::binary);
    ofpbwt.write(reinterpret_cast<char*>(&N), 4);
    ofpbwt.write(reinterpret_cast<char*>(&M), 4);
    std::vector<uint> at(M), a(M), a0(M), a1(M), u(M + 1);
    int i, j, u_, v_, a_;
    for (j = 0; j < N; j++)
    {
        u_ = 0;
        v_ = 0;
        // copy a to at, which stands for previous a at j - 1
        for (i = 0; i < M; i++)
        {
            at[i] = a[i];
        }
        for (i = 0; i < M; i++)
        {
            a_ = j > 0 ? at[i] : i;
            u[i] = u_;
            if (x[j][a_])
                a1[v_++] = a_;
            else
                a0[u_++] = a_;
        }
        u[M] = u_;
        for (i = 0; i < M; i++)
        {
            if (i < u_)
                a[i] = a0[i];
            else
                a[i] = a1[i - u_];
        }
        if (!ofpbwt.write(reinterpret_cast<char*>(&a[0]), a.size() * 4))
            throw std::runtime_error(strerror(errno));
        if (!ofauxu.write(reinterpret_cast<char*>(&u[0]), u.size() * 4))
            throw std::runtime_error(strerror(errno));
    }
}

std::vector<int> pbwt_query(const std::string& vcffile, const std::vector<bool>& z, int s, int L, int Step)
{
    int N, M, i, j;
    std::ifstream ifpbwt(vcffile + ".pbwt", std::ios::binary);
    if (!ifpbwt.read(reinterpret_cast<char*>(&N), 4))
        throw std::runtime_error(strerror(errno));
    if (!ifpbwt.read(reinterpret_cast<char*>(&M), 4))
        throw std::runtime_error(strerror(errno));
    assert(z.size() == N);
    std::ifstream ifauxu(vcffile + ".auxu", std::ios::binary);
    std::vector<int> t(N), u(M + 1), a(M);
    std::vector<int> matches;
    for (i = 0; i < N; i++)
    {
        ifauxu.read(reinterpret_cast<char*>(&u[0]), 4 * (M + 1));
        if (i == 0)
        {
            t[0] = z[0] ? M : 0;
        }
        else
        {
            if (z[i])
                t[i] = u[M] + t[i - 1] - u[t[i - 1]];
            else
                t[i] = u[t[i - 1]];
        }
        if (s < N && i == s)
        {
            ifpbwt.seekg((uint2)s * M * 4 + 8, std::ios_base::beg);
            ifpbwt.read(reinterpret_cast<char*>(&a[0]), 4 * M);
            if (t[s] == M)
            {
                // where we hit the bottom
                for (j = 1; j <= L; j++)
                    matches.push_back(a[M - j]);
            }
            else if (t[s] == 0)
            {
                // where we hit the top;
                for (j = 0; j < L; j++)
                    matches.push_back(a[0 + j]);
            }
            else
            {
                // L haps before z;
                for (j = 1; j <= L; j++)
                    matches.push_back(a[(t[s] - j) > 0 ? (t[s] - j) : 0]);
                // L haps after z;
                for (j = 0; j < L; j++)
                    matches.push_back(a[(t[s] + j) < M - 1 ? (t[s] + j) : (M - 1)]);
            }
            s += Step;
        }
    }
    return matches;
}


int main(int argc, char* argv[])
{
    std::vector<std::string> args(argv + 1, argv + argc);
    if (argc <= 1 || args[0] == "-h" || args[0] == "-help")
    {
        std::cout << "Author: Zilong-Li (zilong.dk@gmail.com)\n"
                  << "Description:\n"
                  << "     build bit pbwt indicies and query the neighbouring haps in pbwt panel given a new hap.\n\n"
                  << "Usage example:\n"
                  << "     " + (std::string)argv[0] + " index -vcf-panel VCF1 -samples ^S1,S2 -region chr1:1-1000 \n"
                  << "     " + (std::string)argv[0] + " query -vcf-panel VCF1 -vcf-query VCF2 -out OUT \n\n"
                  << "Options:\n"
                  << "     -vcf-panel    vcf/bcf file of reference panel\n"
                  << "     -vcf-query    vcf/bcf file of query sequence\n"
                  << "     -sequence     single query sequence eg. 1010010\n"
                  << "     -samples      list of samples to be included or excluded\n"
                  << "     -region       specific region be included\n"
                  << "     -out          output file\n"
                  << "     -L (=2)       number of neighour haplotypes to be selected\n"
                  << "     -S (=8)       make a selection every S SNPs\n"
                  << "     -seed (=1)    seed for generating random number\n"
                  << "     -help         show help message\n"
                  << std::endl;
        return 1;
    }
    std::string vcfpanel, vcfquery, samples = "-", region = "", outfile, sequence;
    int seed = 1, L = 2, Step = 8;
    for (int i = 1; i < args.size(); i++)
    {
        if (args[i] == "-vcf-panel")
            vcfpanel = args[++i];
        if (args[i] == "-vcf-query")
            vcfquery = args[++i];
        if (args[i] == "-sequence")
            sequence = args[++i];
        if (args[i] == "-samples")
            samples = args[++i];
        if (args[i] == "-region")
            region = args[++i];
        if (args[i] == "-out")
            outfile = args[++i];
        if (args[i] == "-seed")
            seed = std::stoi(args[++i]);
        if (args[i] == "-L")
            L = std::stoi(args[++i]);
        if (args[i] == "-S")
            Step = std::stoi(args[++i]);
    }
    if (args[0] == "index")
    {
        pbwt_index(vcfpanel, samples, region);
    }
    else if (args[0] == "query")
    {
        if (!sequence.empty())
        {
            std::vector<bool> z(sequence.size());
            for (int i = 0; i < sequence.size(); i++)
                z[i] = (sequence[i] == '1');
            auto res = pbwt_query(vcfpanel, z, seed, L, Step);
            // print out results
            for (auto& o : res)
                std::cout << o << std::endl;
        }
        else if (!vcfquery.empty())
        {
            // todo
            std::cout << "implementation missing for -vcf-query\n";
        }
        else
        {
            std::cout << "please use either -vcf-query or -sequence option for query\n";
            return 1;
        }
    }
    else
    {
        std::cout << "please run either index or query. use -h to see usage example\n";
        return 1;
    }

    return 0;
}
