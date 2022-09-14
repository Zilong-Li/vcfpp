#include "vcfpp.h"
#include <fstream>

using namespace std;
using namespace vcfpp;

typedef uint32_t uint;

/// build and dump pbwt FM-index to disk;
void pbwt_index(const std::string& vcffile, const std::string& samples, const std::string& region)
{
    BcfReader vcf(vcffile, samples, region);
    BcfRecord var(vcf.header);
    vector<bool> gt;
    vector<vector<bool>> x;
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
    uint N = nsnps, M = nsamples * 2;
    ofstream ofpbwt(vcffile + ".pbwt", std::ios::binary);
    ofstream ofauxu(vcffile + ".auxu", std::ios::binary);
    ofstream ofauxv(vcffile + ".auxv", std::ios::binary);
    ofpbwt.write(reinterpret_cast<char*>(&N), 4);
    ofpbwt.write(reinterpret_cast<char*>(&M), 4);
    vector<uint> at(M), a(M), a0(M), a1(M), u(M + 1), v(M + 1);
    uint i, j, u_, v_, a_;
    for (j = 0; j < N; j++)
    {
        u_ = 0; v_ = 0;
        // copy a to at
        for (i = 0; i < M; i++) { at[i] = a[i]; }
        for (i = 0; i < M; i++)
        {
            a_ = j > 0 ? at[i] : i;
            u[i] = u_; v[i] = v_;
            if (x[j][a_])
            {
                a1[v_] = a_;
                v_++;
            }
            else
            {
                a0[u_] = a_;
                u_++;
            }
        }
        u[M] = u_; v[M] = M;
        for (i = 0; i < M; i++)
        {
            v[i] += u_;
            if (i < u_)
            {
                a[i] = a0[i];
            }
            else
            {
                a[i] = a1[i - u_];
            }
        }
        if (!ofpbwt.write(reinterpret_cast<char*>(&a[0]), a.size() * 4))
            throw runtime_error(strerror(errno));
        if (!ofauxu.write(reinterpret_cast<char*>(&u[0]), u.size() * 4))
            throw runtime_error(strerror(errno));
        if (!ofauxv.write(reinterpret_cast<char*>(&v[0]), v.size() * 4))
            throw runtime_error(strerror(errno));
    }
}

void pbwt_query(const std::string& vcffile, const std::vector<bool>& z)
{
    uint N, M, i, j;
    ifstream ifpbwt(vcffile + ".pbwt", std::ios::binary);
    ifstream ifauxu(vcffile + ".auxu", std::ios::binary);
    ifstream ifauxv(vcffile + ".auxv", std::ios::binary);
    if (!ifpbwt.read(reinterpret_cast<char*>(&N), 4))
        throw runtime_error(strerror(errno));
    if (!ifpbwt.read(reinterpret_cast<char*>(&M), 4))
        throw runtime_error(strerror(errno));
    vector<uint> t(N), v(M + 1), u(M + 1), a(M);
    vector<uint> matches;
    uint s = 1, Step = 32, L = 2;
    for (i = 0; i < N; i++)
    {
        ifauxu.read(reinterpret_cast<char*>(&u[0]), 4 * (M + 1));
        ifauxv.read(reinterpret_cast<char*>(&v[0]), 4 * (M + 1));
        if (i == 0)
        {
            t[0] = z[0] ? M : 0;
        }
        else
        {
            if (z[i])
                t[i] = v[t[i - 1]];
            else
                t[i] = u[t[i - 1]];
        }
        if (s < N && i == s)
        {
            ifpbwt.seekg((s - 1) * M * 4 + 8, ios_base::beg);
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
                    matches.push_back(a[(t[s] + j) > M - 1] ? (t[s] + j) : (M - 1));
            }
            s += Step;
        }
    }
}

int main(int argc, char* argv[])
{

    return 0;
}
