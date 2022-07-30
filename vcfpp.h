/* File: vcfpp.h
** Author: Zilong Li (zilong.dk@gmail.com)
** Copyright (C) 2022
*/

#ifndef VCFPP_H_
#define VCFPP_H_

#include <memory>
#include <string>
#include <type_traits>
#include <vector>

extern "C"
{
#include "htslib/tbx.h"
#include "htslib/vcf.h"
}

namespace vcfpp
{
    bool
    isEndWith(std::string const& s, std::string const& e)
    {
        if (s.length() >= e.length())
        {
            return (0 == s.compare(s.length() - e.length(), e.length(), e));
        }
        else
        {
            return false;
        }
    };

    class BcfHeader
    {
        friend class BcfRecord;
        friend class BcfReader;
        friend class BcfWriter;

    public:
        BcfHeader()
        {
        }

        ~BcfHeader()
        {
            // bcf_hdr_destroy(hdr); // cause double free issue
            bcf_hrec_destroy(hrec);
        }

        // todo : check if the value is valid for vcf specification
        inline void
        AddInfo(const std::string& id, const std::string& number, const std::string& type,
                const std::string& description)
        {
            AddLine("##INFO=<ID=" + id + ",Number=" + number + ",Type=" + type + ",Description=\"" + description +
                    "\">");
        }

        inline void
        AddFormat(const std::string& id, const std::string& number, const std::string& type,
                  const std::string& description)
        {
            AddLine("##FORMAT=<ID=" + id + ",Number=" + number + ",Type=" + type + ",Description=\"" + description +
                    "\">");
        }

        inline void
        AddFilter(const std::string& id, const std::string& number, const std::string& type,
                  const std::string& description)
        {
            AddLine("##FILTER=<ID=" + id + ",Number=" + number + ",Type=" + type + ",Description=\"" + description +
                    "\">");
        }

        inline void
        AddContig(const std::string& id)
        {
            AddLine("##contig=<ID=" + id + ">");
        }

        inline void
        AddLine(const std::string& str)
        {
            ret = bcf_hdr_append(hdr, str.c_str());
            if (ret != 0)
                throw std::runtime_error("could not add " + str + " to header\n");
            ret = bcf_hdr_sync(hdr);
            if (ret != 0)
                throw std::runtime_error("could not add " + str + " to header\n");
        }

        inline void
        AddSample(const std::string& sample)
        {
            bcf_hdr_add_sample(hdr, sample.c_str());
            if (bcf_hdr_sync(hdr) != 0)
            {
                throw std::runtime_error("couldn't add the sample.\n");
            }
        }

        inline std::string
        AsString() const
        {
            kstring_t s = {0, 0, NULL};          // kstring
            if (bcf_hdr_format(hdr, 0, &s) == 0) // append header string to s.s! append!
                return std::string(s.s, s.l);
            else
                throw std::runtime_error("failed to convert formatted header to string");
        }

        std::vector<std::string>
        GetSamples() const
        {
            std::vector<std::string> vec;
            for (int i = 0; i < bcf_hdr_nsamples(hdr); i++)
            {
                vec.emplace_back(std::string(hdr->samples[i]));
            }
            return vec;
        }

        std::vector<std::string>
        GetSeqnames()
        {
            const char** seqs = bcf_hdr_seqnames(hdr, &ret);
            if (ret == 0)
                printf("there is no contig id in the header!\n");
            std::vector<std::string> vec;
            for (int i = 0; i < ret; i++)
            {
                vec.emplace_back(std::string(seqs[i]));
            }
            // TODO: return uninitialized vec may be undefined.
            return vec;
        }

        inline void
        RemoveContig(std::string tag) const
        {

            bcf_hdr_remove(hdr, BCF_HL_CTG, tag.c_str());
        }

        inline void
        RemoveInfo(std::string tag) const
        {
            bcf_hdr_remove(hdr, BCF_HL_INFO, tag.c_str());
        }

        inline void
        RemoveFormat(std::string tag) const
        {
            bcf_hdr_remove(hdr, BCF_HL_FMT, tag.c_str());
        }

        inline void
        RemoveFilter(std::string tag) const
        {
            bcf_hdr_remove(hdr, BCF_HL_FLT, tag.c_str());
        }

        inline void
        SetSamples(const std::string& samples)
        {

            ret = bcf_hdr_set_samples(hdr, samples.c_str(), 0);
            if (ret > 0)
            {
                throw std::runtime_error("the " + std::to_string(ret) + "-th sample are not in the VCF.\n");
            }
            else if (ret == -1)
            {
                throw std::runtime_error("couldn't set samples. something wrong.\n");
            }
        }

        inline void
        SetVersion(const std::string& version) const
        {
            if (bcf_hdr_set_version(hdr, version.c_str()) != 0)
            {
                throw std::runtime_error("couldn't set vcf version correctly.\n");
            }
        }

        int nsamples = 0;

    private:
        bcf_hdr_t* hdr = NULL;   // bcf header
        bcf_hrec_t* hrec = NULL; // populate header
        int ret = 0;
    };

    class BcfRecord
    {
        friend class BcfReader;
        friend class BcfWriter;

    public:
        BcfRecord(BcfHeader& h_) : header(std::make_shared<BcfHeader>(h_))
        {
        }

        ~BcfRecord()
        {
        }

        std::string
        AsString()
        {
            s.s = NULL;
            s.l = 0;
            s.m = 0;
            if (vcf_format(header->hdr, line, &s) == 0)
                return std::string(s.s, s.l);
            else
                throw std::runtime_error("couldn't format current record into a string.\n");
        }

        template <class T>
        typename std::enable_if<std::is_same<T, std::vector<bool>>::value ||
                                    std::is_same<T, std::vector<char>>::value ||
                                    std::is_same<T, std::vector<int>>::value,
                                bool>::type
        GetGenotypes(T& gv)
        {
            ndst = 0;
            ret = bcf_get_genotypes(header->hdr, line, &gts, &ndst);
            if (ret <= 0)
                return 0; // gt not present
            gv.resize(ret);
            nploidy = ret / header->nsamples;
            int i, j, k = 0, nphased = 0;
            // gts = static_cast<int32_t*>(buf);
            for (i = 0; i < header->nsamples; i++)
            {
                for (j = 0; j < nploidy; j++)
                {
                    gv[k++] = bcf_gt_allele(gts[j + i * nploidy]) != 0; // only parse 0 and 1, ie max(nploidy)=2; other
                                                                        // values 2,3... will be converted to 1;
                }
                nphased += (gts[1 + i * nploidy] & 1) == 1;
            }
            if (nphased == header->nsamples)
                isAllPhased = true;
            return 1;
        }

        // return a array for the requested field
        // template <typename T, typename = typename std::enable_if<std::is_same<T,
        // char>::value || std::is_same<T, int>::value || std::is_same<T,
        // float>::value>::type>
        template <typename T, typename S = typename T::value_type>
        typename std::enable_if<std::is_same<T, std::vector<char>>::value || std::is_same<T, std::vector<int>>::value ||
                                    std::is_same<T, std::vector<float>>::value,
                                void>::type
        GetFormat(std::string tag, T& v)
        {
            fmt = bcf_get_fmt(header->hdr, line, tag.c_str());
            shape1 = fmt->n;
            ndst = 0;
            S* buf = NULL;
            if (std::is_same<T, std::vector<int>>::value)
            {
                ret = bcf_get_format_int32(header->hdr, line, tag.c_str(), &buf, &ndst);
            }
            else if (std::is_same<T, std::vector<char>>::value)
            {
                ret = bcf_get_format_char(header->hdr, line, tag.c_str(), &buf, &ndst);
            }
            else if (std::is_same<T, std::vector<float>>::value)
            {
                ret = bcf_get_format_float(header->hdr, line, tag.c_str(), &buf, &ndst);
            }
            if (ret >= 0)
            {
                // user have to check if there is missing in the return v;
                v = std::vector<S>(buf, buf + ret);
            }
            else
            {
                throw std::runtime_error("couldn't parse the " + tag + " format of this variant.\n");
            }
        }

        template <typename T, typename S = typename T::value_type>
        typename std::enable_if<std::is_same<T, std::vector<char>>::value || std::is_same<T, std::vector<int>>::value ||
                                    std::is_same<T, std::vector<float>>::value,
                                void>::type
        GetInfo(std::string tag, T& v)
        {
            info = bcf_get_info(header->hdr, line, tag.c_str());
            S* buf = NULL;
            ndst = 0;
            ret = -1;
            if (info->len > 1)
            {
                if (info->type == BCF_BT_INT8 || info->type == BCF_BT_INT16 || info->type == BCF_BT_INT32)
                {
                    ret = bcf_get_info_int32(header->hdr, line, tag.c_str(), &buf, &ndst);
                }
                else if (info->type == BCF_BT_FLOAT)
                {
                    ret = bcf_get_info_float(header->hdr, line, tag.c_str(), &buf, &ndst);
                }
            }
            if (ret >= 0)
            {
                // user have to check if there is missing in the return v;
                v = std::vector<S>(buf, buf + ret);
            }
            else
            {
                throw std::runtime_error("couldn't parse the " + tag + " format of this variant.\n");
            }
        }

        template <typename T>
        typename std::enable_if<
            std::is_same<T, int>::value || std::is_same<T, float>::value || std::is_same<T, double>::value, void>::type
        GetInfo(std::string tag, T& v)
        {
            info = bcf_get_info(header->hdr, line, tag.c_str());
            if (info->len == 1)
            {
                // scalar value
                if (info->type == BCF_BT_INT8 || info->type == BCF_BT_INT16 || info->type == BCF_BT_INT32)
                {
                    v = info->v1.i;
                }
                else if (info->type == BCF_BT_FLOAT)
                {
                    v = info->v1.f;
                }
            }
        }

        template <typename T>
        typename std::enable_if<std::is_same<T, std::string>::value, void>::type
        GetInfo(std::string tag, T& v)
        {
            info = bcf_get_info(header->hdr, line, tag.c_str());
            if (info->type == BCF_BT_CHAR)
            {
                v = std::string(reinterpret_cast<char*>(info->vptr), info->vptr_len);
            }
        }

        template <typename T>
        typename std::enable_if<std::is_same<T, int>::value || std::is_same<T, float>::value, void>::type
        SetInfo(std::string tag, const T& v)
        {
            ret = -1;
            // bcf_update_info_flag(header->hdr, line, tag.c_str(), NULL, 1);
            int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
            if (std::is_same<T, int>::value)
            {
                if (bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_INT & 0xff))
                    ret = bcf_update_info_int32(header->hdr, line, tag.c_str(), &v, 1);
                else
                    throw std::runtime_error("the given type of tag " + tag + " doesn't match the header");
            }
            else if (std::is_same<T, float>::value)
            {
                // bcf_hrec_set_val
                if (bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_REAL & 0xff))
                    ret = bcf_update_info_float(header->hdr, line, tag.c_str(), &v, 1);
                else
                    throw std::runtime_error("the given type of tag " + tag + " doesn't match the header");
            }
            if (ret < 0)
            {
                throw std::runtime_error("couldn't set " + tag +
                                         " for this variant.\nplease add the tag in header first.\n");
            }
        }

        template <typename T>
        typename std::enable_if<std::is_same<T, std::string>::value || std::is_same<T, std::vector<int>>::value ||
                                    std::is_same<T, std::vector<float>>::value,
                                void>::type
        SetInfo(std::string tag, const T& v)
        {
            ret = -1;
            // bcf_update_info_flag(header->hdr, line, tag.c_str(), NULL, 1);
            int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
            if (std::is_same<T, std::vector<int>>::value)
            {
                if (bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_INT & 0xff))
                    ret = bcf_update_info_int32(header->hdr, line, tag.c_str(), v.data(), v.size());
                else
                    throw std::runtime_error("the given type of tag " + tag + " doesn't match the header");
            }
            else if (std::is_same<T, std::vector<float>>::value)
            {
                if (bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_REAL & 0xff))
                    ret = bcf_update_info_float(header->hdr, line, tag.c_str(), v.data(), v.size());
                else
                    throw std::runtime_error("the given type of tag " + tag + " doesn't match the header");
            }
            else if (std::is_same<T, std::string>::value)
            {
                if (bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_STR & 0xff))
                    ret = bcf_update_info_string(header->hdr, line, tag.c_str(), v.data());
                else
                    throw std::runtime_error("the given type of tag " + tag + " doesn't match the header");
            }
            if (ret < 0)
            {
                throw std::runtime_error("couldn't remove " + tag + " for this variant.\n");
            }
        }

        void
        RemoveInfo(std::string tag)
        {
            ret = -1;
            int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
            if (bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_INT & 0xff))
            {
                ret = bcf_update_info_int32(header->hdr, line, tag.c_str(), NULL, 0);
            }
            else if (bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_REAL & 0xff))
            {
                ret = bcf_update_info_float(header->hdr, line, tag.c_str(), NULL, 0);
            }
            else if (bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_STR & 0xff))
            {
                ret = bcf_update_info_string(header->hdr, line, tag.c_str(), NULL);
            }
            if (ret < 0)
            {
                throw std::runtime_error("couldn't remove " + tag + " for this variant.\n");
            }
        }

        template <class T>
        typename std::enable_if<std::is_same<T, std::vector<bool>>::value ||
                                    std::is_same<T, std::vector<char>>::value ||
                                    std::is_same<T, std::vector<int>>::value,
                                bool>::type
        SetGenotypes(const T& gv, bool phased = false)
        {
            ndst = 0;
            ret = bcf_get_genotypes(header->hdr, line, &gts, &ndst);
            if (ret <= 0)
                return false; // gt not present
            assert(ret == gv.size());
            nploidy = ret / header->nsamples;
            int i, j, k;
            for (i = 0; i < header->nsamples; i++)
            {
                for (j = 0; j < nploidy; j++)
                {
                    k = i * nploidy + j;
                    if (phased)
                        gts[k] = bcf_gt_phased(gv[k]);
                    else
                        gts[k] = bcf_gt_unphased(gv[k]);
                }
            }
            if (bcf_update_genotypes(header->hdr, line, gts, ret) < 0)
                throw std::runtime_error("couldn't set genotypes correctly.\n");
            else
                return true;
        }

        template <typename T>
        typename std::enable_if<std::is_same<T, int32_t>::value || std::is_same<T, float>::value ||
                                    std::is_same<T, double>::value,
                                void>::type
        SetFormat(std::string tag, T v)
        {
            ret = -1;
            int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
            if (std::is_same<T, int32_t>::value)
            {
                if (bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_INT & 0xff))
                    ret = bcf_update_format_int32(header->hdr, line, tag.c_str(), &v, 1);
                else
                    throw std::runtime_error("the given type of tag " + tag + " doesn't match the header");
            }
            else if (std::is_same<T, double>::value)
            {
                float v2 = v;
                if (bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_REAL & 0xff))
                    ret = bcf_update_format_float(header->hdr, line, tag.c_str(), &v2, 1);
                else
                    throw std::runtime_error("the given type of tag " + tag + " doesn't match the header");
            }
            if (ret < 0)
                throw std::runtime_error("couldn't set format " + tag + " correctly.\n");
        }

        template <typename T>
        typename std::enable_if<std::is_same<T, std::string>::value || std::is_same<T, std::vector<char>>::value ||
                                    std::is_same<T, std::vector<int>>::value ||
                                    std::is_same<T, std::vector<float>>::value,
                                void>::type
        SetFormat(std::string tag, const T& v)
        {
            ret = -1;
            int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
            if (std::is_same<T, std::vector<int32_t>>::value)
            {
                if (bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_INT & 0xff))
                    ret = bcf_update_format_int32(header->hdr, line, tag.c_str(), v.data(), v.size());
                else
                    throw std::runtime_error("the given type of tag " + tag + " doesn't match the header");
            }
            else if (std::is_same<T, std::vector<char>>::value || std::is_same<T, std::string>::value)
            {
                if (bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_STR & 0xff))
                    ret = bcf_update_format_char(header->hdr, line, tag.c_str(), v.data(), v.size());
                else
                    throw std::runtime_error("the given type of tag " + tag + " doesn't match the header");
            }
            else if (std::is_same<T, std::vector<float>>::value)
            {
                if (bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_REAL & 0xff))
                    ret = bcf_update_format_float(header->hdr, line, tag.c_str(), v.data(), v.size());
                else
                    throw std::runtime_error("the given type of tag " + tag + " doesn't match the header");
            }
            if (ret < 0)
                throw std::runtime_error("couldn't set format " + tag + " correctly.\n");
        }

        void
        AddLineFromString(const std::string& vcfline)
        {
            std::vector<char> str(vcfline.begin(), vcfline.end());
            str.push_back('\0'); // don't forget string has no \0;
            s.s = &str[0];
            s.l = vcfline.length();
            s.m = vcfline.length();
            ret = vcf_parse(&s, header->hdr, line);
            if (ret > 0)
                throw std::runtime_error("error parsing: " + vcfline + "\n");
            if (line->errcode == BCF_ERR_CTG_UNDEF)
            {
                std::string contig(bcf_hdr_id2name(header->hdr, line->rid));
                hdr_d = bcf_hdr_dup(header->hdr);
                header->hrec = bcf_hdr_id2hrec(hdr_d, BCF_DT_CTG, 0, line->rid);
                if (header->hrec == NULL)
                    throw std::runtime_error("contig" + contig + " unknow and not found in the header.\n");
                ret = bcf_hdr_add_hrec(header->hdr, header->hrec);
                if (ret < 0)
                    throw std::runtime_error("error adding contig " + contig + " to header.\n");
                ret = bcf_hdr_sync(header->hdr);
            }
        }

        bool isAllPhased = false;
        int nploidy = 0;
        int shape1 = 0;
        std::shared_ptr<BcfHeader> header;

    private:
        bcf1_t* line = bcf_init(); // current bcf record
        bcf_hdr_t* hdr_d;          // a dup header by bcf_hdr_dup(header->hdr)
        bcf_fmt_t* fmt = NULL;
        bcf_info_t* info = NULL;
        int32_t* gts = NULL;
        int ndst, ret;
        kstring_t s = {0, 0, NULL}; // kstring
    };

    class BcfReader
    {
    public:
        BcfReader(const std::string& fname_) : fname(fname_), header(BcfHeader())
        {
            fp = hts_open(fname.c_str(), "r");
            header.hdr = bcf_hdr_read(fp);
            header.nsamples = bcf_hdr_nsamples(header.hdr);
        }

        /**
         *  @param samples  samples to include or exclude as a comma-separated string.
         *              LIST        .. select samples in list
         *              ^LIST       .. exclude samples from list
         *              -           .. include all samples
         *              NULL        .. exclude all samples
         */
        BcfReader(const std::string& fname_, const std::string& samples) : fname(fname_)
        {
            fp = hts_open(fname.c_str(), "r");
            header.hdr = bcf_hdr_read(fp);
            header.SetSamples(samples);
            header.nsamples = bcf_hdr_nsamples(header.hdr);
        }

        BcfReader(const std::string& fname_, const std::string& samples, const std::string& region) : fname(fname_)
        {
            fp = hts_open(fname.c_str(), "r");
            header.hdr = bcf_hdr_read(fp);
            header.SetSamples(samples);
            header.nsamples = bcf_hdr_nsamples(header.hdr);
            SetRegion(region);
        }


        ~BcfReader()
        {
            if (fp)
                hts_close(fp);
            if (itr)
                hts_itr_destroy(itr);
        }

        inline int
        SetThreads(int n)
        {
            return hts_set_threads(fp, n);
        }

        // 1. check and load index first
        // 2. query iterval region
        void
        SetRegion(const std::string& region)
        {
            if (isEndWith(fname, "bcf") || isEndWith(fname, "bcf.gz"))
            {
                isBcf = true;
                hidx = bcf_index_load(fname.c_str());
                itr = bcf_itr_querys(hidx, header.hdr, region.c_str());
            }
            else
            {
                isBcf = false;
                tidx = tbx_index_load(fname.c_str());
                assert(tidx != NULL && "error loading tabix index!");
                itr = tbx_itr_querys(tidx, region.c_str());
                assert(itr != NULL && "no interval region found.failed!");
            }
        }

        // TODO reset current region. seek to the first record.
        void Reset();

        bool
        GetNextVariant(BcfRecord& r)
        {
            int ret = -1;
            if (itr != NULL)
            {
                if (isBcf)
                {
                    ret = bcf_itr_next(fp, itr, r.line);
                    return (ret >= 0);
                }
                else
                {
                    int slen = tbx_itr_next(fp, tidx, itr, &s);
                    if (slen > 0)
                    {
                        ret = vcf_parse(&s, r.header->hdr, r.line); // ret > 0, error
                    }
                    return (ret <= 0) && (slen > 0);
                }
            }
            else
            {
                ret = bcf_read(fp, r.header->hdr, r.line);
                return (ret == 0);
            }
        }

        std::string fname;
        bool isBcf;       // if the input file is bcf or vcf;
        BcfHeader header; // bcf header
    private:
        htsFile* fp = NULL;         // hts file
        hts_idx_t* hidx = NULL;     // hts index file
        tbx_t* tidx = NULL;         // .tbi .csi index file for vcf files
        hts_itr_t* itr = NULL;      // hts records iterator
        kstring_t s = {0, 0, NULL}; // kstring
    };

    class BcfWriter
    {
    public:
        BcfWriter(const std::string& fname_) : fname(fname_)
        {
            std::string mode{"w"};
            if (isEndWith(fname, "bcf.gz"))
                mode += "b";
            if (isEndWith(fname, "bcf"))
                mode += "bu";
            if (isEndWith(fname, "vcf.gz"))
                mode += "z";
            fp = hts_open(fname.c_str(), mode.c_str());
        }

        /*!
          @abstract       Open a sequence data (SAM/BAM/CRAM) or variant data
          (VCF/BCF) or possibly-compressed textual line-orientated file
          @param fn       The file name or "-" for stdin/stdout. For indexed files
                          with a non-standard naming, the file name can include the
                          name of the index file delimited with HTS_IDX_DELIM
          @param mode     Mode matching / [rwa][bcefFguxz0-9]* /
          @example
              [rw]b  .. compressed BCF, BAM, FAI
              [rw]bu .. uncompressed BCF
              [rw]z  .. compressed VCF
              [rw]   .. uncompressed VCF
        */

        BcfWriter(const std::string& fname_, const std::string& mode) : fname(fname_)
        {
            fp = hts_open(fname.c_str(), mode.c_str());
        }

        ~BcfWriter()
        {
            hts_close(fp);
            bcf_destroy(b);
        }

        void
        InitalHeader(std::string version = "VCF4.1")
        {
            header = BcfHeader();
            header.hdr = bcf_hdr_init("w");
            header.SetVersion(version);
        }

        // make a copy of given header
        void
        InitalHeader(const BcfHeader& h)
        {
            header = BcfHeader();
            header.hdr = bcf_hdr_dup(h.hdr); // make a copy of given header
            header.nsamples = bcf_hdr_nsamples(header.hdr);
            if (header.hdr == NULL)
                throw std::runtime_error("couldn't copy the header from another vcf.\n");
        }

        void
        WriteLine(const std::string& vcfline)
        {
            if (!isHeaderWritten)
                WriteHeader();
            std::vector<char> line(vcfline.begin(), vcfline.end());
            line.push_back('\0'); // don't forget string has no \0;
            s.s = &line[0];
            s.l = vcfline.length();
            s.m = vcfline.length();
            ret = vcf_parse(&s, header.hdr, b);
            if (ret > 0)
                throw std::runtime_error("error parsing: " + vcfline + "\n");
            if (b->errcode == BCF_ERR_CTG_UNDEF)
            {
                throw std::runtime_error("contig id " + (std::string)bcf_hdr_id2name(header.hdr, b->rid) +
                                         " not found in the header. please run header->AddContig() first.\n");
            }
            ret = bcf_write(fp, header.hdr, b);
            if (ret != 0)
                throw std::runtime_error("error writing: " + vcfline + "\n");
        }

        bool
        WriteHeader()
        {
            ret = bcf_hdr_write(fp, header.hdr);
            if (ret == 0)
                return isHeaderWritten = true;
            else
                return false;
        }

        inline bool
        WriteRecord(BcfRecord& v)
        {
            if (!isHeaderWritten)
                WriteHeader();
            if (bcf_write(fp, v.header->hdr, v.line) < 0)
                return false;
            else
                return true;
        }

        BcfHeader header; // bcf header

    private:
        htsFile* fp = NULL; // hts file
        int ret;
        bcf1_t* b = bcf_init();
        kstring_t s = {0, 0, NULL}; // kstring
        std::string fname;
        bool isHeaderWritten = false;
    };

} // namespace vcfpp

#endif // VCFPP_H_
