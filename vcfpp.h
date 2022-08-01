/*******************************************************************************
 * @file        vcfpp.h
 * @author      Zilong Li
 * @email       zilong.dk@gmail.com
 * @version     v0.0.1
 * @breif       a single C++ file for manipulating VCF
 * @license     MIT
 * Copyright (C) 2022
 * The use of this code is governed by the LICENSE file.
 * @mainpage    something on the mainpage
 ******************************************************************************/

#ifndef VCFPP_H_
#define VCFPP_H_

#include <memory>
#include <string>
#include <type_traits>
#include <vector>

// make sure you have htslib installed
extern "C"
{
#include <htslib/tbx.h>
#include <htslib/vcf.h>
}

namespace vcfpp
{
    bool isEndWith(std::string const& s, std::string const& e)
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

    /**
     * @class BcfHeader
     * @brief Object represents the header in VCF
     * @note  nothing important
     **/
    class BcfHeader
    {
        friend class BcfRecord;
        friend class BcfReader;
        friend class BcfWriter;

    public:
        BcfHeader()
        {
        }

        virtual ~BcfHeader()
        {
        }

        // TODO: check if the value is valid for vcf specification
        inline void addInfo(const std::string& id, const std::string& number, const std::string& type,
                            const std::string& description)
        {
            addLine("##INFO=<ID=" + id + ",Number=" + number + ",Type=" + type + ",Description=\"" + description +
                    "\">");
        }

        inline void addFormat(const std::string& id, const std::string& number, const std::string& type,
                              const std::string& description)
        {
            addLine("##FORMAT=<ID=" + id + ",Number=" + number + ",Type=" + type + ",Description=\"" + description +
                    "\">");
        }

        inline void addFilter(const std::string& id, const std::string& number, const std::string& type,
                              const std::string& description)
        {
            addLine("##FILTER=<ID=" + id + ",Number=" + number + ",Type=" + type + ",Description=\"" + description +
                    "\">");
        }

        inline void addContig(const std::string& id)
        {
            addLine("##contig=<ID=" + id + ">");
        }

        inline void addLine(const std::string& str)
        {
            int ret = 0;
            ret = bcf_hdr_append(hdr, str.c_str());
            if (ret != 0)
                throw std::runtime_error("could not add " + str + " to header\n");
            ret = bcf_hdr_sync(hdr);
            if (ret != 0)
                throw std::runtime_error("could not add " + str + " to header\n");
        }

        inline void addSample(const std::string& sample) const
        {
            bcf_hdr_add_sample(hdr, sample.c_str());
            if (bcf_hdr_sync(hdr) != 0)
            {
                throw std::runtime_error("couldn't add the sample.\n");
            }
        }

        inline std::string asString() const
        {
            kstring_t s = {0, 0, NULL};          // kstring
            if (bcf_hdr_format(hdr, 0, &s) == 0) // append header string to s.s! append!
                return std::string(s.s, s.l);
            else
                throw std::runtime_error("failed to convert formatted header to string");
        }

        std::vector<std::string> getSamples() const
        {
            std::vector<std::string> vec;
            for (int i = 0; i < bcf_hdr_nsamples(hdr); i++)
            {
                vec.push_back(std::string(hdr->samples[i]));
            }
            return vec;
        }

        std::vector<std::string> getSeqnames() const
        {
            int ret = 0;
            const char** seqs = bcf_hdr_seqnames(hdr, &ret);
            if (ret == 0)
                printf("there is no contig id in the header!\n");
            std::vector<std::string> vec;
            for (int i = 0; i < ret; i++)
            {
                vec.push_back(std::string(seqs[i]));
            }
            // TODO: return uninitialized vec may be undefined.
            return vec;
        }

        inline void removeContig(std::string tag) const
        {

            bcf_hdr_remove(hdr, BCF_HL_CTG, tag.c_str());
        }

        inline void removeInfo(std::string tag) const
        {
            bcf_hdr_remove(hdr, BCF_HL_INFO, tag.c_str());
        }

        inline void removeFormat(std::string tag) const
        {
            bcf_hdr_remove(hdr, BCF_HL_FMT, tag.c_str());
        }

        inline void removeFilter(std::string tag) const
        {
            bcf_hdr_remove(hdr, BCF_HL_FLT, tag.c_str());
        }

        inline void setSamples(const std::string& samples) const
        {
            int ret = 0;
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

        inline void setVersion(const std::string& version) const
        {
            bcf_hdr_set_version(hdr, version.c_str());
        }

        inline int nSamples()
        {
            nsamples = bcf_hdr_nsamples(hdr);
            return nsamples;
        }

        int nsamples = 0;

    private:
        bcf_hdr_t* hdr = NULL;   // bcf header
        bcf_hrec_t* hrec = NULL; // populate header
    };

    /**
     * @class BcfRecord
     * @brief Object represents a record in VCF
     * @note  nothing important
     **/
    class BcfRecord
    {
        friend class BcfReader;
        friend class BcfWriter;

    public:
        BcfRecord(BcfHeader& h_) : header(std::make_shared<BcfHeader>(h_))
        {
        }

        virtual ~BcfRecord()
        {
        }

        // TODO resetHeader()
        void resetHeader();

        std::string asString()
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
        getGenotypes(T& gv)
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
                { // only parse 0 and 1, ie max(nploidy)=2; other values 2,3... will be converted to 1;
                    gv[k++] = bcf_gt_allele(gts[j + i * nploidy]) != 0;
                }
                nphased += (gts[1 + i * nploidy] & 1) == 1;
            }
            if (nphased == header->nsamples)
                isAllPhased = true;
            return 1;
        }

        // return a array for the requested field
        template <typename T, typename S = typename T::value_type>
        typename std::enable_if<std::is_same<T, std::vector<char>>::value || std::is_same<T, std::vector<int>>::value ||
                                    std::is_same<T, std::vector<float>>::value,
                                void>::type
        getFormat(std::string tag, T& v)
        {
            fmt = bcf_get_fmt(header->hdr, line, tag.c_str());
            shape1 = fmt->n;
            ndst = 0;
            S* dst = NULL;
            int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
            if (bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_INT & 0xff))
                ret = bcf_get_format_int32(header->hdr, line, tag.c_str(), &dst, &ndst);
            else if (bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_REAL & 0xff))
                ret = bcf_get_format_float(header->hdr, line, tag.c_str(), &dst, &ndst);
            else if (bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_STR & 0xff))
                ret = bcf_get_format_char(header->hdr, line, tag.c_str(), &dst, &ndst);
            if (ret >= 0)
            {
                // user have to check if there is missing in the return v;
                v = std::vector<S>(dst, dst + ret);
            }
            else
            {
                throw std::runtime_error("couldn't parse the " + tag + " format of this variant.\n");
            }
        }

        template <typename T, typename S = typename T::value_type>
        typename std::enable_if<std::is_same<T, std::vector<int>>::value || std::is_same<T, std::vector<float>>::value,
                                void>::type
        getInfo(std::string tag, T& v)
        {
            info = bcf_get_info(header->hdr, line, tag.c_str());
            S* dst = NULL;
            ndst = 0;
            ret = -1;
            if (info->type == BCF_BT_INT8 || info->type == BCF_BT_INT16 || info->type == BCF_BT_INT32)
            {
                ret = bcf_get_info_int32(header->hdr, line, tag.c_str(), &dst, &ndst);
            }
            else if (info->type == BCF_BT_FLOAT)
            {
                ret = bcf_get_info_float(header->hdr, line, tag.c_str(), &dst, &ndst);
            }
            if (ret >= 0)
                v = std::vector<S>(dst, dst + ret); // user have to check if there is missing in the return v;
            else
                throw std::runtime_error("couldn't parse the " + tag + " format of this variant.\n");
        }

        template <typename T>
        typename std::enable_if<
            std::is_same<T, int>::value || std::is_same<T, float>::value || std::is_same<T, double>::value, void>::type
        getInfo(std::string tag, T& v)
        {
            info = bcf_get_info(header->hdr, line, tag.c_str());
            // scalar value
            if (info->len == 1)
            {
                if (info->type == BCF_BT_INT8 || info->type == BCF_BT_INT16 || info->type == BCF_BT_INT32)
                {
                    v = info->v1.i;
                }
                else if (info->type == BCF_BT_FLOAT)
                {
                    v = info->v1.f;
                }
            }
            else
            {
                throw std::runtime_error(tag + " has multiple values. please feed a vector instead.\n");
            }
        }

        template <typename T>
        typename std::enable_if<std::is_same<T, std::string>::value, void>::type getInfo(std::string tag, T& v)
        {
            info = bcf_get_info(header->hdr, line, tag.c_str());
            if (info->type == BCF_BT_CHAR)
                v = std::string(reinterpret_cast<char*>(info->vptr), info->vptr_len);
            else
                throw std::runtime_error(tag + " has to be of string type\n");
        }

        template <typename T>
        typename std::enable_if<std::is_same<T, int>::value || std::is_same<T, float>::value, void>::type
        setInfo(std::string tag, const T& v)
        {
            ret = -1;
            // bcf_hrec_set_val
            // bcf_update_info_flag(header->hdr, line, tag.c_str(), NULL, 1);
            int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
            if (bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_INT & 0xff))
                ret = bcf_update_info_int32(header->hdr, line, tag.c_str(), &v, 1);
            else if (bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_REAL & 0xff))
                ret = bcf_update_info_float(header->hdr, line, tag.c_str(), &v, 1);
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
        setInfo(std::string tag, const T& v)
        {
            ret = -1;
            // bcf_update_info_flag(header->hdr, line, tag.c_str(), NULL, 1);
            int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
            if (bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_INT & 0xff))
                ret = bcf_update_info_int32(header->hdr, line, tag.c_str(), v.data(), v.size());
            else if (bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_REAL & 0xff))
                ret = bcf_update_info_float(header->hdr, line, tag.c_str(), v.data(), v.size());
            else if (bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_STR & 0xff))
                ret = bcf_update_info_string(header->hdr, line, tag.c_str(), v.data());
            if (ret < 0)
            {
                throw std::runtime_error("couldn't set " + tag +
                                         " for this variant.\nplease add the tag in header first.\n");
            }
        }

        void removeInfo(std::string tag)
        {
            ret = -1;
            int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
            if (bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_INT & 0xff))
                ret = bcf_update_info_int32(header->hdr, line, tag.c_str(), NULL, 0);
            else if (bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_REAL & 0xff))
                ret = bcf_update_info_float(header->hdr, line, tag.c_str(), NULL, 0);
            else if (bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_STR & 0xff))
                ret = bcf_update_info_string(header->hdr, line, tag.c_str(), NULL);
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
        setGenotypes(const T& gv, bool phased = false)
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
        setFormat(std::string tag, T v)
        {
            ret = -1;
            float v2 = v;
            int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
            if (bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_INT & 0xff))
                ret = bcf_update_format_int32(header->hdr, line, tag.c_str(), &v, 1);
            else if (bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_REAL & 0xff))
                ret = bcf_update_format_float(header->hdr, line, tag.c_str(), &v2, 1);
            if (ret < 0)
                throw std::runtime_error("couldn't set format " + tag + " correctly.\n");
        }

        template <typename T>
        typename std::enable_if<std::is_same<T, std::string>::value || std::is_same<T, std::vector<char>>::value ||
                                    std::is_same<T, std::vector<int>>::value ||
                                    std::is_same<T, std::vector<float>>::value,
                                void>::type
        setFormat(std::string tag, const T& v)
        {
            ret = -1;
            int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
            if (bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_INT & 0xff))
                ret = bcf_update_format_int32(header->hdr, line, tag.c_str(), v.data(), v.size());
            else if (bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_STR & 0xff))
                ret = bcf_update_format_char(header->hdr, line, tag.c_str(), v.data(), v.size());
            else if (bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_REAL & 0xff))
                ret = bcf_update_format_float(header->hdr, line, tag.c_str(), v.data(), v.size());
            if (ret < 0)
                throw std::runtime_error("couldn't set format " + tag + " correctly.\n");
        }

        void addLineFromString(const std::string& vcfline)
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

        inline bool isNoneMissing() const
        {
            return noneMissing;
        }

        inline bool isSV() const
        {
            if (bcf_get_info(header->hdr, line, "SVTYPE") == NULL)
                return false;
            else
                return true;
        }

        inline bool isIndel() const
        {
            // REF has multiple allels
            if (REF().length() > 1 && !isSV())
                return true;
            for (int i = 1; i < line->n_allele; i++)
            {
                std::string alt(line->d.allele[i]);
                if (alt == ".")
                    return true;
                if (alt.length() != REF().length() && !isSV())
                    return true;
            }
            return false;
        }

        inline bool isDeletion() const
        {
            if (ALT().length() > 1)
                return false;
            if (!isIndel())
                return false;
            if (ALT().length() == 0)
                return true;
            else if (ALT()[0] == '.')
                return true;
            if (REF().length() > ALT().length())
                return true;
            return false;
        }

        inline bool isSNP() const
        {
            // REF has multiple allels
            if (REF().length() > 1)
                return false;
            for (int i = 1; i < line->n_allele; i++)
            {
                std::string snp(line->d.allele[i]);
                if (!(snp == "A" || snp == "C" || snp == "G" || snp == "T"))
                {
                    return false;
                }
            }
            return true;
        }

        inline std::string CHROM() const
        {
            return std::string(bcf_hdr_id2name(header->hdr, line->rid));
        }

        // 0-based
        inline int64_t POS() const
        {
            return line->pos;
        }

        // 0-based start of all type of variants
        inline int64_t Start() const
        {
            return line->pos;
        }

        // for SV variants
        inline int64_t End() const
        {
            return line->pos + line->rlen;
        }

        inline std::string REF() const
        {
            return std::string(line->d.allele[0]);
        }

        inline std::string ALT() const
        {
            std::string s;
            for (int i = 1; i < line->n_allele; i++)
            {
                s += std::string(line->d.allele[i]) + ",";
            }
            if (s.length() > 1)
                s.pop_back();
            return s;
        }

        inline float QUAL()
        {
            if (bcf_float_is_missing(line->qual))
            {
                noneMissing = false;
                return bcf_float_missing;
            }
            else
            {
                return line->qual;
            }
        }

        inline std::string FILTER()
        {
            if (line->d.n_flt == 0)
            {
                return ".";
            }
            else if (line->d.n_flt == 1)
            {
                return std::string(bcf_hdr_int2id(header->hdr, BCF_DT_ID, line->d.flt[0]));
            }
            else
            {
                std::string s;
                for (int i = 1; i < line->d.n_flt; i++)
                {
                    s += std::string(bcf_hdr_int2id(header->hdr, BCF_DT_ID, line->d.flt[i])) + ",";
                }
                s.pop_back();
                return s;
            }
        }

        inline bool allPhased() const
        {
            return isAllPhased;
        }

        inline int ploidy() const
        {
            return nploidy;
        }

        inline std::tuple<int, int> shapeOfQuery() const
        {
            return std::make_tuple(header->nsamples, shape1);
        }

    private:
        std::shared_ptr<BcfHeader> header;
        bcf1_t* line = bcf_init(); // current bcf record
        bcf_hdr_t* hdr_d;          // a dup header by bcf_hdr_dup(header->hdr)
        bcf_fmt_t* fmt = NULL;
        bcf_info_t* info = NULL;
        int32_t* gts = NULL;
        int ndst, ret;
        kstring_t s = {0, 0, NULL}; // kstring
        bool noneMissing = true;    // whenever parsing a tag have to reset this variable
        bool isAllPhased = false;
        int nploidy = 0;
        int shape1 = 0;
    };

    /**
     * @class BcfReader
     * @brief Stream in variants from vcf/bcf file or stdin
     * @note  nothing important
     **/
    class BcfReader
    {
    public:
        /**
         *  @brief construct a vcf/bcf reader from file.
         *  @param fname_   the input vcf/bcf with suffix vcf(.gz) or bcf(.gz)
         */
        BcfReader(const std::string& fname_) : fname(fname_)
        {
            fp = hts_open(fname.c_str(), "r");
            header.hdr = bcf_hdr_read(fp);
            header.nsamples = bcf_hdr_nsamples(header.hdr);
        }

        /**
         *  @brief construct a vcf/bcf reader with subset samples
         *  @param fname_   the input vcf/bcf with suffix vcf(.gz) or bcf(.gz)
         *  @param samples  LIST samples to include or exclude as a comma-separated string. \n
         *                  LIST : select samples in list \n
         *                  ^LIST : exclude samples from list \n
         *                  "-" : include all samples \n
         *                  "" : exclude all samples
         */
        BcfReader(const std::string& fname_, const std::string& samples) : fname(fname_)
        {
            fp = hts_open(fname.c_str(), "r");
            header.hdr = bcf_hdr_read(fp);
            header.setSamples(samples);
            header.nsamples = bcf_hdr_nsamples(header.hdr);
        }

        /**
         *  @brief construct a vcf/bcf reader with subset samples in target region
         *  @param fname_   the input vcf/bcf with suffix vcf(.gz) or bcf(.gz)
         *  @param samples  LIST samples to include or exclude as a comma-separated string. \n
         *                  LIST : select samples in list \n
         *                  ^LIST : exclude samples from list \n
         *                  "-" : include all samples \n
         *                  "" : exclude all samples
         *  @param region samtools-like region "chr:start-end"
         */
        BcfReader(const std::string& fname_, const std::string& samples, const std::string& region) : fname(fname_)
        {
            fp = hts_open(fname.c_str(), "r");
            header.hdr = bcf_hdr_read(fp);
            header.setSamples(samples);
            header.nsamples = bcf_hdr_nsamples(header.hdr);
            setRegion(region);
        }


        virtual ~BcfReader()
        {
            if (fp)
                hts_close(fp);
            if (itr)
                hts_itr_destroy(itr);
        }

        inline int setThreads(int n)
        {
            return hts_set_threads(fp, n);
        }

        // 1. check and load index first
        // 2. query iterval region
        void setRegion(const std::string& region)
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
        void reset();

        bool getNextVariant(BcfRecord& r)
        {
            int ret = -1;
            if (itr != NULL)
            {
                if (isBcf)
                {
                    ret = bcf_itr_next(fp, itr, r.line);
                    bcf_unpack(r.line, BCF_UN_ALL);
                    return (ret >= 0);
                }
                else
                {
                    int slen = tbx_itr_next(fp, tidx, itr, &s);
                    if (slen > 0)
                    {
                        ret = vcf_parse(&s, r.header->hdr, r.line); // ret > 0, error
                        bcf_unpack(r.line, BCF_UN_ALL);
                    }
                    return (ret <= 0) && (slen > 0);
                }
            }
            else
            {
                ret = bcf_read(fp, r.header->hdr, r.line);
                // unpack record immediately. not lazy
                bcf_unpack(r.line, BCF_UN_ALL);
                return (ret == 0);
            }
        }

        BcfHeader header; // bcf header

    private:
        htsFile* fp = NULL;         // hts file
        hts_idx_t* hidx = NULL;     // hts index file
        tbx_t* tidx = NULL;         // .tbi .csi index file for vcf files
        hts_itr_t* itr = NULL;      // hts records iterator
        kstring_t s = {0, 0, NULL}; // kstring
        std::string fname;
        bool isBcf; // if the input file is bcf or vcf;
    };

    /**
     * @class BcfWriter
     * @brief Stream out variants to vcf/bcf file or stdout
     * @note  nothing important
     **/
    class BcfWriter
    {
    public:
        /**
         * @brief          Open VCF/BCF file for writing. The format is infered from file's suffix
         * @param fname    The file name or "-" for stdin/stdout. For indexed files
         */
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

        /**
         * @brief          Open VCF/BCF file for writing using given mode
         * @param fname    The file name or "-" for stdin/stdout. For indexed files
         * @param mode     Mode matching \n
         *                 [w]b  .. compressed BCF \n
         *                 [w]bu .. uncompressed BCF \n
         *                 [w]z  .. compressed VCF \n
         *                 [w]   .. uncompressed VCF
         */
        BcfWriter(const std::string& fname_, const std::string& mode) : fname(fname_)
        {
            fp = hts_open(fname.c_str(), mode.c_str());
        }

        virtual ~BcfWriter()
        {
            hts_close(fp);
            bcf_destroy(b);
        }

        void initalHeader(std::string version = "VCF4.1")
        {
            header = BcfHeader();
            header.hdr = bcf_hdr_init("w");
            header.setVersion(version);
        }

        // make a copy of given header
        void initalHeader(const BcfHeader& h)
        {
            header = BcfHeader();
            header.hdr = bcf_hdr_dup(h.hdr); // make a copy of given header
            header.nsamples = bcf_hdr_nsamples(header.hdr);
            if (header.hdr == NULL)
                throw std::runtime_error("couldn't copy the header from another vcf.\n");
        }

        void writeLine(const std::string& vcfline)
        {
            if (!isHeaderWritten)
                writeHeader();
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

        bool writeHeader()
        {
            ret = bcf_hdr_write(fp, header.hdr);
            if (ret == 0)
                return isHeaderWritten = true;
            else
                return false;
        }

        inline bool writeRecord(BcfRecord& v)
        {
            if (!isHeaderWritten)
                writeHeader();
            if (bcf_write(fp, v.header->hdr, v.line) < 0)
                return false;
            else
                return true;
        }

        BcfHeader header = BcfHeader(); // bcf header


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
