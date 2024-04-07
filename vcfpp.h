/*******************************************************************************
 * @file        https://github.com/Zilong-Li/vcfpp/vcfpp.h
 * @author      Zilong Li
 * @email       zilong.dk@gmail.com
 * @version     v0.3.6
 * @breif       a single C++ file for manipulating VCF
 * Copyright (C) 2022-2023.The use of this code is governed by the LICENSE file.
 ******************************************************************************/

/*! \mainpage The documentation of the single C++ file *vcfpp.h* for manipulating VCF/BCF
 *
 * \section intro_sec Introduction
 *
 * This project https://github.com/Zilong-Li/vcfpp introduces a single C++ file as interface to the basic
 * htslib, which can be easily included in a C++ program for scripting high-performance genomic analyses.
 *
 * - vcfpp.BcfHeader keeps track of the header information in  VCF/BCF
 * - vcfpp.BcfRecord keeps track of the variants information in VCF/BCF
 * - vcfpp.BcfReader streams in variants from VCF/BCF file or stdin
 * - vcfpp.BcfWriter streams out variants to VCF/BCF file or stdout
 *
 * \section install_sec Installation
 *
 * - <EM> include "vcfpp.h" </EM> to your program and compile it by <EM> g++ my.cpp -std=c++11 -Wall -I. -lhts
 * - </EM>
 * - make sure you have https://github.com/samtools/htslib installed on your system and the it is in your
 * environment.
 *
 *
 * \copyright Copyright (C) 2022 Zilong Li . This project is governed by the LICENSE file in
 * https://github.com/Zilong-Li/vcfpp.
 *
 *
 */

#pragma once

#ifndef VCFPP_H_
#    define VCFPP_H_

#    include <iostream>
#    include <memory>
#    include <string>
#    include <type_traits>
#    include <vector>

// make sure you have htslib installed
extern "C"
{
#    include <htslib/kstring.h>
#    include <htslib/tbx.h>
#    include <htslib/vcf.h>
#    include <htslib/vcfutils.h>
}

namespace vcfpp
{
template<typename T>
using isValidFMT =
    typename std::enable_if<std::is_same<T, std::string>::value || std::is_same<T, std::vector<char>>::value
                                || std::is_same<T, std::vector<int>>::value
                                || std::is_same<T, std::vector<float>>::value,
                            bool>::type;

template<typename T>
using isValidInfo =
    typename std::enable_if<std::is_same<T, std::string>::value || std::is_same<T, std::vector<int>>::value
                                || std::is_same<T, std::vector<float>>::value,
                            bool>::type;

template<typename T>
using isInfoVector = typename std::enable_if<std::is_same<T, std::vector<int>>::value
                                                 || std::is_same<T, std::vector<float>>::value,
                                             bool>::type;

template<typename T>
using isScalar = typename std::enable_if<std::is_same<T, int>::value || std::is_same<T, float>::value
                                             || std::is_same<T, double>::value,
                                         bool>::type;

template<typename T>
using isString = typename std::enable_if<std::is_same<T, std::string>::value, void>::type;

template<typename T>
using isValidGT = typename std::enable_if<std::is_same<T, std::vector<bool>>::value
                                              || std::is_same<T, std::vector<char>>::value,
                                          bool>::type;

template<typename T>
using isFormatVector = typename std::enable_if<std::is_same<T, std::vector<float>>::value
                                                   || std::is_same<T, std::vector<char>>::value
                                                   || std::is_same<T, std::vector<int>>::value,
                                               bool>::type;

template<typename T>
isScalar<T> isMissing(T const & v)
{
    if(std::is_same<T, float>::value)
    {
        return bcf_float_is_missing(v);
    }
    else if(std::is_same<T, int>::value)
    {
        return bcf_int32_missing(v);
    }
}

// test if a string is end with another string
inline bool isEndWith(std::string const & s, std::string const & e)
{
    if(s.length() >= e.length())
    {
        return (0 == s.compare(s.length() - e.length(), e.length(), e));
    }
    else
    {
        return false;
    }
}

// string split by separator
inline std::vector<std::string> split_string(const std::string & s, const std::string & separators)
{
    std::vector<std::string> ret;
    bool is_seperator[256] = {false};
    for(auto & ch : separators)
    {
        is_seperator[(unsigned int)ch] = true;
    }
    int begin = 0;
    for(int i = 0; i <= (int)s.size(); i++)
    {
        if(is_seperator[(uint8_t)s[i]] || i == (int)s.size())
        {
            ret.push_back(std::string(s.begin() + begin, s.begin() + i));
            begin = i + 1;
        }
    }
    return ret;
}

// deleter for htsFile
struct htsFile_close
{
    void operator()(htsFile * x)
    {
        if(x) hts_close(x);
    }
};

// deleter for variant
struct variant_close
{
    void operator()(bcf1_t * x)
    {
        if(x) bcf_destroy(x);
    }
};


// deleter for bcf header
struct bcf_hdr_close
{
    void operator()(bcf_hdr_t * x)
    {
        if(x) bcf_hdr_destroy(x);
    }
};
    
/**
 * @class BcfHeader
 * @brief Object represents header of the VCF, offering methods to access and modify the tags
 * @note  BcfHeader has 3 friends, BcfReader, BcfWriter and BcfRecord.
 **/
class BcfHeader
{
    friend class BcfRecord;
    friend class BcfReader;
    friend class BcfWriter;

  private:
    bcf_hdr_t * hdr = NULL; // bcf header
    bcf_hrec_t * hrec = NULL; // populate header

  public:
    BcfHeader() {}

    ~BcfHeader()
    {
        if(hrec) bcf_hrec_destroy(hrec);
        if(hdr) bcf_hdr_destroy(hdr);
    }

    /** @brief stream out the header */
    friend std::ostream & operator<<(std::ostream & out, const BcfHeader & h)
    {
        out << h.asString();
        return out;
    }

    // TODO: check if the value is valid for vcf specification

    /** @brief add INFO field to header
     *  @param id tag name in INFO field
     *  @param number Number of elements
     *  @param type Integer, Floast or String
     *  @param description Description of the tag
     *  */
    inline void addINFO(const std::string & id,
                        const std::string & number,
                        const std::string & type,
                        const std::string & description)
    {
        addLine("##INFO=<ID=" + id + ",Number=" + number + ",Type=" + type + ",Description=\"" + description
                + "\">");
    }

    /** @brief add FORMAT field to header
     *  @param id tag name in FORMAT field
     *  @param number Number of elements
     *  @param type Integer, Floast or String
     *  @param description Description of the tag
     *  */
    inline void addFORMAT(const std::string & id,
                          const std::string & number,
                          const std::string & type,
                          const std::string & description)
    {
        addLine("##FORMAT=<ID=" + id + ",Number=" + number + ",Type=" + type + ",Description=\"" + description
                + "\">");
    }

    /**
     * @brief add one FILTER line to header
     * @param id FILTER name
     * @param description Description of the FILTER
     * */
    inline void addFILTER(const std::string & id, const std::string & description)
    {
        addLine("##FILTER=<ID=" + id + ",Description=\"" + description + "\">");
    }

    /** @brief add contig to header
     *  @param id contig or chromosome name
     *  */
    inline void addContig(const std::string & id)
    {
        addLine("##contig=<ID=" + id + ">");
    }

    /**
     * @brief add one line to header
     * */
    inline void addLine(const std::string & str)
    {
        int ret = 0;
        ret = bcf_hdr_append(hdr, str.c_str());
        if(ret != 0) throw std::runtime_error("could not add " + str + " to header\n");
        ret = bcf_hdr_sync(hdr);
        if(ret != 0) throw std::runtime_error("could not add " + str + " to header\n");
    }

    /** @brief add one sample name to header */
    inline void addSample(const std::string & sample) const
    {
        bcf_hdr_add_sample(hdr, sample.c_str());
        if(bcf_hdr_sync(hdr) != 0)
        {
            throw std::runtime_error("couldn't add the sample.\n");
        }
    }

    /** @brief return header as a string */
    inline std::string asString() const
    {
        kstring_t s = {0, 0, NULL}; // kstring
        if(bcf_hdr_format(hdr, 0, &s) == 0) // append header string to s.s! append!
        {
            std::string out(s.s, s.l);
            free(s.s);
            return out;
        }
        else
            throw std::runtime_error("failed to convert formatted header to string");
    }

    /** @brief return all samples in a string vector */
    std::vector<std::string> getSamples() const
    {
        std::vector<std::string> vec;
        for(int i = 0; i < bcf_hdr_nsamples(hdr); i++)
        {
            vec.push_back(std::string(hdr->samples[i]));
        }
        return vec;
    }

    /** @brief return all contig/chromosome names in a string vector */
    std::vector<std::string> getSeqnames() const
    {
        int ret = 0;
        const char ** seqs = bcf_hdr_seqnames(hdr, &ret);
        if(ret == 0) printf("there is no contig id in the header!\n");
        std::vector<std::string> vec;
        for(int i = 0; i < ret; i++)
        {
            vec.push_back(std::string(seqs[i]));
        }
        free(seqs);
        return vec;
    }

    /**
     *  @brief get the type of a given tag
     *  @param tag in the FORMAT
     *  @return 1: int; 2: float; 3: string; 0: error;
     */
    inline int getFormatType(std::string tag) const
    {
        int tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, tag.c_str());
        if(tag_id < 0) return 0;
        if(bcf_hdr_id2type(hdr, BCF_HL_FMT, tag_id) == (BCF_HT_INT & 0xff))
            return 1;
        else if(bcf_hdr_id2type(hdr, BCF_HL_FMT, tag_id) == (BCF_HT_REAL & 0xff))
            return 2;
        else if(bcf_hdr_id2type(hdr, BCF_HL_FMT, tag_id) == (BCF_HT_STR & 0xff))
            return 3;
        else
            return 0;
    }

    /** @brief remove a contig tag from header */
    inline void removeContig(std::string tag) const
    {

        bcf_hdr_remove(hdr, BCF_HL_CTG, tag.c_str());
    }

    /** @brief remove a INFO tag from header */
    inline void removeInfo(std::string tag) const
    {
        bcf_hdr_remove(hdr, BCF_HL_INFO, tag.c_str());
    }

    /** @brief remove a FORMAT tag from header */
    inline void removeFormat(std::string tag) const
    {
        bcf_hdr_remove(hdr, BCF_HL_FMT, tag.c_str());
    }

    /** @brief remove a FILTER tag from header */
    inline void removeFilter(std::string tag) const
    {
        bcf_hdr_remove(hdr, BCF_HL_FLT, tag.c_str());
    }

    /**
     * @brief explicitly set samples to be extracted
     * @param samples samples to include or exclude  as a comma-separated string
     * */
    inline void setSamples(const std::string & samples) const
    {
        int ret = 0;
        ret = bcf_hdr_set_samples(hdr, samples.c_str(), 0);
        if(ret != 0)
        {
            throw std::runtime_error("the " + std::to_string(ret)
                                     + "-th sample are not in the VCF.\nparameter samples:" + samples);
        }
    }

    /** @brief set the VCF version */
    inline void setVersion(const std::string & version) const
    {
        bcf_hdr_set_version(hdr, version.c_str());
    }

    /** @brief return the number of samples in the VCF */
    inline int nSamples() const
    {
        return bcf_hdr_nsamples(hdr);
    }
};

/**
 * @class BcfRecord
 * @brief Object represents a variant record in the VCF, offering methods to access and modify fields.
 * @note  BcfRecord has to be associated with a BcfHeader object and needs to be filled in by calling
 *BcfReader.getNextVariant function.
 **/
class BcfRecord
{
    friend class BcfReader;
    friend class BcfWriter;

  private:
    BcfHeader* header;
    std::shared_ptr<bcf1_t> line = std::shared_ptr<bcf1_t>(bcf_init(), variant_close()); // variant
    bcf_hdr_t * hdr_d; // a dup header by bcf_hdr_dup(header->hdr)
    bcf_fmt_t * fmt = NULL;
    bcf_info_t * info = NULL;
    int32_t * gts = NULL;
    int ndst, ret, nsamples;
    bool noneMissing = true; // whenever parsing a tag have to reset this variable
    bool isAllPhased = false;
    int nploidy = 0; // the number of ploidy
    int nvalues = 0; // the number of values for a tag in FORMAT

  public:
    /// if there is "." in GT for the sample, then it's coded as missing (TRUE)
    std::vector<char> isGenoMissing;

  public:
    /// empty constructor. call init() afterwards
    BcfRecord() {}

    /// constructor with a given BcfHeader object
    BcfRecord(BcfHeader & h)
    {
        init(h);
    }

    ~BcfRecord() {}

    /// initilize a BcfRecord object by pointing to another BcfHeader object
    void init(BcfHeader & h)
    {
        header = &h;
        if(!header->hdr) throw std::runtime_error("please initial header first\n");
        nsamples = header->nSamples();
        typeOfGT.resize(nsamples);
        gtPhase.resize(nsamples, 0);
    }

    /** @brief stream out the variant */
    friend std::ostream & operator<<(std::ostream & out, const BcfRecord & v)
    {
        out << v.asString();
        return out;
    }

    /** @brief return current variant as raw string */
    inline std::string asString() const
    {
        kstring_t s = {0, 0, NULL}; // kstring
        if(vcf_format(header->hdr, line.get(), &s) == 0)
        {
            std::string out(s.s, s.l);
            free(s.s);
            return out;
        }
        else
            throw std::runtime_error("couldn't format current record into a string.\n");
    }

    /**
     * @brief fill in the input vector with genotypes of 0 and 1. only works for ploidy<=2. Genotypes with
     * missing allele is coded as heterozygous
     * @param v valid input includes vector<bool> and vector<char> type
     * @return bool
     * @note  use isNoneMissing() to check if all genotypes are with no missingness. Alternatively, one can
     * use vector<int> as the input type as noted in the other overloading function getGenotypes().
     * */
    template<typename T>
    isValidGT<T> getGenotypes(T & v)
    {
        ndst = 0;
        ret = bcf_get_genotypes(header->hdr, line.get(), &gts, &ndst);
        if(ret <= 0) throw std::runtime_error("genotypes not present");
        // if nploidy is not set manually. find the max nploidy using the first variant (eg. 2) resize v as
        // max(nploidy)
        // * nsamples (ret) NOTE: if ret == nsamples and only one sample then haploid
        if(ret != nploidy * nsamples && nploidy > 0)
        {
            // rare case if nploidy is set manually. eg. only one sample. the first variant is '1' but the
            // second is '1/0'. ret = 1 but nploidy should be 2
            v.resize(nploidy * nsamples);
        }
        else
        {
            v.resize(ret);
            nploidy = ret / nsamples;
        }
        // work with nploidy == 1, haploid?
        isGenoMissing.assign(nsamples, 0);
        int i, j, nphased = 0;
        noneMissing = true;
        fmt = bcf_get_fmt(header->hdr, line.get(), "GT");
        int nploidy_cur = ret / nsamples; // requires nploidy_cur <= nploidy
        for(i = 0; i < nsamples; i++)
        {
            // check and fill in typeOfGT; only supports SNPs now. check vcfstats.c for inspiration
            typeOfGT[i] = bcf_gt_type(fmt, i, NULL, NULL);
            if(typeOfGT[i] == GT_UNKN)
            {
                noneMissing = false; // set missing as het, user should check if isNoneMissing();
                isGenoMissing[i] = 1;
                v[i * nploidy + 0] = 1;
                for(j = 1; j < nploidy_cur; j++) v[i * nploidy + j] = 0;
                continue;
            }

            for(j = 0; j < nploidy_cur; j++)
            {
                // TODO: right now only parse 0 and 1, ie max(nploidy)=2; other values 2,3... will be set to 1
                v[i * nploidy + j] = bcf_gt_allele(gts[j + i * nploidy_cur]) != 0;
            }
            if(nploidy == 2)
            {
                gtPhase[i] = (gts[1 + i * nploidy_cur] & 1) == 1;
                nphased += gtPhase[i];
            }
        }
        if(nphased == nsamples)
            isAllPhased = true;
        else
            isAllPhased = false;
        return true;
    }

    /**
     * @brief fill in the input vector with genotype values, 0, 1 or -9 (missing).
     * @param v valid input is vector<int> type
     * @return bool
     * @note this function provides full capability to handl all kinds of genotypes in multi-ploidy data with
     * the cost of more spae than the other function. missing allele is set as -9.
     * */
    bool getGenotypes(std::vector<int> & v)
    {
        ndst = 0;
        ret = bcf_get_genotypes(header->hdr, line.get(), &gts, &ndst);
        if(ret <= 0)
            throw std::runtime_error(
                "genotypes not present. make sure you initilized the variant object first\n");
        v.resize(ret);
        isGenoMissing.assign(nsamples, 0);
        nploidy = ret / nsamples;
        int i, j, nphased = 0;
        noneMissing = true;
        for(i = 0; i < nsamples; i++)
        {
            int32_t * ptr = gts + i * nploidy;
            int is_phased = 0;
            for(j = 0; j < nploidy; j++)
            {
                // if true, the sample has smaller ploidy
                if(ptr[j] == bcf_int32_vector_end) break;
                // missing allele
                if(bcf_gt_is_missing(ptr[j]))
                {
                    noneMissing = false;
                    isGenoMissing[i] = 1;
                    v[i * nploidy + j] = -9;
                    continue;
                }
                v[i * nploidy + j] = bcf_gt_allele(ptr[j]);
                is_phased += bcf_gt_is_phased(ptr[j]);
            }
            if(nploidy == is_phased)
            {
                gtPhase[i] = true;
                nphased += gtPhase[i];
            }
        }
        if(nphased == nsamples)
            isAllPhased = true;
        else
            isAllPhased = false;
        return true;
    }

    /**
     * @brief get tag value in FORMAT
     * @param tag valid tag name in FORMAT column declared in the VCF header
     * @param v valid input include vector of float, char, int type
     * @return bool
     * */
    template<typename T, typename S = typename T::value_type>
    isFormatVector<T> getFORMAT(std::string tag, T & v)
    {
        fmt = bcf_get_fmt(header->hdr, line.get(), tag.c_str());
        if(!fmt) throw std::runtime_error("there is no " + tag + " in FORMAT of this variant.\n");
        nvalues = fmt->n;
        ndst = 0;
        S * dst = NULL;
        int tagid = header->getFormatType(tag);
        if(tagid == 1)
            ret = bcf_get_format_int32(header->hdr, line.get(), tag.c_str(), &dst, &ndst);
        else if(tagid == 2)
            ret = bcf_get_format_float(header->hdr, line.get(), tag.c_str(), &dst, &ndst);
        else if(tagid == 3)
            ret = bcf_get_format_char(header->hdr, line.get(), tag.c_str(), &dst, &ndst);
        else
            throw std::runtime_error("can not find the type of " + tag + " in the header file.\n");
        if(ret >= 0)
        {
            // user have to check if there is missing in the return v;
            v = std::vector<S>(dst, dst + ret);
            free(dst);
            return true;
        }
        else
        {
            throw std::runtime_error("failed to parse the " + tag + " format of this variant.\n");
        }
    }

    /**
     * @brief get tag value in FORMAT
     * @param tag valid tag name in FORMAT column declared in the VCF header
     * @param v vector of string
     * @return bool
     * */
    bool getFORMAT(std::string tag, std::vector<std::string> & v)
    {
        fmt = bcf_get_fmt(header->hdr, line.get(), tag.c_str());
        if(!fmt) throw std::runtime_error("there is no " + tag + " in FORMAT for this variant of ID=" + ID());
        nvalues = fmt->n;
        // if ndst < (fmt->n+1)*nsmpl; then realloc is involved
        ret = -1, ndst = 0;
        char ** dst = NULL;
        if(header->getFormatType(tag) == 3)
            ret = bcf_get_format_string(header->hdr, line.get(), tag.c_str(), &dst, &ndst);
        if(ret > 0)
        {
            v.clear();
            for(int i = 0; i < nsamples; i++)
            {
                // Ugly: GT field is considered to be a string by the VCF header but BCF represents it as INT.
                v.emplace_back(dst[i]);
            };
            free(dst);
            return true;
        }
        else
        {
            throw std::runtime_error("couldn't parse the " + tag + " format of this variant.\n");
        }
    }

    /**
     * @brief get tag value in INFO
     * @param tag valid tag name in INFO column declared in the VCF header
     * @param v valid input include vector of float, int type
     * @return bool
     * */
    template<typename T, typename S = typename T::value_type>
    isInfoVector<T> getINFO(std::string tag, T & v)
    {
        info = bcf_get_info(header->hdr, line.get(), tag.c_str());
        if(!info) throw std::runtime_error("there is no " + tag + " tag in INFO of this variant.\n");
        S * dst = NULL;
        ndst = 0;
        ret = -1;
        if(info->type == BCF_BT_INT8 || info->type == BCF_BT_INT16 || info->type == BCF_BT_INT32)
        {
            ret = bcf_get_info_int32(header->hdr, line.get(), tag.c_str(), &dst, &ndst);
        }
        else if(info->type == BCF_BT_FLOAT)
        {
            ret = bcf_get_info_float(header->hdr, line.get(), tag.c_str(), &dst, &ndst);
        }
        if(ret >= 0)
            v = std::vector<S>(dst, dst + ret); // user have to check if there is missing in the return v;
        else
            throw std::runtime_error("couldn't parse the " + tag + " format of this variant.\n");
        free(dst);
        return true;
    }

    /**
     * @brief get tag value in INFO
     * @param tag valid tag name in INFO column declared in the VCF header
     * @param v valid input include scalar value of float or int type
     * @return bool
     * */
    template<typename T>
    isScalar<T> getINFO(std::string tag, T & v)
    {
        info = bcf_get_info(header->hdr, line.get(), tag.c_str());
        if(!info) throw std::runtime_error("there is no " + tag + " tag in INFO of this variant.\n");
        // scalar value
        if(info->len == 1)
        {
            if(info->type == BCF_BT_INT8 || info->type == BCF_BT_INT16 || info->type == BCF_BT_INT32)
            {
                v = info->v1.i;
            }
            else if(info->type == BCF_BT_FLOAT)
            {
                v = info->v1.f;
            }
            return true;
        }
        else
        {
            throw std::runtime_error(tag + " has multiple values. please feed a vector instead.\n");
        }
    }

    /**
     * @brief get tag value in INFO
     * @param tag valid tag name in INFO column declared in the VCF header
     * @param v valid input is std::string
     * @return bool
     * */
    template<typename T>
    isString<T> getINFO(std::string tag, T & v)
    {
        info = bcf_get_info(header->hdr, line.get(), tag.c_str());
        if(!info) throw std::runtime_error("there is no " + tag + " tag in INFO of this variant.\n");
        if(info->type == BCF_BT_CHAR)
            v = std::string(reinterpret_cast<char *>(info->vptr), info->vptr_len);
        else
            throw std::runtime_error(tag + " has to be of string type\n");
    }

    /**
     * @brief set tag value for INFO
     * @param tag valid tag name in INFO column declared in the VCF header
     * @param v valid input include scalar value of float or int type
     * @return bool
     * */
    template<typename T>
    isScalar<T> setINFO(std::string tag, const T & v)
    {
        ret = -1;
        // bcf_hrec_set_val
        // bcf_update_info_flag(header.hdr, line, tag.c_str(), NULL, 1);
        int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
        if(bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_INT & 0xff))
            ret = bcf_update_info_int32(header->hdr, line.get(), tag.c_str(), &v, 1);
        else if(bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_REAL & 0xff))
        {
            float v2 = static_cast<float>(v);
            ret = bcf_update_info_float(header->hdr, line.get(), tag.c_str(), &v2, 1);
        }
        if(ret < 0)
        {
            throw std::runtime_error("couldn't set " + tag
                                     + " for this variant.\nplease add the tag in header first.\n");
        }
        else
        {
            return true;
        }
    }

    /**
     * @brief set tag value for INFO
     * @param tag valid tag name in INFO column declared in the VCF header
     * @param v valid input include vector<int> vector<float> std::string
     * @return boolean
     * */
    template<typename T>
    isValidInfo<T> setINFO(std::string tag, const T & v)
    {
        ret = -1;
        // bcf_update_info_flag(header.hdr, line, tag.c_str(), NULL, 1);
        int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
        if(bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_INT & 0xff))
            ret = bcf_update_info_int32(header->hdr, line.get(), tag.c_str(), v.data(), v.size());
        else if(bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_REAL & 0xff))
            ret = bcf_update_info_float(header->hdr, line.get(), tag.c_str(), v.data(), v.size());
        else if(bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_STR & 0xff))
            ret = bcf_update_info_string(header->hdr, line.get(), tag.c_str(), v.data());
        if(ret < 0)
            throw std::runtime_error("couldn't set " + tag
                                     + " for this variant.\nplease add the tag in header first.\n");
        else
            return true;
    }

    /// remove the given tag from INFO of the variant
    void removeINFO(std::string tag)
    {
        ret = -1;
        int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
        if(bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_INT & 0xff))
            ret = bcf_update_info_int32(header->hdr, line.get(), tag.c_str(), NULL, 0);
        else if(bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_REAL & 0xff))
            ret = bcf_update_info_float(header->hdr, line.get(), tag.c_str(), NULL, 0);
        else if(bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_STR & 0xff))
            ret = bcf_update_info_string(header->hdr, line.get(), tag.c_str(), NULL);
        if(ret < 0)
        {
            throw std::runtime_error("couldn't remove " + tag + " for this variant.\n");
        }
    }

    /**
     * @brief set genotypes from scratch even if genotypes not present
     * @param v  the genotypes of vector<int> type
     * @return bool
     * */
    bool setGenotypes(const std::vector<int> & v)
    {
        // bcf_gt_type
        int i, j, k;
        nploidy = v.size() / nsamples;
        gts = (int32_t *)malloc(v.size() * sizeof(int32_t));
        for(i = 0; i < nsamples; i++)
        {
            for(j = 0; j < nploidy; j++)
            {
                k = i * nploidy + j;
                if(v[k] == -9 || v[k] == bcf_int32_missing)
                    gts[k] = bcf_gt_missing;
                else if(gtPhase[i])
                    gts[k] = bcf_gt_phased(v[k]);
                else
                    gts[k] = bcf_gt_unphased(v[k]);
            }
        }
        if(bcf_update_genotypes(header->hdr, line.get(), gts, v.size()) < 0)
            throw std::runtime_error("couldn't set genotypes correctly.\n");
        return true;
    }

    /**
     * @brief set phasing status for all diploid samples using given vector
     * @param v valid input includes vector<char>
     * */
    void setPhasing(const std::vector<char> & v)
    {
        assert((int)v.size() == nsamples);
        gtPhase = v;
    }

    /// remove the given tag from FORMAT of the variant
    void removeFORMAT(std::string tag)
    {
        ret = -1;
        int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
        if(bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_INT & 0xff))
            ret = bcf_update_format_int32(header->hdr, line.get(), tag.c_str(), NULL, 0);
        else if(bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_STR & 0xff))
            ret = bcf_update_format_char(header->hdr, line.get(), tag.c_str(), NULL, 0);
        else if(bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_REAL & 0xff))
            ret = bcf_update_format_float(header->hdr, line.get(), tag.c_str(), NULL, 0);
        if(ret < 0) throw std::runtime_error("couldn't remove " + tag + " correctly.\n");
    }

    /**
     * @brief set tag values for all samples in FORMAT using given vector
     * @param tag valid tag name in FORMAT column declared in the VCF header
     * @param v valid input include vector<int>, vector<float>, vector<char>, std::string
     * @return bool
     * */
    template<typename T>
    isValidFMT<T> setFORMAT(std::string tag, const T & v)
    {
        ret = -1;
        int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
        if(bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_INT & 0xff))
            ret = bcf_update_format_int32(header->hdr, line.get(), tag.c_str(), v.data(), v.size());
        else if(bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_STR & 0xff))
            ret = bcf_update_format_char(header->hdr, line.get(), tag.c_str(), v.data(), v.size());
        else if(bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_REAL & 0xff))
            ret = bcf_update_format_float(header->hdr, line.get(), tag.c_str(), v.data(), v.size());
        if(ret < 0) throw std::runtime_error("couldn't set format " + tag + " correctly.\n");
        return true;
    }

    /**
     * @brief set tag for a single sample in FORMAT using given singular value. this works only when there is
     * one sample in the vcf
     * @param tag valid tag name in FORMAT column declared in the VCF header
     * @param v valid input include int, float or double
     * @return bool
     * */
    template<typename T>
    isScalar<T> setFORMAT(std::string tag, const T & v)
    {
        ret = -1;
        float v2 = v;
        int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
        if(bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_INT & 0xff))
            ret = bcf_update_format_int32(header->hdr, line.get(), tag.c_str(), &v, 1);
        else if(bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_REAL & 0xff))
            ret = bcf_update_format_float(header->hdr, line.get(), tag.c_str(), &v2, 1);
        if(ret < 0) throw std::runtime_error("couldn't set format " + tag + " correctly.\n");
        return true;
    }

    /** @brief add one variant record from given string*/
    void addLineFromString(const std::string & vcfline)
    {
        kstring_t s = {0, 0, NULL};
        kputsn(vcfline.c_str(), vcfline.length(), &s);
        ret = vcf_parse1(&s, header->hdr, line.get());
        free(s.s);
        if(ret > 0) throw std::runtime_error("error parsing: " + vcfline + "\n");
        if(line->errcode == BCF_ERR_CTG_UNDEF)
        {
            std::string contig(bcf_hdr_id2name(header->hdr, line->rid));
            hdr_d = bcf_hdr_dup(header->hdr);
            header->hrec = bcf_hdr_id2hrec(hdr_d, BCF_DT_CTG, 0, line->rid);
            if(header->hrec == NULL)
                throw std::runtime_error("contig" + contig + " unknow and not found in the header.\n");
            ret = bcf_hdr_add_hrec(header->hdr, header->hrec);
            if(ret < 0) throw std::runtime_error("error adding contig " + contig + " to header.\n");
            ret = bcf_hdr_sync(header->hdr);
        }
    }

    /** @brief if all samples have non missing values for the tag in FORMAT */
    inline bool isNoneMissing() const
    {
        return noneMissing;
    }

    /** @brief return boolean value indicates if current variant is Structual Variant or not */
    inline bool isSV() const
    {
        if(bcf_get_info(header->hdr, line.get(), "SVTYPE") == NULL)
            return false;
        else
            return true;
    }

    /** @brief return boolean value indicates if current variant is exclusively INDEL */
    inline bool isIndel() const
    {
        // REF has multiple allels
        if(REF().length() > 1 && !isSV()) return true;
        for(int i = 1; i < line->n_allele; i++)
        {
            std::string alt(line->d.allele[i]);
            if(alt == ".") return true;
            if(alt.length() != REF().length() && !isSV()) return true;
        }
        return false;
    }

    /** @brief return boolean value indicates if current variant is multiallelic sites */
    inline bool isMultiAllelics() const
    {
        if(line->n_allele <= 2) return false;
        return true;
    }

    /** @brief return boolean value indicates if current variant is exclusively multiallelic SNP sites */
    inline bool isMultiAllelicSNP() const
    {
        // skip REF with length > 1, i.e. INDEL
        if(REF().length() > 1 || line->n_allele <= 2) return false;
        for(int i = 1; i < line->n_allele; i++)
        {
            std::string snp(line->d.allele[i]);
            if(snp.length() != 1)
            {
                return false;
            }
        }
        return true;
    }

    /** @brief return boolean value indicates if current variant is exclusively biallelic SNP. Note ALT=* are
     * skipped */
    inline bool isSNP() const
    {
        // REF and ALT have multiple allels
        if(REF().length() > 1 || line->n_allele > 2) return false;
        std::string snp(line->d.allele[1]);
        if(!(snp == "A" || snp == "C" || snp == "G" || snp == "T"))
        {
            return false;
        }
        return true;
    }

    /** @brief return boolean value indicates if current variant has SNP type defined in vcf.h (htslib>=1.16)
     */
    inline bool hasSNP() const
    {
        int type = bcf_has_variant_types(line.get(), VCF_SNP, bcf_match_overlap);
        if(type < 0) throw std::runtime_error("something wrong with variant type\n");
        if(type == 0) return false;
        return true;
    }

    /// return boolean value indicates if current variant has INDEL type defined in vcf.h (htslib>=1.16)
    inline bool hasINDEL() const
    {
        int type = bcf_has_variant_types(line.get(), VCF_INDEL, bcf_match_overlap);
        if(type < 0) throw std::runtime_error("something wrong with variant type\n");
        if(type == 0) return false;
        return true;
    }

    /// return boolean value indicates if current variant has INS type defined in vcf.h (htslib>=1.16)
    inline bool hasINS() const
    {
        int type = bcf_has_variant_types(line.get(), VCF_INS, bcf_match_overlap);
        if(type < 0) throw std::runtime_error("something wrong with variant type\n");
        if(type == 0) return false;
        return true;
    }

    /// return boolean value indicates if current variant has DEL type defined in vcf.h (htslib>=1.16)
    inline bool hasDEL() const
    {
        int type = bcf_has_variant_types(line.get(), VCF_DEL, bcf_match_overlap);
        if(type < 0) throw std::runtime_error("something wrong with variant type\n");
        if(type == 0) return false;
        return true;
    }

    /// return boolean value indicates if current variant has MNP type defined in vcf.h (htslib>=1.16)
    inline bool hasMNP() const
    {
        int type = bcf_has_variant_types(line.get(), VCF_MNP, bcf_match_overlap);
        if(type < 0) throw std::runtime_error("something wrong with variant type\n");
        if(type == 0) return false;
        return true;
    }

    /// return boolean value indicates if current variant has BND type defined in vcf.h (htslib>=1.16)
    inline bool hasBND() const
    {
        int type = bcf_has_variant_types(line.get(), VCF_BND, bcf_match_overlap);
        if(type < 0) throw std::runtime_error("something wrong with variant type\n");
        if(type == 0) return false;
        return true;
    }

    /// return boolean value indicates if current variant has OTHER type defined in vcf.h (htslib>=1.16)
    inline bool hasOTHER() const
    {
        int type = bcf_has_variant_types(line.get(), VCF_OTHER, bcf_match_overlap);
        if(type < 0) throw std::runtime_error("something wrong with variant type\n");
        if(type == 0) return false;
        return true;
    }

    /// return boolean value indicates if current variant has OVERLAP type defined in vcf.h (htslib>=1.16)
    inline bool hasOVERLAP() const
    {
        int type = bcf_has_variant_types(line.get(), VCF_OVERLAP, bcf_match_overlap);
        if(type < 0) throw std::runtime_error("something wrong with variant type\n");
        if(type == 0) return false;
        return true;
    }

    /** @brief return CHROM name */
    inline std::string CHROM() const
    {
        return std::string(bcf_hdr_id2name(header->hdr, line->rid));
    }

    /** @brief return ID field */
    inline std::string ID() const
    {
        return std::string(line->d.id);
    }

    /** @brief return 1-base position */
    inline int64_t POS() const
    {
        return line->pos + 1;
    }

    /** @brief modify CHROM value */
    inline void setCHR(const char * chr)
    {
        line->rid = bcf_hdr_name2id(header->hdr, chr);
    }

    /** @brief modify position given 1-based value */
    inline void setPOS(int64_t p)
    {
        line->pos = p - 1;
    }

    /** @brief update ID */
    inline void setID(const char * s)
    {
        bcf_update_id(header->hdr, line.get(), s);
    }

    /** @brief set REF and ALT alleles given a string seperated by comma */
    inline void setRefAlt(const char * alleles_string)
    {
        bcf_update_alleles_str(header->hdr, line.get(), alleles_string);
    }

    /** @brief return 0-base start of the variant (can be any type) */
    inline int64_t Start() const
    {
        return line->pos;
    }

    /** @brief return 0-base end of the variant (can be any type) */
    inline int64_t End() const
    {
        return line->pos + line->rlen;
    }

    /** @brief return raw REF alleles as string */
    inline std::string REF() const
    {
        return std::string(line->d.allele[0]);
    }

    /** @brief swap REF and ALT for biallelic SNP */
    inline void swap_REF_ALT()
    {
        if(!isMultiAllelicSNP()) std::swap(line->d.allele[0], line->d.allele[1]);
    }

    /** @brief return raw ALT alleles as string */
    inline std::string ALT() const
    {
        std::string s;
        for(int i = 1; i < line->n_allele; i++)
        {
            s += std::string(line->d.allele[i]) + ",";
        }
        if(s.length() > 1) s.pop_back();
        return s;
    }

    /** @brief return QUAL value */
    inline float QUAL()
    {
        if(bcf_float_is_missing(line->qual))
        {
            noneMissing = false;
            return bcf_float_missing;
        }
        else
        {
            return line->qual;
        }
    }

    /** @brief return raw FILTER column as string */
    inline std::string FILTER()
    {
        if(line->d.n_flt == 0)
        {
            return ".";
        }
        else if(line->d.n_flt == 1)
        {
            return std::string(bcf_hdr_int2id(header->hdr, BCF_DT_ID, line->d.flt[0]));
        }
        else
        {
            std::string s;
            for(int i = 1; i < line->d.n_flt; i++)
            {
                s += std::string(bcf_hdr_int2id(header->hdr, BCF_DT_ID, line->d.flt[i])) + ",";
            }
            s.pop_back();
            return s;
        }
    }

    /** @brief return raw INFO column as string. recommend to use getINFO for specific tag. */
    inline std::string allINFO()
    {
        int32_t max_dt_id = header->hdr->n[BCF_DT_ID];
        kstring_t * s = (kstring_t *)calloc(1, sizeof(kstring_t));
        if(line->n_info)
        {
            int first = 1;
            for(int i = 0; i < line->n_info; ++i)
            {
                bcf_info_t * z = &line->d.info[i];
                if(!z->vptr) continue;
                if(!first) kputc(';', s);
                first = 0;
                if(z->key < 0 || z->key >= max_dt_id || header->hdr->id[BCF_DT_ID][z->key].key == NULL)
                    throw std::runtime_error("Invalid BCF and wrong INFO tag");
                kputs(header->hdr->id[BCF_DT_ID][z->key].key, s);
                if(z->len <= 0) continue;
                kputc('=', s);
                if(z->len == 1)
                {
                    switch(z->type)
                    {
                        case BCF_BT_INT8:
                            if(z->v1.i == bcf_int8_missing)
                                kputc('.', s);
                            else
                                kputw(z->v1.i, s);
                            break;
                        case BCF_BT_INT16:
                            if(z->v1.i == bcf_int16_missing)
                                kputc('.', s);
                            else
                                kputw(z->v1.i, s);
                            break;
                        case BCF_BT_INT32:
                            if(z->v1.i == bcf_int32_missing)
                                kputc('.', s);
                            else
                                kputw(z->v1.i, s);
                            break;
                        case BCF_BT_INT64:
                            if(z->v1.i == bcf_int64_missing)
                                kputc('.', s);
                            else
                                kputll(z->v1.i, s);
                            break;
                        case BCF_BT_FLOAT:
                            if(bcf_float_is_missing(z->v1.f))
                                kputc('.', s);
                            else
                                kputd(z->v1.f, s);
                            break;
                        case BCF_BT_CHAR:
                            kputc(z->v1.i, s);
                            break;
                        default:
                            throw std::runtime_error("Unexpected type in INFO");
                    }
                }
                else
                    bcf_fmt_array(s, z->len, z->type, z->vptr);
            }
            if(first) kputc('.', s);
        }
        else
            kputc('.', s);
        std::string out = std::string(s->s, s->l);
        free(s->s);
        free(s);
        return out;
    }

    /** @brief return boolean value indicates if genotypes of all samples are phased */
    inline bool allPhased() const
    {
        return isAllPhased;
    }

    /** @brief return the number of ploidy of current variant */
    inline int ploidy() const
    {
        return nploidy;
    }

    /** @brief in a rare case, one may want to set the number of ploidy manually */
    inline void setPloidy(int v)
    {
        nploidy = v;
    }

    /**
     * @brief vector of nsamples length. keep track of the type of genotype (one of GT_HOM_RR, GT_HET_RA,
     *        GT_HOM_AA, GT_HET_AA, GT_HAPL_R, GT_HAPL_A or GT_UNKN).
     * @note GT_HOM_RR 0 \n
     *       GT_HOM_AA 1 \n
     *       GT_HET_RA 2 \n
     *       GT_HET_AA 3 \n
     *       GT_HAPL_R 4 \n
     *       GT_HAPL_A 5 \n
     *       GT_UNKN   6 \n
     * */
    std::vector<char> typeOfGT;

    /** @brief vector of nsamples length. keep track of the phasing status of each sample */
    std::vector<char> gtPhase;
};

/**
 * @class BcfReader
 * @brief Stream in variants from compressed/uncompressed VCF/BCF file or stdin
 **/
class BcfReader
{
  private:
    std::shared_ptr<htsFile> fp; // hts file
    hts_idx_t * hidx = NULL; // hts index file
    tbx_t * tidx = NULL; // .tbi .csi index file for vcf files
    hts_itr_t * itr = NULL; // hts records iterator
    kstring_t s = {0, 0, NULL}; // kstring
    std::string fname;
    bool isBcf; // if the input file is bcf or vcf;

  public:
    /// a BcfHeader object
    BcfHeader header;
    /// number of samples in the VCF
    int nsamples;
    /// number of samples in the VCF
    std::vector<std::string> SamplesName;

    /// Construct an empty BcfReader
    BcfReader() {}

    /**
     *  @brief construct a vcf/bcf reader from file.
     *  @param file   the input vcf/bcf with suffix vcf(.gz)/bcf(.gz) or stdin "-"
     */
    BcfReader(const std::string & file) : fname(file)
    {
        open(file);
    }

    /**
     *  @brief construct a vcf/bcf reader with subset samples
     *  @param file   the input vcf/bcf with suffix vcf(.gz)/bcf(.gz) or stdin "-"
     *  @param region samtools-like region "chr:start-end", skip if empty
     */
    BcfReader(const std::string & file, const std::string & region) : fname(file)
    {
        open(file);
        if(!region.empty()) setRegion(region);
        SamplesName = header.getSamples();
    }

    /**
     *  @brief construct a vcf/bcf reader with subset samples in target region
     *  @param file   the input vcf/bcf with suffix vcf(.gz) or bcf(.gz)
     *  @param region samtools-like region "chr:start-end", skip if empty
     *  @param samples  LIST samples to include or exclude as a comma-separated string. \n
     *                  LIST : select samples in list \n
     *                  ^LIST : exclude samples from list \n
     *                  "-" : include all samples \n
     *                  "" : exclude all samples
     */
    BcfReader(const std::string & file, const std::string & region, const std::string & samples) : fname(file)
    {
        open(file);
        if(!region.empty()) setRegion(region);
        if(!samples.empty()) setSamples(samples);
    }

    /// return a BcfHeader object
    const BcfHeader & getHeader() const
    {
        return header;
    }

    /// open a VCF/BCF/STDIN file for streaming in
    void open(const std::string & file)
    {
        fname = file;
        fp = std::shared_ptr<htsFile>(hts_open(file.c_str(), "r"), htsFile_close());
        header.hdr = bcf_hdr_read(fp.get());
        nsamples = bcf_hdr_nsamples(header.hdr);
        SamplesName = header.getSamples();
    }

    /// return if the file is opened successfully
    bool isOpen() const
    {
        if(fp != NULL)
            return true;
        else
            return false;
    }

    /// close the BcfReader object.
    void close()
    {
        free(s.s);
        if(fp) fp.reset();
        if(itr) hts_itr_destroy(itr);
    }

    ~BcfReader() {}

    /** @brief set the number of threads to use */
    inline int setThreads(int n)
    {
        return hts_set_threads(fp.get(), n);
    }

    /** @brief get the number of records of given region */
    uint64_t getVariantsCount(BcfRecord & r, const std::string & region)
    {
        uint64_t c{0};
        while(getNextVariant(r)) c++;
        setRegion(region); // reset the region
        return c;
    }

    /**
     * @brief explicitly stream to specific samples
     * @param samples the string is bcftools-like format, which is comma separated list of samples to include
     * (or exclude with "^" prefix).
     * */
    void setSamples(const std::string & samples)
    {
        header.setSamples(samples);
        nsamples = bcf_hdr_nsamples(header.hdr);
        SamplesName = header.getSamples();
    }

    /**
     * @brief explicitly stream to specific region
     * @param region the string is samtools-like format which is chr:start-end
     * */
    void setRegion(const std::string & region)
    {
        // 1. check and load index first
        // 2. query iterval region
        // 3. if region is empty, use "."
        if(isEndWith(fname, "bcf") || isEndWith(fname, "bcf.gz"))
        {
            isBcf = true;
            hidx = bcf_index_load(fname.c_str());
            if(itr) bcf_itr_destroy(itr); // reset current region.
            if(region.empty())
                itr = bcf_itr_querys(hidx, header.hdr, ".");
            else
                itr = bcf_itr_querys(hidx, header.hdr, region.c_str());
        }
        else
        {
            isBcf = false;
            tidx = tbx_index_load(fname.c_str());
            assert(tidx != NULL && "error loading tabix index!");
            if(itr) tbx_itr_destroy(itr); // reset current region.
            if(region.empty())
                itr = tbx_itr_querys(tidx, ".");
            else
                itr = tbx_itr_querys(tidx, region.c_str());
            assert(itr != NULL && "no interval region found.failed!");
        }
    }

    /** @brief read in the next variant
     *  @param r r is a BcfRecord object to be filled in. */
    bool getNextVariant(BcfRecord & r)
    {
        int ret = -1;
        if(itr != NULL)
        {
            if(isBcf)
            {
                ret = bcf_itr_next(fp.get(), itr, r.line.get());
                bcf_unpack(r.line.get(), BCF_UN_ALL);
                return (ret >= 0);
            }
            else
            {
                int slen = tbx_itr_next(fp.get(), tidx, itr, &s);
                if(slen > 0)
                {
                    ret = vcf_parse1(&s, r.header->hdr, r.line.get()); // ret > 0, error
                    bcf_unpack(r.line.get(), BCF_UN_ALL);
                }
                return (ret <= 0) && (slen > 0);
            }
        }
        else
        {
            ret = bcf_read(fp.get(), r.header->hdr, r.line.get());
            // unpack record immediately. not lazy
            bcf_unpack(r.line.get(), BCF_UN_ALL);
            return (ret == 0);
        }
    }
};

/**
 * @class BcfWriter
 * @brief Stream out variants to compressed/uncompressed VCF/BCF file or stdout
 **/
class BcfWriter
{
  private:
    std::shared_ptr<htsFile> fp; // hts file
    std::shared_ptr<bcf1_t> b = std::shared_ptr<bcf1_t>(bcf_init(), variant_close()); // variant
    int ret;
    bool isHeaderWritten = false;
    const BcfHeader * hp;

  public:
    /// header object initialized by initalHeader
    BcfHeader header;

    /// Construct an empty BcfWriter
    BcfWriter() {}

    /**
     * @brief          Open VCF/BCF file for writing. The format is infered from file's suffix
     * @param fname    The file name or "-" for stdin/stdout. For indexed files
     * @param version  The output header is constructed with the internal template given a specific version
     */
    BcfWriter(const std::string & fname, std::string version = "VCF4.1")
    {
        open(fname);
        initalHeader(version);
        initalHeader(header);
    }

    /**
     * @brief          Open VCF/BCF file for writing. The format is infered from file's suffix
     * @param fname    The file name or "-" for stdin/stdout. For indexed files
     * @param h        The output header is pointing to the given BcfHeader object
     */
    BcfWriter(const std::string & fname, const BcfHeader & h)
    {
        open(fname);
        initalHeader(h);
    }

    /**
     * @brief          Open VCF/BCF file for writing using given mode
     * @param fname    The file name or "-" for stdin/stdout. For indexed files
     * @param version  The output header is constructed with the internal template given a specific version
     * @param mode     Mode matching \n
     *                 [w]b  .. compressed BCF \n
     *                 [w]bu .. uncompressed BCF \n
     *                 [w]z  .. compressed VCF \n
     *                 [w]   .. uncompressed VCF
     */
    BcfWriter(const std::string & fname, const std::string & version, const std::string & mode)
    {
        open(fname, mode);
        initalHeader(version);
        initalHeader(header);
    }

    /**
     * @brief          Open VCF/BCF file for writing using given mode
     * @param fname    The file name or "-" for stdin/stdout. For indexed files
     * @param h        The output header is pointing to the given BcfHeader object
     * @param mode     Mode matching \n
     *                 [w]b  .. compressed BCF \n
     *                 [w]bu .. uncompressed BCF \n
     *                 [w]z  .. compressed VCF \n
     *                 [w]   .. uncompressed VCF
     */
    BcfWriter(const std::string & fname, const BcfHeader & h, const std::string & mode)
    {
        open(fname, mode);
        initalHeader(h);
    }

    ~BcfWriter() {}

    /**
     * @brief          Open VCF/BCF file for writing. The format is infered from file's suffix
     * @param fname    The file name or "-" for stdin/stdout. For indexed files
     */
    void open(const std::string & fname)
    {
        std::string mode{"w"};
        if(isEndWith(fname, "bcf.gz")) mode += "b";
        if(isEndWith(fname, "bcf")) mode += "bu";
        if(isEndWith(fname, "vcf.gz")) mode += "z";
        fp = std::shared_ptr<htsFile>(hts_open(fname.c_str(), mode.c_str()), htsFile_close());
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
    void open(const std::string & fname, const std::string & mode)
    {
        fp = std::shared_ptr<htsFile>(hts_open(fname.c_str(), mode.c_str()), htsFile_close());
    }

    /// close the BcfWriter object.
    void close()
    {
        if(!isHeaderWritten) writeHeader();
        if(b) b.reset();
        if(fp) fp.reset();
    }

    /// initial a VCF header using the internal template given a specific version. VCF4.1 is the default
    void initalHeader(std::string version = "VCF4.1")
    {
        header.hdr = bcf_hdr_init("w");
        header.setVersion(version);
    }

    /// initial a VCF header by pointing to header of another VCF
    void initalHeader(const BcfHeader & h)
    {
        hp = &h;
    }

    /// copy a string to a vcf line
    void writeLine(const std::string & vcfline)
    {
        if(!isHeaderWritten && !writeHeader()) throw std::runtime_error("could not write header\n");
        kstring_t s = {0, 0, NULL};
        kputsn(vcfline.c_str(), vcfline.length(), &s);
        ret = vcf_parse1(&s, hp->hdr, b.get());
        free(s.s);
        if(ret > 0) throw std::runtime_error("error parsing: " + vcfline + "\n");
        if(b->errcode == BCF_ERR_CTG_UNDEF)
        {
            throw std::runtime_error("contig id " + (std::string)bcf_hdr_id2name(hp->hdr, b->rid)
                                     + " not found in the header. please run header->AddContig() first.\n");
        }
        ret = bcf_write(fp.get(), hp->hdr, b.get());
        if(ret != 0) throw std::runtime_error("error writing: " + vcfline + "\n");
    }

    /// streams out the header
    bool writeHeader()
    {
        ret = bcf_hdr_write(fp.get(), hp->hdr);
        if(ret == 0)
            return isHeaderWritten = true;
        else
            return false;
    }

    /// streams out the given variant of BcfRecord type
    inline bool writeRecord(BcfRecord & v)
    {
        if(!isHeaderWritten) writeHeader();
        if(bcf_write(fp.get(), v.header->hdr, v.line.get()) < 0)
            return false;
        else
            return true;
    }
};

} // namespace vcfpp

#endif // VCFPP_H_
