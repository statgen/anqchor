#ifndef __POS_LOCI_H
#define __POS_LOCI_H

// this is very similar to genomeLocus, except that chrom is not a part of the key

#include <vector>
#include <cstdint>
#include <set>
#include <map>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <climits>
#include <algorithm>

#include "qgen_error.h"
#include "hts_utils.h"

// A single genomic int32_terval
class posLocus
{
public:
    int32_t beg1; // includes 1-based, excludes 0-based
    int32_t end0; // excludes 1-based, includes 0-based

    posLocus(int32_t b, int32_t e) : beg1(b), end0(e) {}

    // compare between genomeLocus
    // l : rhs argument of the same type
    // returns true iff this < l
    inline bool operator<(const posLocus &l) const
    {
        if (beg1 == l.beg1)
        {
            return end0 < l.end0;
        }
        else
            return beg1 < l.beg1;
    }

    // compare between genomeLocus
    // l : rhs argument of the same type
    // returns true iff this == l
    inline bool operator==(const posLocus &l) const
    {
        return ((beg1 == l.beg1) && (end0 == l.end0));
    }

    // returns length of the interval
    int32_t length() const { return end0 - beg1 + 1; }

    // returns the total overlapping intervals
    int32_t overlapBases(int32_t _beg1, int32_t _end0) const
    {
        if ((beg1 <= end0) && (_beg1 <= end0))
        {
            return (_beg1 < end0 ? _beg1 : end0) - (beg1 < end0 ? end0 : beg1) + 1;
        }
        else
            return 0;
    }

    int32_t overlapBases(const posLocus &l) const
    {
        return overlapBases(l.beg1, l.end0);
    }

    bool overlaps(int32_t _beg1, int32_t _end0) const
    {
        if ((beg1 <= _end0) && (_beg1 <= end0))
            return true;
        else
            return false;
    }

    // check overlap with other locus
    bool overlaps(const posLocus &l) const
    {
        return overlaps(l.beg1, l.end0);
    }

    // merge two locus if possible
    bool merge(const posLocus &l)
    {
        if ((beg1 - 1 <= l.end0) && (l.beg1 - 1 <= end0))
        {
            if (l.beg1 < beg1)
                beg1 = l.beg1;
            if (l.end0 > end0)
                end0 = l.end0;
            return true;
        }
        else
        {
            return false;
        }
    }

    // check whether the interval contains a particular position in 1-based coordinate
    bool contains1(int32_t pos1 = INT_MAX) const
    {
        return ((pos1 >= beg1) && (pos1 <= end0));
    }

    // check whether the interval contains a particular position in 0-based coordinate
    bool contains0(int32_t pos0) const { return contains1(pos0 + 1); }

    //   // parse a string in [chr]:[beg1]-[end0] format
    //   static bool parseRegion(const char* region, std::string& chrom, int32_t& beg1, int32_t& end0) {
    //     char buf[255];
    //     strcpy(buf,region);
    //     const char* pcolon = strchr(region,':');
    //     const char* pminus = strchr(pcolon+1,'-');
    //     if ( pcolon == NULL )
    //       error("Cannot parse %s in posLocus::posLocus()");
    //     chrom.assign(region,0,pcolon-region);
    //     beg1 = atoi(pcolon+1);
    //     if ( pminus == NULL ) end0 = INT_MAX;
    //     else {
    //       end0 = atoi(pminus+1);
    //       if ( end0 == 0 ) end0 = INT_MAX;
    //     }
    //     return true;
    //   }

    //   // parse a string in [chr]:[beg1]-[end0] format
    //   static bool parseBegLenStrand(const char* region, std::string& chrom, int32_t& beg1, int32_t& end0, bool& fwdStrand) {
    //     char buf[255];
    //     strcpy(buf,region);
    //     const char* pcolon1 = strchr(region,':');
    //     const char* pcolon2 = strchr(pcolon1+1,':');
    //     const char* pcolon3 = strchr(pcolon2+1,':');
    //     if ( pcolon1 == NULL )
    //       error("Cannot parse %s in genomeLocus::genomeLocus()");
    //     chrom.assign(region,0,pcolon1-region);
    //     beg1 = atoi(pcolon1+1);
    //     end0 = beg1 + atoi(pcolon2+1) - 1;
    //     fwdStrand = ((pcolon3 != NULL) && (pcolon3[1] == '-')) ? false : true;
    //     return true;
    //   }
};

// Collection of genomic locus
class posLoci
{
public:
    std::set<posLocus> loci;
    std::set<posLocus>::iterator it;
    bool overlapResolved;
    int32_t maxLength;

    posLoci() : overlapResolved(false), maxLength(0) {}

    posLoci(int32_t beg1, int32_t end0) : overlapResolved(false), maxLength(0)
    {
        add(beg1, end0);
        resolveOverlaps();
    }

    // posLoci(const char *reg) : overlapResolved(false), maxLength(0)
    // {
    //     add(reg);
    //     resolveOverlaps();
    // }

    // functions for iterating each locus
    inline void rewind() { it = loci.begin(); }
    inline bool next()
    {
        ++it;
        return (it != loci.end());
    }
    inline bool isend() { return (it == loci.end()); }
    inline const posLocus &currentLocus() { return (*it); }

    // check the size
    bool empty() { return loci.empty(); }

    bool clear()
    {
        if (loci.empty())
            return false;

        loci.clear();
        it = loci.begin();
        overlapResolved = false;
        maxLength = 0;
        return true;
    }

    int32_t numLocus() const { return (int32_t)loci.size(); }

    // add a locus
    bool add(int32_t beg1, int32_t end0)
    {
        overlapResolved = false;
        if (end0 - beg1 + 1 > maxLength)
            maxLength = end0 - beg1 + 1;
        std::pair<std::set<posLocus>::iterator, bool> ret = loci.insert(posLocus(beg1, end0));
        it = ret.first;
        return ret.second;
    }

    // Resolve overlapping int32_tervals
    int32_t resolveOverlaps()
    {
        if (!overlapResolved)
        {
            std::set<posLocus>::iterator it;
            std::set<posLocus>::iterator prev;
            int32_t numMerged = 0;
            for (it = loci.begin(); it != loci.end(); ++it)
            {
                if (it != loci.begin())
                {
                    if (prev->overlaps(*it))
                    {
                        // if overlaps, erase both and insert merged one
                        posLocus locus = *prev;
                        locus.merge(*it);
                        if ((int32_t)locus.length() > maxLength)
                            maxLength = locus.length();
                        loci.erase(it);
                        loci.erase(prev);
                        prev = it = loci.insert(locus).first;
                        ++numMerged;
                    }
                    else
                    {
                        prev = it;
                    }
                }
                else
                {
                    prev = it;
                }
            }
            overlapResolved = true;
            return numMerged;
        }
        else
        {
            return 0;
        }
        return 0;
    }

    unsigned long totalLength() const
    {
        // resolveOverlaps();
        unsigned long sz = 0;
        std::set<posLocus>::iterator it2;
        for (it2 = loci.begin(); it2 != loci.end(); ++it2)
        {
            sz += it2->length();
        }
        return sz;
    }

    bool moveTo(int32_t pos1 = INT_MAX)
    {
        if (loci.empty())
            return false;
        if ((it != loci.end()) && it->contains1(pos1))
            return true;

        posLocus locus(pos1, pos1);
        it = loci.lower_bound(locus);
        if (it == loci.begin())
        { // do nothing
            // notice("beg");
            return (it->contains1(pos1));
        }
        else if (it == loci.end())
        {
            // notice("end");
            std::set<posLocus>::iterator i = it;
            --i;
            if (i->contains1(pos1))
            {
                it = i;
                return true;
            }
            else
            {
                rewind();
                return false;
            }
        }
        else
        {
            // notice("mid");
            if (it->contains1(pos1))
                return true;
            else
            {
                std::set<posLocus>::iterator i = it;
                --i;
                if (i->contains1(pos1))
                {
                    it = i;
                    return true;
                }
                else
                {
                    rewind();
                    return false;
                }
            }
        }
    }

    bool contains1(int32_t pos1)
    {
        if (loci.empty())
            return false;
        posLocus locus(pos1, pos1);
        std::set<posLocus>::iterator it2 = loci.lower_bound(locus);
        if (it2 != loci.begin())
            --it2;
        while (it2 != loci.end() && (it2->beg1 <= pos1))
        {
            if (it2->end0 >= pos1)
                return true;
            ++it2;
        }
        return false;
    }

    bool overlaps(int32_t beg1, int32_t end0)
    {
        if (loci.empty())
            return false;

        posLocus locus(overlapResolved ? beg1 : beg1 - maxLength, overlapResolved ? beg1 : beg1 - maxLength);
        if (loci.empty())
            return false;
        std::set<posLocus>::iterator it2 = loci.lower_bound(locus);
        if (it2 != loci.begin())
            --it2;
        while (it2 != loci.end() && (it2->beg1 <= end0))
        {
            if ((it2->beg1 <= end0) && (beg1 <= it2->end0))
                return true;
            ++it2;
        }
        return false;
    }

    bool contains(int32_t beg1, int32_t end0)
    {
        if (loci.empty())
            return false;

        resolveOverlaps();
        posLocus locus(beg1 - maxLength, beg1 - maxLength);
        std::set<posLocus>::iterator it2 = loci.lower_bound(locus);
        if (it2 != loci.begin())
            --it2;
        while (it2 != loci.end() && (it2->beg1 <= end0))
        {
            if ((it2->beg1 <= beg1) && (end0 <= it2->end0))
                return true;
            ++it2;
        }
        return false;
    }
};

// Collection of position loci
template <class T>
class posLocusMap
{
public:
    std::map<posLocus, T> loci;
    typename std::map<posLocus, T>::iterator it;
    int32_t maxLength;

    posLocusMap() : maxLength(0) { it = loci.end(); }
    posLocusMap(int32_t beg1, int32_t end0, const T &val) : maxLength(0)
    {
        add(beg1, end0, val);
        it = loci.end();
    }

    // functions for iterating each locus
    void rewind() { it = loci.begin(); }
    bool next()
    {
        ++it;
        return (it != loci.end());
    }
    bool isend() { return (it == loci.end()); }
    const posLocus &currentLocus() { return (*it); }

    // check the size
    bool empty() { return loci.empty(); }

    bool clear()
    {
        if (loci.empty())
            return false;
        loci.clear();
        it = loci.begin();
        maxLength = 0;
        return true;
    }

    int32_t numLocus() const { return (int32_t)loci.size(); }

    // add a locus
    bool add(int32_t beg1, int32_t end0, const T &val)
    {
        if (end0 - beg1 + 1 > maxLength)
            maxLength = end0 - beg1 + 1;
        std::pair<typename std::map<posLocus, T>::iterator, bool> ret =
            loci.insert(std::pair<posLocus, T>(posLocus(beg1, end0), val));
        it = ret.first;
        return ret.second;
    }

    unsigned long totalLength() const
    {
        unsigned long sz = 0;
        typename std::map<posLocus, T>::iterator it2;
        for (it2 = loci.begin(); it2 != loci.end(); ++it2)
        {
            sz += it2->first.length();
        }
        return sz;
    }

    bool moveTo(int32_t pos1 = INT_MAX)
    {
        if (loci.empty())
            return false;

        if ((it != loci.end()) && it->first.contains1(pos1))
            return true;

        posLocus locus(pos1, pos1);
        it = loci.lower_bound(locus);
        if (it == loci.begin())
        { // do nothing
            return (it->first.contains1(pos1));
        }
        else if (it == loci.end())
        {
            typename std::map<posLocus, T>::iterator i = it;
            --i;
            if (i->first.contains1(pos1))
            {
                it = i;
                return true;
            }
            else
            {
                rewind();
                return false;
            }
        }
        else
        {
            if (it->first.contains1(pos1))
                return true;
            else
            {
                typename std::map<posLocus, T>::iterator i = it;
                --i;
                if (i->first.contains1(pos1))
                {
                    it = i;
                    return true;
                }
                else
                {
                    rewind();
                    return false;
                }
            }
        }
    }

    bool contains1(int32_t pos1)
    {
        posLocus locus(pos1, pos1);
        typename std::map<posLocus, T>::iterator it2 = loci.lower_bound(locus);
        if (it2 != loci.begin())
            --it2;
        while (it2 != loci.end() && (it2->first.beg1 <= pos1))
        {
            if (it2->first.end0 >= pos1)
                return true;
            ++it2;
        }
        return false;
    }

    bool overlaps(int32_t beg1, int32_t end0)
    {
        posLocus locus(beg1 - maxLength, beg1 - maxLength);
        if (loci.empty())
            return false;
        typename std::map<posLocus, T>::iterator it2 = loci.lower_bound(locus);
        if (it2 != loci.begin())
            --it2;
        while (it2 != loci.end() && (it2->first.beg1 <= end0))
        {
            if ((it2->first.beg1 <= end0) && (beg1 <= it2->first.end0))
                return true;
            ++it2;
        }
        return false;
    }

    bool contains(int32_t beg1, int32_t end0)
    {
        if (loci.empty())
            return false;

        posLocus locus(beg1 - maxLength, beg1 - maxLength);
        typename std::map<posLocus, T>::iterator it2 = loci.lower_bound(locus);
        if (it2 != loci.begin())
            --it2;
        while (it2 != loci.end() && (it2->first.beg1 <= end0))
        {
            if ((it2->first.beg1 <= beg1) && (end0 <= it2->first.end0))
                return true;
            ++it2;
        }
        return false;
    }
};

#endif // POS_LOCI_H
