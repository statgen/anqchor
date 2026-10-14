#ifndef __QPGEN_UTILS_H
#define __QPGEN_UTILS_H

#include <vector>
#include <string>

#include "Eigen/Dense"
#include "qgenlib/qgen_error.h"
#include "qgenlib/hts_utils.h"
#include "qpgen.h"

// identify overlapping IDS from list1 and list2, and store sorted overlapping IDs in overlapping_ids
// Field delimiter of user-supplied tables (phenotype, covariate and PRS matrices, pair files,
// sample, trait and variant lists): a tab by default, or a single space with --space-delim.
// The two are never mixed, so with the default a field may contain spaces.
#define QPGEN_SPACE_DELIM_HELP "Input tables (phenotype/covariate matrices, pair files, sample/trait/variant lists) are delimited by single spaces instead of tabs (default: tab-delimited)"
inline char table_delimiter(bool space_delim) { return space_delim ? ' ' : '\t'; }

// --cat-cov-names / --max-cat-levels, shared by pair-assoc, rect-assoc and region-assoc
#define QPGEN_CAT_COV_NAMES_HELP "Comma-separated covariates (column names in --cov) to treat as categorical, as REGENIE's --catCovarList. Their values may be any strings; each covariate with L levels among the analyzed samples is replaced by L-1 indicator columns (the first level in sorted order is the reference)"
#define QPGEN_MAX_CAT_LEVELS_HELP "Maximum number of levels of a categorical covariate, as REGENIE's --maxCatLevels (default: 10)"
#define QPGEN_COV_NAMES_HELP "Comma-separated covariates (column names in --cov) to use, as REGENIE's --covarColList; ranges such as PC{1:10} are expanded. Covariates in --cat-cov-names are always used. Default: all columns"
// splits a comma-separated list, dropping empty items
inline std::vector<std::string> split_comma_list(const std::string& s) {
    std::vector<std::string> out;
    size_t beg = 0;
    while ( beg <= s.size() ) {
        size_t end = s.find(',', beg);
        if ( end == std::string::npos ) end = s.size();
        if ( end > beg ) out.push_back(s.substr(beg, end - beg));
        beg = end + 1;
    }
    return out;
}

// split_comma_list() with REGENIE-style ranges expanded: "PC{1:3}" -> PC1, PC2, PC3
inline std::vector<std::string> expand_name_ranges(const std::string& s) {
    std::vector<std::string> out;
    for(const std::string& item : split_comma_list(s)) {
        const size_t lb = item.find('{'), colon = item.find(':', lb), rb = item.find('}', lb);
        if ( lb == std::string::npos || colon == std::string::npos || rb == std::string::npos || !( lb < colon && colon < rb ) ) {
            out.push_back(item);
            continue;
        }
        const int beg = atoi(item.substr(lb + 1, colon - lb - 1).c_str());
        const int end = atoi(item.substr(colon + 1, rb - colon - 1).c_str());
        if ( end < beg ) error("Invalid range %s: the end is smaller than the start", item.c_str());
        for(int i = beg; i <= end; ++i) out.push_back(item.substr(0, lb) + std::to_string(i) + item.substr(rb + 1));
    }
    return out;
}

int32_t identify_overlapping_ids(const std::vector<std::string>& list1, const std::vector<std::string>& list2, std::vector<std::string>& overlapping_ids);

// index overlapping IDs from list1 and list2, and store sorted indices in index1 and index2
int32_t index_overlapping_ids(const std::vector<std::string>& list1, const std::vector<std::string>& list2, std::vector<int32_t>& index1, std::vector<int32_t>& index2);

void fill_list_indices(const std::vector<std::string>& list, const std::map<std::string, int32_t>& list2idx, std::vector<int32_t>& indices);

int32_t rebuild_id2index_map(
    const std::vector<std::string>& ids,
    std::map<std::string, int32_t>& id2idx_map
);

#endif // __QPGEN_UTILS_H