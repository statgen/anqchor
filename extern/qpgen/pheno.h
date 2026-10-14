#ifndef __ASSOC_PHENO_H
#define __ASSOC_PHENO_H

#include <set>
#include <string>

#include "qgenlib/qgen_error.h"
#include "qgenlib/hts_utils.h"
#include "qgenlib/genome_loci.h"
#include "qpgen.h"
#include "Eigen/Dense"

class PhenoMatrix {
public:
    Eigen::MatrixXd pheno_mat;
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> pheno_mask;
    std::vector<std::string> pheno_ids;
    std::vector<genomeLocus> pheno_loci;
    genomeLocusMap<int32_t> pheno_locusmap;
    std::vector<std::string> samp_ids;
    std::string format;
    bool has_missing;
    bool has_loci;
    std::set<std::string> missing_strs; 
    std::map<std::string, int32_t> pheno_id2idx;
    std::map<std::string, int32_t> samp_id2idx;
    // Categorical covariates (--cat-cov-names): their columns hold text levels in cat_strs (one
    // string per sample, "" when missing) and 0 in pheno_mat until expand_categorical() turns
    // them into indicator columns. cat_strs is empty when no column is categorical; otherwise it
    // has one entry per column, empty for numeric columns.
    std::set<std::string> cat_ids;
    std::vector<std::vector<std::string> > cat_strs;
    // --cov-names: when non-empty, other (non-categorical) columns are read as text without being
    // parsed, and keep_covariate_columns() then drops them
    std::set<std::string> selected_ids;
    bool is_text_column(const std::string& id) const {
        return cat_ids.count(id) > 0 || ( !selected_ids.empty() && selected_ids.count(id) == 0 );
    }

    PhenoMatrix() : has_missing(false), has_loci(false) {}
    ~PhenoMatrix() {}
    bool is_missing(const std::string& str) { return missing_strs.find(str) != missing_strs.end(); }
    bool is_missing(const char* str) { return is_missing(std::string(str)); }
    bool add_missing_str(const std::string& str) { return missing_strs.insert(str).second; }
    void reset_missing_strs() { missing_strs.clear(); }
    // register comma-separated missing strings (e.g. "NA,NaN,."); returns the number of new strings added
    int32_t add_missing_strs(const std::string& comma_separated_strs) {
        int32_t n_added = 0;
        size_t beg = 0;
        while ( beg <= comma_separated_strs.size() ) {
            size_t end = comma_separated_strs.find(',', beg);
            if ( end == std::string::npos ) end = comma_separated_strs.size();
            if ( end > beg && add_missing_str(comma_separated_strs.substr(beg, end - beg)) ) ++n_added;
            beg = end + 1;
        }
        return n_added;
    }
    // refresh has_missing from the current mask (e.g. after subsetting samples or phenotypes)
    void recompute_has_missing() { has_missing = ( pheno_mask.size() > 0 ) && !pheno_mask.all(); }
    // impute missing cells column-wise with the mean (as_half_min=false) or half of the minimum (as_half_min=true)
    // of the observed values, mark them observed, and refresh has_missing. Columns with no observed value are left
    // missing. Returns the number of imputed cells.
    int64_t impute_missing(bool as_half_min);
    
    // Keeps only the listed columns and the categorical ones, in file order (--cov-names). Aborts on
    // a listed name that is not a column of file. No-op for an empty list.
    void keep_covariate_columns(const std::vector<std::string>& cov_names, const char* file);
    // the --cov-names columns, so that the other columns are not parsed (call before load_pheno_matrix)
    void set_selected_columns(const std::vector<std::string>& ids) { selected_ids = std::set<std::string>(ids.begin(), ids.end()); }
    // Columns to read as categorical (call before load_pheno_matrix); each must be in the file
    void set_categorical_ids(const std::vector<std::string>& ids) { cat_ids = std::set<std::string>(ids.begin(), ids.end()); }
    bool has_categorical() const { for(const std::vector<std::string>& v : cat_strs) if ( !v.empty() ) return true; return false; }
    // Replaces each categorical column by L-1 indicator columns <id>_<level> for its L observed
    // levels (sorted; the first is the reference), on the current samples. A missing value stays
    // missing in every indicator. Aborts when a column has more than max_levels levels. Returns
    // the number of indicator columns.
    int32_t expand_categorical(int32_t max_levels);

    bool load_pheno_matrix(const char* pheno_file, const char* pheno_format, const char delim = '\t'); // delim: see table_delimiter()
    int32_t subset_sample_ids(const std::vector<std::string>& samp_ids);
    int32_t subset_pheno_ids(const std::vector<std::string>& pheno_ids);
    bool subset_sample_pheno_ids(const std::vector<std::string>& samp_ids, const std::vector<std::string>& pheno_ids);
    int32_t subset_sample_indices(const std::vector<int32_t>& samp_indices);
    int32_t subset_pheno_indices(const std::vector<int32_t>& pheno_indices);
    bool subset_sample_pheno_indices(const std::vector<int32_t>& samp_indices, const std::vector<int32_t>& pheno_indices);

    bool sample_ids_sorted() const;
    bool pheno_ids_sorted() const;
};

#endif // __ASSOC_PHENO_H