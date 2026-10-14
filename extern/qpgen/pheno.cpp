#include "pheno.h"
#include <strings.h>
#include "qgenlib/qgen_error.h"
#include "qgenlib/tsv_reader.h"
#include "qpgen_utils.h"

#include <cstring>
#include <cstdlib>
#include <cctype>
#include <cmath>

namespace {

// Parses one phenotype/covariate value. A value that is not a finite number, and is not a
// registered missing string, aborts instead of being read as 0 (atof() would silently turn
// "NA", "XX" or "male" into 0).
double parse_pheno_value(const char* s, const char* pheno_file, int32_t line_no, int32_t field_no)
{
    char* end = NULL;
    const double v = strtod(s, &end);
    while ( end != NULL && *end != '\0' && isspace((unsigned char)*end) ) ++end; // allow trailing '\r'
    if ( end == s || end == NULL || *end != '\0' || !std::isfinite(v) ) {
        error("Value '%s' at line %d, field %d of %s is not a finite number. Code categorical values (e.g. sex) as numbers, and mark missing values with a string registered as missing (--missing-str, where the command supports it)",
              s, line_no, field_no, pheno_file);
    }
    return v;
}

// Aborts unless the header has a field beyond the index columns. A header that splits into too
// few fields almost always means the file uses another delimiter than the one in effect.
void check_header_fields(int32_t nfields, int32_t icol_data, char delim, const char* pheno_file)
{
    if ( nfields <= icol_data ) {
        const char* how = ( delim == '\t' ) ? "tabs" : ( ( delim == ' ' ) ? "single spaces" : "whitespace" );
        error("The header of %s has %d field(s) when split on %s, but at least %d are expected. Input tables are tab-delimited unless --space-delim is given; check the file's delimiter",
              pheno_file, nfields, how, icol_data + 1);
    }
}

std::vector<std::string> pick_strs(const std::vector<std::string>& v, const std::vector<int32_t>& idx) {
    std::vector<std::string> out;
    out.reserve(idx.size());
    for(int32_t i : idx) out.push_back(v[i]);
    return out;
}

} // anonymous namespace

bool PhenoMatrix::subset_sample_pheno_indices(const std::vector<int32_t>& samp_indices, const std::vector<int32_t>& pheno_indices) {
    if ( samp_indices.size() == this->samp_ids.size() && pheno_indices.size() == this->pheno_ids.size() ) {
        bool change_needed = false;
        for ( int i = 0; i < samp_indices.size(); ++i ) {
            if ( samp_indices[i] != i ) {
                change_needed = true;
                break;
            }
        }
        if ( !change_needed ) {
            for ( int j = 0; j < pheno_indices.size(); ++j ) {
                if ( pheno_indices[j] != j ) {
                    change_needed = true;
                    break;
                }
            }
        }
        if ( !change_needed ) {
            return true;
        }
    }

    std::vector<std::string> samp_ids_sub;
    std::vector<std::string> pheno_ids_sub;
    Eigen::MatrixXd pheno_mat_sub(samp_indices.size(), pheno_indices.size());
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> pheno_mask_sub(samp_indices.size(), pheno_indices.size());
    std::vector<genomeLocus> pheno_loci_sub;
    for ( size_t i = 0; i < samp_indices.size(); ++i ) {
        int32_t sidx = samp_indices[i];
        if ( sidx < 0 || sidx >= (int32_t)samp_ids.size() ) {
            error("Sample index %d is out of bounds [0,%d)", sidx, (int32_t)samp_ids.size());
        }
        samp_ids_sub.push_back( samp_ids[sidx] );
        for ( size_t j = 0; j < pheno_indices.size(); ++j ) {
            int32_t pidx = pheno_indices[j];
            if ( pidx < 0 || pidx >= (int32_t)this->pheno_ids.size() ) {
                error("Phenotype index %d is out of bounds [0,%d)", pidx, (int32_t)this->pheno_ids.size());
            }
            if ( i == 0 ) {
                pheno_ids_sub.push_back( this->pheno_ids[pidx] );
                if ( has_loci ) {
                    pheno_loci_sub.push_back( this->pheno_loci[pidx] );
                }
            }
            pheno_mat_sub(i,j) = pheno_mat(sidx, pidx);
            pheno_mask_sub(i,j) = pheno_mask(sidx, pidx);
        }
    }
    if ( !cat_strs.empty() ) {
        std::vector<std::vector<std::string> > cat_sub(pheno_indices.size());
        for(size_t j = 0; j < pheno_indices.size(); ++j) {
            if ( !cat_strs[pheno_indices[j]].empty() ) cat_sub[j] = pick_strs(cat_strs[pheno_indices[j]], samp_indices);
        }
        cat_strs.swap(cat_sub);
    }
    pheno_mat = pheno_mat_sub;
    pheno_mask = pheno_mask_sub;
    this->samp_ids = samp_ids_sub;
    this->pheno_ids = pheno_ids_sub;
    if ( has_loci ) {
        this->pheno_loci = pheno_loci_sub; 
    }
    rebuild_id2index_map(samp_ids, samp_id2idx);
    rebuild_id2index_map(pheno_ids, pheno_id2idx);
    recompute_has_missing();
    return true;
}

int32_t PhenoMatrix::subset_pheno_indices(const std::vector<int32_t>& pheno_indices) {
    if ( pheno_indices.size() == this->pheno_ids.size() ) {
        bool change_needed = false;
        for ( int j = 0; j < pheno_indices.size(); ++j ) {
            if ( pheno_indices[j] != j ) {
                change_needed = true;
                break;
            }
        }
        if ( !change_needed ) {
            return (int32_t)this->pheno_ids.size();
        }
    }

    std::vector<std::string> pheno_ids_sub;
    Eigen::MatrixXd pheno_mat_sub(pheno_mat.rows(), pheno_indices.size());
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> pheno_mask_sub(pheno_mat.rows(), pheno_indices.size());
    std::vector<genomeLocus> pheno_loci_sub;
    for ( size_t j = 0; j < pheno_indices.size(); ++j ) {
        int32_t idx = pheno_indices[j];
        if ( idx < 0 || idx >= (int32_t)pheno_ids.size() ) {
            error("Phenotype index %d is out of bounds [0,%d)", idx, (int32_t)this->pheno_ids.size());
        }
        pheno_ids_sub.push_back( pheno_ids[idx] );
        if ( has_loci ) {
            pheno_loci_sub.push_back( pheno_loci[idx] );
        }
        pheno_mat_sub.col(j) = pheno_mat.col( idx );
        pheno_mask_sub.col(j) = pheno_mask.col( idx );
    }
    if ( !cat_strs.empty() ) {
        std::vector<std::vector<std::string> > cat_sub(pheno_indices.size());
        for(size_t j = 0; j < pheno_indices.size(); ++j) cat_sub[j] = cat_strs[pheno_indices[j]];
        cat_strs.swap(cat_sub);
    }
    pheno_mat = pheno_mat_sub;
    pheno_mask = pheno_mask_sub;
    this->pheno_ids = pheno_ids_sub;
    if ( has_loci ) {
        this->pheno_loci = pheno_loci_sub;
    }

    rebuild_id2index_map(pheno_ids, pheno_id2idx);
    recompute_has_missing();

    return (int32_t)this->pheno_ids.size();
}

int32_t PhenoMatrix::subset_sample_indices(const std::vector<int32_t>& samp_indices) {
    if ( samp_indices.size() == this->samp_ids.size() ) {
        bool change_needed = false;
        for ( int i = 0; i < samp_indices.size(); ++i ) {
            if ( samp_indices[i] != i ) {
                change_needed = true;
                break;
            }
        }
        if ( !change_needed ) {
            return (int32_t)this->samp_ids.size();
        }
    }

    std::vector<std::string> samp_ids_sub;
    Eigen::MatrixXd pheno_mat_sub(samp_indices.size(), pheno_mat.cols());
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> pheno_mask_sub(samp_indices.size(), pheno_mat.cols());
    for ( size_t i = 0; i < samp_indices.size(); ++i ) {
        int32_t idx = samp_indices[i];
        if ( idx < 0 || idx >= (int32_t)samp_ids.size() ) {
            error("Sample index %d is out of bounds [0,%d)", idx, (int32_t)samp_ids.size());
        }
        samp_ids_sub.push_back( samp_ids[idx] );
        pheno_mat_sub.row(i) = pheno_mat.row( idx );
        pheno_mask_sub.row(i) = pheno_mask.row( idx );
    }
    for(std::vector<std::string>& v : cat_strs) {
        if ( !v.empty() ) v = pick_strs(v, samp_indices);
    }
    pheno_mat = pheno_mat_sub;
    pheno_mask = pheno_mask_sub;
    this->samp_ids = samp_ids_sub;

    rebuild_id2index_map(samp_ids, samp_id2idx);
    recompute_has_missing();

    return (int32_t)this->samp_ids.size();
}

bool PhenoMatrix::subset_sample_pheno_ids(const std::vector<std::string>& samp_ids, const std::vector<std::string>& pheno_ids) {
    std::vector<int32_t> phe_sample_indices;
    std::vector<int32_t> phe_pheno_indices;
    for(int32_t i = 0; i < samp_ids.size(); ++i ) {
        const std::string& sid = samp_ids[i];
        std::map<std::string, int32_t>::const_iterator it = samp_id2idx.find( sid );
        if ( it != samp_id2idx.end() ) {
            phe_sample_indices.push_back( it->second );
        }
        else {
            error("Sample ID '%s' not found in phenotype matrix", sid.c_str());
        }
    }
    for(int32_t i = 0; i < pheno_ids.size(); ++i ) {
        const std::string& pid = pheno_ids[i];
        std::map<std::string, int32_t>::const_iterator it = pheno_id2idx.find( pid );
        if ( it != pheno_id2idx.end() ) {
            phe_pheno_indices.push_back( it->second );
        }
        else {
            error("Phenotype ID '%s' not found in phenotype matrix", pid.c_str());
        }
    }
    return subset_sample_pheno_indices(phe_sample_indices, phe_pheno_indices);
}


int32_t PhenoMatrix::subset_sample_ids(const std::vector<std::string>& samp_ids) {
    std::vector<int32_t> phe_sample_indices;
    for(int32_t i = 0; i < samp_ids.size(); ++i ) {
        const std::string& sid = samp_ids[i];
        std::map<std::string, int32_t>::const_iterator it = samp_id2idx.find( sid );
        if ( it != samp_id2idx.end() ) {
            phe_sample_indices.push_back( it->second );
        }
        else {
            error("Sample ID '%s' not found in phenotype matrix", sid.c_str());
        }
    }
    return subset_sample_indices( phe_sample_indices );
} 

int32_t PhenoMatrix::subset_pheno_ids(const std::vector<std::string>& pheno_ids) {
    std::vector<int32_t> phe_pheno_indices;
    for(int32_t i = 0; i < pheno_ids.size(); ++i ) {
        const std::string& pid = pheno_ids[i];
        std::map<std::string, int32_t>::const_iterator it = pheno_id2idx.find( pid );
        if ( it != pheno_id2idx.end() ) {
            phe_pheno_indices.push_back( it->second );
        }
        else {
            error("Phenotype ID '%s' not found in phenotype matrix", pid.c_str());
        }
    }
    return subset_pheno_indices( phe_pheno_indices );
}

int64_t PhenoMatrix::impute_missing(bool as_half_min) {
    int64_t n_imputed = 0;
    int32_t n_cols_all_missing = 0;
    int32_t n_cols_nonpositive_min = 0;
    for ( int32_t j = 0; j < pheno_mat.cols(); ++j ) {
        double sum = 0.0, minv = 0.0;
        int32_t n_obs = 0;
        for ( int32_t i = 0; i < pheno_mat.rows(); ++i ) {
            if ( !pheno_mask(i, j) ) continue;
            double v = pheno_mat(i, j);
            if ( n_obs == 0 || v < minv ) minv = v;
            sum += v;
            ++n_obs;
        }
        if ( n_obs == 0 ) { ++n_cols_all_missing; continue; }
        if ( n_obs == pheno_mat.rows() ) continue;
        if ( as_half_min && minv <= 0.0 ) ++n_cols_nonpositive_min;
        double fill = as_half_min ? 0.5 * minv : sum / (double)n_obs;
        for ( int32_t i = 0; i < pheno_mat.rows(); ++i ) {
            if ( !pheno_mask(i, j) ) {
                pheno_mat(i, j) = fill;
                pheno_mask(i, j) = true;
                ++n_imputed;
            }
        }
    }
    if ( n_cols_all_missing > 0 ) {
        warning("%d phenotypes have no observed values and were left missing", n_cols_all_missing);
    }
    if ( n_cols_nonpositive_min > 0 ) {
        warning("%d phenotypes have a non-positive observed minimum, so half of the minimum is not below the observed values; half-minimum imputation assumes positive measurements", n_cols_nonpositive_min);
    }
    recompute_has_missing();
    return n_imputed;
}

bool PhenoMatrix::load_pheno_matrix(const char* pheno_file, const char* pheno_format, const char delim) {
    format.assign(pheno_format);
    if ( format.compare("regenie") == 0 ||
         format.compare("tsv-sample-row") == 0 ) { // sample in rows, phenos in columns
        notice("Loading sample-row format phenotype matrix from %s", pheno_file);
        int32_t icol_indid = -1;
        int32_t icol_data = -1;
        if ( format.compare("regenie") == 0 ) {
            icol_indid = 1;
            icol_data = 2;
        }
        else {
            icol_indid = 0;
            icol_data = 1;
        }
        tsv_reader tr(pheno_file);
        if ( delim != '\0' ) {
            tr.delimiter = delim; // '\0' would split on any whitespace
        }
        // read header
        int32_t nlines = 0;
        int32_t nfields = -1;
        int32_t nrow_est = 100;
        int32_t ntraits = -1;
        int32_t nsamps = 0;
        while ( tr.read_line() ) {
            if ( nlines == 0 ) {
                nfields = tr.nfields;
                check_header_fields(nfields, icol_data, delim, pheno_file);
                for ( int32_t i = icol_data; i < nfields; ++i ) {
                    pheno_ids.push_back( std::string(tr.str_field_at(i)) );
                }
                ntraits = nfields - icol_data;
                notice("Loading %d traits from phenotype matrix", ntraits);
                pheno_mat.resize(nrow_est, ntraits);
                pheno_mask.resize(nrow_est, ntraits);
                if ( !cat_ids.empty() || !selected_ids.empty() ) {
                    cat_strs.assign(ntraits, std::vector<std::string>());
                    for(int32_t j = 0; j < ntraits; ++j) {
                        if ( is_text_column(pheno_ids[j]) ) cat_strs[j].reserve(nrow_est);
                    }
                }
                ++nlines;
            }
            else {
                if ( nfields != tr.nfields ) {
                    error("Inconsistent number of fields at line %d in phenotype file %s", nlines + 1, pheno_file);
                }
                // check if expansion is needed
                if ( nsamps >= nrow_est ) {
                    nrow_est = nrow_est * 2;
                    notice("Expanding phenotype matrix to %d samples", nrow_est);
                    pheno_mat.conservativeResize(nrow_est, ntraits);
                    pheno_mask.conservativeResize(nrow_est, ntraits);
                }
                samp_ids.push_back( std::string(tr.str_field_at(icol_indid)) );
                for ( int32_t i = icol_data; i < nfields; ++i ) {
                    const bool is_cat = !cat_strs.empty() && is_text_column(pheno_ids[i - icol_data]);
                    if ( is_missing( tr.str_field_at(i) ) ) {
                        pheno_mask(nsamps, i - icol_data) = false;
                        pheno_mat(nsamps, i - icol_data) = 0;
                        has_missing = true;
                        if ( is_cat ) cat_strs[i - icol_data].push_back("");
                    }
                    else if ( is_cat ) { // a categorical level: any string
                        pheno_mask(nsamps, i - icol_data) = true;
                        pheno_mat(nsamps, i - icol_data) = 0;
                        cat_strs[i - icol_data].push_back(tr.str_field_at(i));
                    }
                    else {
                        pheno_mask(nsamps, i - icol_data) = true;
                        pheno_mat(nsamps, i - icol_data) = parse_pheno_value(tr.str_field_at(i), pheno_file, nlines + 1, i + 1);
                    }
                }
                ++nlines;
                ++nsamps;
            }
        }
        notice("Resizing phenotype matrix for %d samples and %d traits", nsamps, ntraits);
        // resize to the actual number of samples
        pheno_mat.conservativeResize(nsamps, ntraits);
        pheno_mask.conservativeResize(nsamps, ntraits);

        rebuild_id2index_map(samp_ids, samp_id2idx);
        rebuild_id2index_map(pheno_ids, pheno_id2idx);
        for(const std::string& id : cat_ids) {
            if ( pheno_id2idx.find(id) == pheno_id2idx.end() ) {
                error("Categorical covariate %s (--cat-cov-names) is not in %s", id.c_str(), pheno_file);
            }
        }
        return true;
    }
    else if ( format.compare("tsv-sample-col") == 0 ||
              format.compare("tensorqtl") == 0 ) {
        // sample in columns, phenos in rows
        int32_t icol_pheid = -1;
        int32_t icol_chrom = -1;
        int32_t icol_beg = -1;
        int32_t icol_end = -1;
        int32_t icol_data = -1;
        if ( format.compare("tensorqtl") == 0 ) {
            icol_chrom = 0;
            icol_beg = 1;
            icol_end = 2;
            has_loci = true;
            icol_pheid = 3;
            icol_data = 4;
        }
        else {
            icol_pheid = 0;
            icol_data = 1;
        }
        tsv_reader tr(pheno_file);
        if ( delim != '\0' ) {
            tr.delimiter = delim; // '\0' would split on any whitespace
        }
        // read header
        int32_t nlines = 0;
        int32_t nfields = -1;
        int32_t ncol_est = 100;
        int32_t nsamps = -1;
        int32_t ntraits = 0;
        while ( tr.read_line() ) {
            if ( nlines == 0 ) {
                nfields = tr.nfields;
                check_header_fields(nfields, icol_data, delim, pheno_file);
                if ( has_loci ) {
                    // header names are matched case-insensitively, with an optional leading '#'
                    // (tensorQTL writes "#chr", other BED producers write "#Chr" or "#CHROM")
                    const char* h_chr = tr.str_field_at(0);
                    if ( h_chr[0] == '#' ) ++h_chr;
                    if ( ( strcasecmp("chr", h_chr) != 0 && strcasecmp("chrom", h_chr) != 0 ) ||
                         strcasecmp("start", tr.str_field_at(1)) != 0 ||
                         strcasecmp("end", tr.str_field_at(2)) != 0 ) {
                        error("Expected chr,start,end columns for TensorQTL format in pheno file %s", pheno_file);
                    }
                }
                for ( int32_t i = icol_data; i < nfields; ++i ) {
                    samp_ids.push_back( std::string(tr.str_field_at(i)) );
                }
                nsamps = nfields - icol_data;
                pheno_mat.resize(nsamps, ncol_est);
                pheno_mask.resize(nsamps, ncol_est);
                nlines++;
            }
            else {
                if ( nfields != tr.nfields ) {
                    error("Inconsistent number of fields at line %d in phenotype file %s", nlines + 1, pheno_file);
                }
                // check if expansion is needed
                if ( ntraits >= ncol_est ) {
                    ncol_est = ncol_est * 2;
                    pheno_mat.conservativeResize(nsamps, ncol_est);
                    pheno_mask.conservativeResize(nsamps, ncol_est);
                }
                pheno_ids.push_back( std::string(tr.str_field_at(icol_pheid)) );
                if ( has_loci ) {
                    const char* chrom = tr.str_field_at(icol_chrom);
                    int32_t beg1 = tr.int_field_at(icol_beg);
                    int32_t end0 = tr.int_field_at(icol_end);
                    pheno_loci.push_back( genomeLocus(chrom, beg1, end0) );
                    pheno_locusmap.add(chrom, beg1, end0, ntraits);
                }
                const bool is_cat = is_text_column(pheno_ids.back());
                if ( !cat_ids.empty() || !selected_ids.empty() ) {
                    cat_strs.push_back(std::vector<std::string>());
                    if ( is_cat ) cat_strs.back().reserve(nsamps);
                }
                for ( int32_t i = icol_data; i < nfields; ++i ) {
                    if ( is_missing( tr.str_field_at(i) ) ) {
                        pheno_mask(i - icol_data, ntraits) = false;
                        pheno_mat(i - icol_data, ntraits) = 0;
                        has_missing = true;
                        if ( is_cat ) cat_strs.back().push_back("");
                    }
                    else if ( is_cat ) { // a categorical level: any string
                        pheno_mask(i - icol_data, ntraits) = true;
                        pheno_mat(i - icol_data, ntraits) = 0;
                        cat_strs.back().push_back(tr.str_field_at(i));
                    }
                    else {
                        pheno_mask(i - icol_data, ntraits) = true;
                        pheno_mat(i - icol_data, ntraits) = parse_pheno_value(tr.str_field_at(i), pheno_file, nlines + 1, i + 1);
                    }
                }
                ++nlines;
                ++ntraits;
            }
        }
        // resize to the actual number of samples
        pheno_mat.conservativeResize(nsamps, ntraits);
        pheno_mask.conservativeResize(nsamps, ntraits);

        rebuild_id2index_map(samp_ids, samp_id2idx);
        rebuild_id2index_map(pheno_ids, pheno_id2idx);
        for(const std::string& id : cat_ids) {
            if ( pheno_id2idx.find(id) == pheno_id2idx.end() ) {
                error("Categorical covariate %s (--cat-cov-names) is not in %s", id.c_str(), pheno_file);
            }
        }
        return true;
    }
    else {
        error("Unsupported phenotype format: %s", pheno_format);
        return false;
    }
}

void PhenoMatrix::keep_covariate_columns(const std::vector<std::string>& cov_names, const char* file) {
    if ( cov_names.empty() ) return;
    std::set<std::string> keep(cov_names.begin(), cov_names.end());
    for(const std::string& id : cov_names) {
        if ( pheno_id2idx.find(id) == pheno_id2idx.end() ) {
            error("Covariate %s (--cov-names) is not in %s", id.c_str(), file);
        }
    }
    keep.insert(cat_ids.begin(), cat_ids.end());
    std::vector<int32_t> idx;
    for(int32_t j = 0; j < (int32_t)pheno_ids.size(); ++j) {
        if ( keep.count(pheno_ids[j]) ) idx.push_back(j);
    }
    notice("Using %zu of %zu covariates in %s (--cov-names and --cat-cov-names)", idx.size(), pheno_ids.size(), file);
    subset_pheno_indices(idx);
}

int32_t PhenoMatrix::expand_categorical(int32_t max_levels) {
    if ( !has_categorical() ) { cat_strs.clear(); return 0; }
    const Eigen::Index n = pheno_mat.rows();
    std::vector<std::string> new_ids;
    std::vector<Eigen::VectorXd> new_cols;
    std::vector<Eigen::Matrix<bool, Eigen::Dynamic, 1> > new_masks;
    int32_t n_indicators = 0;
    for(size_t j = 0; j < pheno_ids.size(); ++j) {
        if ( cat_strs[j].empty() ) { // numeric column
            new_ids.push_back(pheno_ids[j]);
            new_cols.push_back(pheno_mat.col(j));
            new_masks.push_back(pheno_mask.col(j));
            continue;
        }
        std::set<std::string> levels;
        for(Eigen::Index i = 0; i < n; ++i) {
            if ( pheno_mask(i, j) ) levels.insert(cat_strs[j][i]);
        }
        if ( (int32_t)levels.size() > max_levels ) {
            error("Categorical covariate %s has %zu levels among the analyzed samples, more than --max-cat-levels %d",
                  pheno_ids[j].c_str(), levels.size(), max_levels);
        }
        std::string level_list;
        for(const std::string& l : levels) level_list += ( level_list.empty() ? "" : ", " ) + l;
        if ( levels.size() < 2 ) {
            warning("Categorical covariate %s has %zu level(s) among the analyzed samples (%s) and is dropped",
                    pheno_ids[j].c_str(), levels.size(), level_list.c_str());
            continue;
        }
        notice("Categorical covariate %s: %zu levels (%s); reference %s; %zu indicator columns",
               pheno_ids[j].c_str(), levels.size(), level_list.c_str(), levels.begin()->c_str(), levels.size() - 1);
        for(std::set<std::string>::const_iterator it = std::next(levels.begin()); it != levels.end(); ++it) {
            Eigen::VectorXd col(n);
            for(Eigen::Index i = 0; i < n; ++i) col(i) = ( pheno_mask(i, j) && cat_strs[j][i] == *it ) ? 1.0 : 0.0;
            new_ids.push_back(pheno_ids[j] + "_" + *it);
            new_cols.push_back(col);
            new_masks.push_back(pheno_mask.col(j));
            ++n_indicators;
        }
    }
    pheno_mat.resize(n, new_cols.size());
    pheno_mask.resize(n, new_cols.size());
    for(size_t j = 0; j < new_cols.size(); ++j) {
        pheno_mat.col(j) = new_cols[j];
        pheno_mask.col(j) = new_masks[j];
    }
    pheno_ids.swap(new_ids);
    has_loci = false;
    pheno_loci.clear();
    cat_strs.clear();
    rebuild_id2index_map(pheno_ids, pheno_id2idx);
    recompute_has_missing();
    return n_indicators;
}

bool PhenoMatrix::sample_ids_sorted() const {
    for ( size_t i = 1; i < samp_ids.size(); ++i ) {
        if ( samp_ids[i] < samp_ids[i-1] ) {
            return false;
        }
    }
    return true;
}

bool PhenoMatrix::pheno_ids_sorted() const {
    for ( size_t i = 1; i < pheno_ids.size(); ++i ) {
        if ( pheno_ids[i] < pheno_ids[i-1] ) {
            return false;
        }
    }
    return true;
}