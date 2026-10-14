#include "qpgen.h"
#include "qgenlib/qgen_utils.h"
#include <cstring>

// load PLINK files in a flexible manner.
// First, it searches for .pgen/.psam/.pvar files
// If not found, it searches for .bed/.bim/.fam files
bool PlinkReader::prep_prefix(const char*  _prefix) {
    // check if .pgen/.psam/.pvar files exist
    std::string str_prefix(_prefix);
    std::string _pgenf = str_prefix + ".pgen";
    std::string _psamf = str_prefix + ".psam";
    std::string _pvarf = str_prefix + ".pvar";

    if (check_file_existence(_pgenf.c_str()) && check_file_existence(_pvarf.c_str()) && check_file_existence(_psamf.c_str())) {
        return prep_pgen(_pgenf.c_str(), _pvarf.c_str(), _psamf.c_str());
    }
    // check if .bed/.bim/.fam files exist
    std::string _bedf = str_prefix + ".bed";
    std::string _bimf = str_prefix + ".bim";
    std::string _famf = str_prefix + ".fam";
    if (check_file_existence(_bedf.c_str()) && check_file_existence(_bimf.c_str()) && check_file_existence(_famf.c_str())) {
        return prep_bed(_bedf.c_str(), _bimf.c_str(), _famf.c_str());
    }
    // if none of the above, return false
    return false;
}

bool PlinkReader::prep_pgen(const char* _pgenf, const char* _pvarf, const char* _psamf) {
    // read the sample file first
    if ( !load_psam(_psamf) ) {
        error("Cannot parse the sample info at %s", _psamf);
    }

    mode_bed = false;   
    pgenf.assign(_pgenf);
    psamf.assign(_psamf);
    pvarf.assign(_pvarf);

    // load the variant file, but only partially
    //if ( !read_pvar(_pvarf) ) {
    //    error("Cannot parse the first variant at %s", _pvarf);
    //}

    // DO NOT LOAD the genotype file yet    
    return true;
}

bool PlinkReader::prep_bed(const char* _bedf, const char* _bimf, const char* _famf) {
    // read the sample file first
    if ( !load_fam(_famf) ) {
        error("Cannot parse the sample info at %s", _famf);
    }

    mode_bed = true;
    bedf.assign(_bedf);
    bimf.assign(_bimf);
    famf.assign(_famf);

    // load the variant file, but only partially
    //if ( !read_bim(_bimf) ) {
    //    error("Cannot parse the first variant at %s", _bimf);
    //}

    // DO NOT LOAD the genotype file yet    
    return true;
}

bool PlinkReader::load_psam(const char* _psamf) {
    tsv_reader tr_ind;

    if ( !tr_ind.open(_psamf) ) { // return false if file cannot be opened
        return false;
    }

    bool has_header = false;
    uint32_t nsamps = 0;
    int32_t idx_fid = -1, idx_iid = -1, idx_pat = -1, idx_mat = -1, idx_sex = -1, idx_pheno = -1;
    while( tr_ind.read_line() ) {
        if ( ( nsamps == 0 ) && ( !has_header ) ) { // first line to process
            if ( tr_ind.str_field_at(0)[0] == '#' ) {  // If header column exists
                for(int32_t i=0; i < tr_ind.nfields; ++i) {
                    const char* str = tr_ind.str_field_at(i) + ( i == 0 ? 1 : 0 );
                    if ( strcmp("FID", str) == 0 ) { idx_fid = i; }
                    else if ( strcmp("IID", str) == 0 ) { idx_iid = i; }
                    else if ( strcmp("PAT", str) == 0 ) { idx_pat = i; }
                    else if ( strcmp("MAT", str) == 0 ) { idx_mat = i; }
                    else if ( strcmp("SEX", str) == 0 ) { idx_sex = i; }
                    else if ( strcmp("PHENO", str) == 0 ) { idx_pheno = i; }
                }
                has_header = true;

                //notice("%d %d %d %d %d %d", idx_fid, idx_iid, idx_pat, idx_mat, idx_sex, idx_pheno);
                continue;
            }
            else {  // if header column does not exist, use the default order
                idx_fid = tr_ind.nfields > 0 ? 0 : -1;
                idx_iid = tr_ind.nfields > 1 ? 1 : -1;
                idx_pat = tr_ind.nfields > 2 ? 2 : -1;
                idx_mat = tr_ind.nfields > 3 ? 3 : -1;
                idx_sex = tr_ind.nfields > 4 ? 4 : -1;
                idx_pheno = tr_ind.nfields > 5 ? 5 : -1;
            }
        }

        // create a sample 
        plink_samp_t samp;
        if ( idx_fid >= 0 ) { samp.famID.assign(tr_ind.str_field_at(idx_fid)); } 
        if ( idx_iid >= 0 ) { samp.indID.assign(tr_ind.str_field_at(idx_iid)); }
        if ( idx_pat >= 0 ) { samp.dadID.assign(tr_ind.str_field_at(idx_pat)); }
        if ( idx_mat >= 0 ) { samp.momID.assign(tr_ind.str_field_at(idx_mat)); }
        if ( idx_sex >= 0 ) { samp.sex = tr_ind.int_field_at(idx_sex); }
        if ( idx_pheno >= 0 ) { samp.pheno = tr_ind.double_field_at(idx_pheno); } 

        // add the sample to the list
        samps.push_back(samp);
        // samples are matched by IID (as in the phenotype files), whatever the FID
        ++nsamps;
        if ( !samp2idx.insert(std::make_pair(samp.indID, (int32_t)nsamps)).second ) {
            error("Duplicate IID %s in %s. Samples are matched by IID, so IIDs must be unique", samp.indID.c_str(), _psamf);
        }
        //error("%s %s %s %d %d", iid.c_str(), samp.famID.c_str(), samp.indID.c_str(), idx_fid, idx_sex);
        //++nsamps;
    }

    // fill in the sample indices to load
    for(uint32_t i=0; i < nsamps; ++i) {
        samp_idx.push_back((int32_t)(i+1));
    }

    return nsamps > 0;
}

bool PlinkReader::load_fam(const char* _famf) {
    tsv_reader tr_ind;

    if ( !tr_ind.open(_famf) ) { // return false if file cannot be opened
        return false;
    }

    uint32_t nsamps = 0;
    while( tr_ind.read_line() ) {
        // create a sample 
        plink_samp_t samp;

        samp.famID.assign(tr_ind.str_field_at(0));
        samp.indID.assign(tr_ind.str_field_at(1));
        if ( tr_ind.nfields > 2 ) { samp.dadID.assign(tr_ind.str_field_at(2)); }
        if ( tr_ind.nfields > 3 ) { samp.momID.assign(tr_ind.str_field_at(3)); }
        if ( tr_ind.nfields > 4 ) { samp.sex = tr_ind.int_field_at(4); }
        if ( tr_ind.nfields > 5 ) { samp.pheno = tr_ind.double_field_at(5); }

        // add the sample to the list
        samps.push_back(samp);
        // samples are matched by IID (as in the phenotype files), whatever the FID
        ++nsamps;
        if ( !samp2idx.insert(std::make_pair(samp.indID, (int32_t)nsamps)).second ) {
            error("Duplicate IID %s in %s. Samples are matched by IID, so IIDs must be unique", samp.indID.c_str(), _famf);
        }
        //++nsamps;
    }

    // fill in the sample indices to load
    for(uint32_t i=0; i < nsamps; ++i) {
        samp_idx.push_back((int32_t)(i+1));
    }
    return true;
}

void PlinkReader::set_filter_sample_id(std::vector<std::string>& samp_ids, bool exclude) {
    samp_idx.clear();

    std::set<uint32_t> idxset;
    for(int32_t i=0; i < samp_ids.size(); ++i) {
        std::map<std::string, int32_t>::iterator it = samp2idx.find(samp_ids[i]);
        if ( it != samp2idx.end() ) {
            idxset.insert(it->second);
        }
    }

    for(uint32_t i=0; i < samps.size(); ++i) {
        if ( exclude ) {
            if ( idxset.find(i+1) == idxset.end() ) {
                samp_idx.push_back(i+1);
            }
        }
        else {
            if ( idxset.find(i+1) != idxset.end() ) {
                samp_idx.push_back(i+1);
            }
        }
    }

    if ( samp_idx.empty() ) {
        error("No samples to be included after subsetting to %zu", samp_ids.size());
    }
}

bool PlinkReader::stream_pvar() {
    if ( nvars == 0 ) {
        if ( !tr_var.open(pvarf.c_str()) ) { // return false if file cannot be opened
            error("Cannot open %s", pvarf.c_str());
            return false;
        } 
    }

    // ignore headers
    int32_t idx_chrom = -1, idx_pos = -1, idx_id = -1, idx_ref = -1, idx_alt = -1, idx_qual = -1, idx_cm = -1, idx_filter = -1, idx_format = -1, idx_info = -1;
    while( tr_var.read_line() ) {
        const char* s = tr_var.str_field_at(0);
        if ( s[0] != '#' ) {
            break;
        }
        else if ( strncmp(s, "#CHROM", 6) != 0 ) {  // meta line
            continue;
        }

        // this is a header line that contains column names
        for(int32_t i=0; i < tr_var.nfields; ++i) {
            const char* str = tr_var.str_field_at(i) + ( i == 0 ? 1 : 0 );
            if ( strcmp("CHROM", str) == 0 ) { idx_chrom = i; }
            else if ( strcmp("POS", str) == 0 ) { idx_pos = i; }
            else if ( strcmp("ID", str) == 0 ) { idx_id = i; }
            else if ( strcmp("REF", str) == 0 ) { idx_ref = i; }
            else if ( strcmp("ALT", str) == 0 ) { idx_alt = i; }
            else if ( strcmp("QUAL", str) == 0 ) { idx_qual = i; }
            else if ( strcmp("FILTER", str) == 0 ) { idx_filter = i; }
            else if ( strcmp("INFO", str) == 0 ) { idx_info = i; }
            else if ( strcmp("FORMAT", str) == 0 ) { idx_format = i; }
            else if ( strcmp("CM", str) == 0 ) { idx_cm = i; }
        }
    }

    if ( tr_var.fields == NULL ) return false;

    // read variant info
    plink_var_t var;
    var.schrom.assign(tr_var.str_field_at(idx_chrom));
    var.pos = tr_var.int_field_at(idx_pos);
    var.vid.assign(tr_var.str_field_at(idx_id));
    var.ref.assign(tr_var.str_field_at(idx_ref));
    split(var.alts, ",", tr_var.str_field_at(idx_alt));
    var.cM = tr_var.double_field_at(idx_cm);
    vars.push_back(var);
    ++nvars;

    return true;
}

bool PlinkReader::stream_bim() {
    if ( nvars == 0 ) {
        if ( !tr_var.open(bimf.c_str()) ) { // return false if file cannot be opened
            error("Cannot open %s", bimf.c_str());
            return false;
        } 
    }

    if ( tr_var.read_line() ) {
        plink_var_t var;
        var.schrom.assign(tr_var.str_field_at(0));
        var.vid.assign(tr_var.str_field_at(1));
        var.cM = tr_var.double_field_at(2);
        var.pos = tr_var.int_field_at(3);
        var.alts.push_back(tr_var.str_field_at(4));
        var.ref.assign(tr_var.str_field_at(5));
        vars.push_back(var);
        ++nvars;
        return true;
    }
    else {
        return false;
    }
}

bool PlinkReader::stream_bed_genos() {
    // read the variant position first
    if ( !stream_bim() ) {
        return false;
    }

    // load the genotypes if needed
    if ( !pgr_loaded ) {
        pgr.Load(bedf, (int32_t)samps.size(), samp_idx, nthreads);
        pgr_loaded = true;
    }

    // read the genotypes, read as integers
    pgr.ReadIntHardcalls(int_buf, 0, nvars-1, 0);
    return true;
}

bool PlinkReader::stream_pgen_genos() {
    // read the variant position first
    if ( !stream_pvar() ) {
        return false;
    }

    // load the genotypes if needed
    if ( !pgr_loaded ) {
        pgr.Load(pgenf, (int32_t)samps.size(), samp_idx, nthreads);
        pgr_loaded = true;
    }

    // read the genotypes, read as integers
    pgr.ReadIntHardcalls(int_buf, 0, nvars-1, 0);
    return true;
}

int32_t PlinkReader::get_variant_idx_from_cpra(const std::vector<std::string>& var_cpra, std::vector<int32_t>& variant_idx) {
    notice("bimf = %s", bimf.c_str());
    notice("pvarf = %s", pvarf.c_str());
    notice("mode_bed = %d", mode_bed);


    variant_idx.resize(var_cpra.size());
    std::fill(variant_idx.begin(), variant_idx.end(), -1);

    // create a map from CPRA to variant index
    std::map<std::string, int32_t> cpra2idx;
    for(int32_t i=0; i < var_cpra.size(); ++i) {
        cpra2idx[var_cpra[i]] = i;
    }

    // scan all variants from the 
    if ( mode_bed ) {
        return bim_cpra2idx(cpra2idx, variant_idx);
    } else {
        return pvar_cpra2idx(cpra2idx, variant_idx);
    }
}

int32_t PlinkReader::bim_cpra2idx(std::map<std::string, int32_t>& cpra2idx, std::vector<int32_t>& variant_idx) {
    tsv_reader tr_bim(bimf.c_str());
    char buf[65535];
    int32_t nv = 0, nmatch = 0;
    while ( tr_bim.read_line() ) {
        snprintf(buf, 65535, "%s:%s:%s:%s", tr_bim.str_field_at(0), tr_bim.str_field_at(3), tr_bim.str_field_at(5), tr_bim.str_field_at(4));
        if ( tr_bim.int_field_at(3) == 271129 ) 
            notice("buf = %s, cpra2idx.size() = %zu", buf, cpra2idx.size());
        std::map<std::string, int32_t>::iterator it = cpra2idx.find(buf);
        if ( it != cpra2idx.end() ) {
            variant_idx[it->second] = nv;
            ++nmatch;
        }
        ++nv;
    }
    tr_bim.close();
    return nmatch;
}

int32_t PlinkReader::pvar_cpra2idx(std::map<std::string, int32_t>& cpra2idx, std::vector<int32_t>& variant_idx) {
    tsv_reader tr_pvar(pvarf.c_str());
    char buf[65535];
    int32_t nv = 0, nmatch = 0;
    int32_t idx_chrom = -1, idx_pos = -1, idx_id = -1, idx_ref = -1, idx_alt = -1, idx_qual = -1, idx_cm = -1, idx_filter = -1, idx_format = -1, idx_info = -1;
    while ( tr_pvar.read_line() ) {
        const char* s = tr_pvar.str_field_at(0);
        if ( s[0] == '#' ) {
            if ( strncmp(s, "#CHROM", 6) == 0 ) {  // meta line
                // this is a header line that contains column names
                for(int32_t i=0; i < tr_pvar.nfields; ++i) {
                    const char* str = tr_pvar.str_field_at(i) + ( i == 0 ? 1 : 0 );
                    if ( strcmp("CHROM", str) == 0 ) { idx_chrom = i; }
                    else if ( strcmp("POS", str) == 0 ) { idx_pos = i; }
                    else if ( strcmp("ID", str) == 0 ) { idx_id = i; }
                    else if ( strcmp("REF", str) == 0 ) { idx_ref = i; }
                    else if ( strcmp("ALT", str) == 0 ) { idx_alt = i; }
                    else if ( strcmp("QUAL", str) == 0 ) { idx_qual = i; }
                    else if ( strcmp("FILTER", str) == 0 ) { idx_filter = i; }
                    else if ( strcmp("INFO", str) == 0 ) { idx_info = i; }
                    else if ( strcmp("FORMAT", str) == 0 ) { idx_format = i; }
                    else if ( strcmp("CM", str) == 0 ) { idx_cm = i; }
                }
                if ( ( idx_chrom < 0 ) || ( idx_pos < 0 ) || ( idx_ref < 0 ) || ( idx_alt < 0 ) ) {
                    error("Cannot find the required columns in %s", pvarf.c_str());
                }
            }
        }
        else {
            snprintf(buf, 65535, "%s:%s:%s:%s", tr_pvar.str_field_at(idx_chrom), tr_pvar.str_field_at(idx_pos), tr_pvar.str_field_at(idx_ref), tr_pvar.str_field_at(idx_alt));
            //error("%s %d %d %d %d", buf, idx_chrom, idx_pos, idx_ref, idx_alt);
            std::map<std::string, int32_t>::iterator it = cpra2idx.find(buf);
            if ( it != cpra2idx.end() ) {
                variant_idx[it->second] = nv;
                ++nmatch;
            }
            ++nv;
        }
    }
    tr_pvar.close();
    return nmatch;
}

bool PlinkReader::get_bed_genos_at(int32_t var_idx) {
    // load the genotypes if needed
    if ( !pgr_loaded ) {
        pgr.Load(bedf, (int32_t)samps.size(), samp_idx, nthreads);
        pgr_loaded = true;
    }

    // read the genotypes, read as integers
    pgr.ReadIntHardcalls(int_buf, 0, var_idx, 0);
    return true;
}

bool PlinkReader::get_pgen_genos_at(int32_t var_idx) {
    // load the genotypes if needed
    if ( !pgr_loaded ) {
        pgr.Load(pgenf, (int32_t)samps.size(), samp_idx, nthreads);
        pgr_loaded = true;
    }

    // read the genotypes, read as integers
    pgr.ReadIntHardcalls(int_buf, 0, var_idx, 0);
    return true;
}

bool PgenIdxReader::prep_pgen(const char* _pgenf, const char* _pivarf, const char* _psamf) {
    // read the sample file first
    if ( !load_psam(_psamf) ) {
        error("Cannot parse the sample info at %s", _psamf);
    }

    pgenf.assign(_pgenf);
    psamf.assign(_psamf);
    pivarf.assign(_pivarf);

    // DO NOT LOAD the genotype file yet    
    return true;
}

bool PgenIdxReader::load_psam(const char* _psamf) {
    tsv_reader tr_ind;

    if ( !tr_ind.open(_psamf) ) { // return false if file cannot be opened
        return false;
    }

    bool has_header = false;
    uint32_t nsamps = 0;
    int32_t idx_fid = -1, idx_iid = -1, idx_pat = -1, idx_mat = -1, idx_sex = -1, idx_pheno = -1;
    samps.clear();
    samp2idx.clear();
    while( tr_ind.read_line() ) {
        if ( ( nsamps == 0 ) && ( !has_header ) ) { // first line to process
            if ( tr_ind.str_field_at(0)[0] == '#' ) {  // If header column exists
                for(int32_t i=0; i < tr_ind.nfields; ++i) {
                    const char* str = tr_ind.str_field_at(i) + ( i == 0 ? 1 : 0 );
                    if ( strcmp("FID", str) == 0 ) { idx_fid = i; }
                    else if ( strcmp("IID", str) == 0 ) { idx_iid = i; }
                    else if ( strcmp("PAT", str) == 0 ) { idx_pat = i; }
                    else if ( strcmp("MAT", str) == 0 ) { idx_mat = i; }
                    else if ( strcmp("SEX", str) == 0 ) { idx_sex = i; }
                    else if ( strcmp("PHENO", str) == 0 ) { idx_pheno = i; }
                }
                has_header = true;

                //notice("%d %d %d %d %d %d", idx_fid, idx_iid, idx_pat, idx_mat, idx_sex, idx_pheno);
                continue;
            }
            else {  // if header column does not exist, use the default order
                idx_fid = tr_ind.nfields > 0 ? 0 : -1;
                idx_iid = tr_ind.nfields > 1 ? 1 : -1;
                idx_pat = tr_ind.nfields > 2 ? 2 : -1;
                idx_mat = tr_ind.nfields > 3 ? 3 : -1;
                idx_sex = tr_ind.nfields > 4 ? 4 : -1;
                idx_pheno = tr_ind.nfields > 5 ? 5 : -1;
            }
        }

        // create a sample 
        plink_samp_t samp;
        if ( idx_fid >= 0 ) { samp.famID.assign(tr_ind.str_field_at(idx_fid)); } 
        if ( idx_iid >= 0 ) { samp.indID.assign(tr_ind.str_field_at(idx_iid)); }
        if ( idx_pat >= 0 ) { samp.dadID.assign(tr_ind.str_field_at(idx_pat)); }
        if ( idx_mat >= 0 ) { samp.momID.assign(tr_ind.str_field_at(idx_mat)); }
        if ( idx_sex >= 0 ) { samp.sex = tr_ind.int_field_at(idx_sex); }
        if ( idx_pheno >= 0 ) { samp.pheno = tr_ind.double_field_at(idx_pheno); } 

        // add the sample to the list
        samps.push_back(samp);
        // samples are matched by IID (as in the phenotype files), whatever the FID
        ++nsamps;
        if ( !samp2idx.insert(std::make_pair(samp.indID, (int32_t)nsamps)).second ) {
            error("Duplicate IID %s in %s. Samples are matched by IID, so IIDs must be unique", samp.indID.c_str(), _psamf);
        }
        //error("%s %s %s %d %d", iid.c_str(), samp.famID.c_str(), samp.indID.c_str(), idx_fid, idx_sex);
        //++nsamps;
    }

    // fill in the sample indices to load
    for(uint32_t i=0; i < nsamps; ++i) {
        samp_idx.push_back((int32_t)(i+1));
    }

    return nsamps > 0;
}

void PgenIdxReader::subset_sample_indices(const std::vector<int32_t>& samp_indices, bool exclude) {
    if ( exclude ) {
        samp_idx.clear();
        std::set<uint32_t> idxset;
        for(int32_t i=0; i < samp_indices.size(); ++i) {
            idxset.insert(samp_indices[i]-1); // convert to 0-based index
        }
        for(int32_t i=0; i < samps.size(); ++i) {
            if ( idxset.find(i) == idxset.end() ) { // if not in the exclude set
                samp_idx.push_back(i+1); // convert to 1-based index
            }
        }
    }
    else {
        samp_idx = samp_indices; // copy the indices
    }
    // an empty subset would make the reader load every sample, misaligning the genotypes
    // with the phenotypes
    if ( samp_idx.empty() ) {
        error("No samples are left after subsetting the samples of %s", psamf.c_str());
    }
}


void PgenIdxReader::subset_sample_ids(const std::vector<std::string>& samp_ids, bool exclude) {
    std::vector<int32_t> idx;
    int32_t n_miss = 0;
    for(int32_t i=0; i < samp_ids.size(); ++i) {
        std::map<std::string, int32_t>::iterator it = samp2idx.find(samp_ids[i]);
        if ( it != samp2idx.end() ) {
            idx.push_back(it->second); // convert to 1-based index
        }
        else {
            ++n_miss;
        }
    }        
    if ( n_miss > 0 ) {
        warning("%d of %zu requested samples are not found in the sample file %s", n_miss, samp_ids.size(), psamf.c_str());
    }
    subset_sample_indices(idx, exclude);
}

bool PgenIdxReader::read_pos(const char* chrom, int32_t pos) {
    // load the variant file if not loaded
    if ( !pivar_loaded ) {
        if ( !tr_pivar.open(pivarf.c_str()) ) { // return false if file cannot be opened
            error("Cannot open %s", pivarf.c_str());
            return false;
        } 
        pivar_loaded = true;
    }

    notice("Jumping to %s:%d in %s, icol_pivar_idx: %d", chrom, pos, pivarf.c_str(), icol_pivar_idx);

    tr_pivar.jump_to(chrom, pos);
    while( tr_pivar.read_line() ) { // find the variant
        if ( tr_pivar.nfields < icol_pivar_idx ) {
            error("Invalid pvar file format at %s", pivarf.c_str());
            return false;
        }
        const char* chrom2 = tr_pivar.str_field_at(0);
        int32_t pos2 = tr_pivar.int_field_at(1);
        last_pivar_chrom.assign(chrom2);
        last_pivar_pos = pos2;
        if ( strcmp(chrom2, chrom) != 0 ) { // something is wrong. No variant found 
            return false;           
        }
        if ( pos2 < pos ) {
            continue;
        }
        // found the variant
        cur_var.schrom.assign(chrom2);
        cur_var.pos = pos2;
        cur_var.vid.assign(tr_pivar.str_field_at(2));
        cur_var.ref.assign(tr_pivar.str_field_at(3));
        split(cur_var.alts, ",", tr_pivar.str_field_at(4));
        cur_var_idx = tr_pivar.int_field_at(icol_pivar_idx)-1; 
        return true;
    }
    return false;
}

bool PgenIdxReader::read_pivar(const char* cpra) {
    // load the variant file if not loaded
    if ( !pivar_loaded ) {
        if ( !tr_pivar.open(pivarf.c_str()) ) { // return false if file cannot be opened
            error("Cannot open %s", pivarf.c_str());
            return false;
        } 
        pivar_loaded = true;
    }

    if ( cpra != NULL ) {
        cpra_t cpra_obj(cpra);
        // Decide between streaming forward and a tabix jump from the last pvar line actually
        // consumed, not from cur_var (which only changes on a successful lookup): a failed
        // lookup leaves the cursor past the requested position, and a repeated request or a
        // request at the same position must re-read lines that were already consumed.
        if ( last_pivar_pos < 0 ||
             last_pivar_chrom.compare(cpra_obj.chrom) != 0 ||
             cpra_obj.pos <= last_pivar_pos ||
             cpra_obj.pos - last_pivar_pos > jump_thres_bp ) {
            tr_pivar.jump_to(cpra_obj.chrom.c_str(), cpra_obj.pos);
        }
        while( tr_pivar.read_line() ) { // find the variant
            if ( tr_pivar.nfields <= icol_pivar_idx ) {
                error("Invalid pvar file format at %s", pivarf.c_str());
                return false;
            }
            const char* chrom = tr_pivar.str_field_at(0);
            int32_t pos = tr_pivar.int_field_at(1);
            last_pivar_chrom.assign(chrom);
            last_pivar_pos = pos;
            const char* ref = tr_pivar.str_field_at(3);
            const char* alts = tr_pivar.str_field_at(4);
            if ( cpra_obj.chrom.compare(chrom) == 0 && 
                 pos == cpra_obj.pos && 
                 cpra_obj.ref.compare(ref) == 0 &&
                 cpra_obj.alts.compare(alts) == 0 ) {
                // found the variant
                cur_var.schrom.assign(chrom);
                cur_var.pos = pos;
                cur_var.vid.assign(tr_pivar.str_field_at(2));
                cur_var.ref.assign(ref);
                split(cur_var.alts, ",", alts);
                cur_var_idx = tr_pivar.int_field_at(icol_pivar_idx)-1; 
                return true;
            }
            else if ( cpra_obj.chrom.compare(chrom) != 0 || pos > cpra_obj.pos ) {
                return false; // variant not found
            }
        }
    }
    else { // read the next variant
        while( tr_pivar.read_line() ) { // find the variant
            const char* chrom = tr_pivar.str_field_at(0);
            if ( chrom[0] == '#' ) continue;

            if ( tr_pivar.nfields <= icol_pivar_idx ) {
                error("Invalid pvar file format at %s", pivarf.c_str());
                return false;
            }

            int32_t pos = tr_pivar.int_field_at(1);
            last_pivar_chrom.assign(chrom);
            last_pivar_pos = pos;
            const char* ref = tr_pivar.str_field_at(3);
            const char* alts = tr_pivar.str_field_at(4);

            cur_var.schrom.assign(chrom);
            cur_var.pos = pos;
            cur_var.vid.assign(tr_pivar.str_field_at(2));
            cur_var.ref.assign(ref);
            split(cur_var.alts, ",", alts);
            cur_var_idx = tr_pivar.int_field_at(icol_pivar_idx)-1; 
            return true;
        }
    }
    return false; // reached end of file without finding the next variant
}

// Reads the current (or specified) variant, opportunistically keeping the
// sparse (difflist) representation for hardcall-only variants. When the variant
// is read sparsely, int_buf is NOT expanded; consumers should use is_sparse() /
// get_sparse_*() instead. When read densely, int_buf (or dbl_buf for dosages) is
// fully populated as usual.
bool PgenIdxReader::get_genos_sparse(int32_t var_idx) {
    // load the genotypes if needed
    if ( !pgen_loaded ) {
        notice("Loading pgen file %s with %zu/%zu samples", pgenf.c_str(), samp_idx.size(), samps.size());
        pgr.Load(pgenf, (int32_t)samps.size(), samp_idx, nthreads);
        pgen_loaded = true;
        dosage_present = pgr.DosagePresent();
        dbl_buf = (double*)malloc(sizeof(double) * samp_idx.size());
    }

    int32_t vidx = var_idx < 0 ? cur_var_idx : var_idx;

    if ( dosage_present ) {
        // dosage variants have no sparse hardcall representation; read densely
        cur_is_sparse = false;
        if ( mean_impute_dosage ) {
            // missing dosages are replaced by the variant mean dosage
            pgr.ReadMeanimpute(dbl_buf, (size_t)samp_idx.size(), 0, vidx, 0);
        }
        else {
            pgr.Read(dbl_buf, (size_t)samp_idx.size(), 0, vidx, 0);
        }
    }
    else {
        // returns true (sparse: sparse_* filled) or false (dense: int_buf filled)
        cur_is_sparse = pgr.ReadMaybeSparseHardcalls(0, vidx, &sparse_common_geno, sparse_sample_idxs, sparse_genos, int_buf);
    }
    return true;
}

bool PgenIdxReader::get_genos(int32_t var_idx) {
    bool ret = get_genos_sparse(var_idx);
    // backward-compatible behavior: always provide the fully expanded int_buf
    if ( ret && cur_is_sparse ) {
        int32_t n = (int32_t)samp_idx.size();
        if ( (int32_t)int_buf.size() != n ) int_buf.resize(n);
        std::fill(int_buf.begin(), int_buf.end(), sparse_common_geno);
        for(size_t k=0; k < sparse_sample_idxs.size(); ++k) {
            int_buf[sparse_sample_idxs[k]] = sparse_genos[k];
        }
        cur_is_sparse = false; // int_buf now holds the dense representation
    }
    return ret;
}

bool PgenIdxReader::sample_ids_sorted() const {
    // sort based on indID only, since famID can be empty or "0"
    for(int32_t i=1; i < samps.size(); ++i) {
        if ( samps[i-1].indID.compare(samps[i].indID) > 0 ) {
            return false; 
        }
    }
    return true;
}

bool MultiPgenIdxReader::prep_pgen_list(const char* listf, const char* pgen_suffix, const char* pivar_suffix, const char* psam_suffix) {
    // read the pgen list file
    tsv_reader tr(listf);
    // add one pgen file at a time
    // while ensuring that the sample ids are consistent
    // return false if any error occurs
    while ( tr.read_line() ) {
        const char* chrom = tr.str_field_at(0);
        int32_t beg = tr.int_field_at(1);
        int32_t end = tr.int_field_at(2);
        if ( tr.nfields == 4 ) { // assume that the PLINK prefix is given as input
            std::string prefix = tr.str_field_at(3);
            std::string pgenf = prefix + pgen_suffix;
            std::string pvarf = prefix + pivar_suffix;
            std::string psamf = prefix + psam_suffix;
            add_pgen(chrom, beg, end, pgenf.c_str(), pvarf.c_str(), psamf.c_str());
        }
        else if ( tr.nfields == 6 ) { // CHROM BEG END PGEN PSAM PIVAR, as documented in --help and docs/
            std::string pgenf = tr.str_field_at(3);
            std::string psamf = tr.str_field_at(4);
            std::string pvarf = tr.str_field_at(5);
            add_pgen(chrom, beg, end, pgenf.c_str(), pvarf.c_str(), psamf.c_str());
        }
        else {
            error("Invalid pgen list file format at %s. Expecting 4 or 6 fields, got %zu", listf, tr.nfields);
            return false;
        }
    }

    // check if loci are overlapping
    genomeLocus prev_locus("", 0, 0);
    bool is_beginning = true;
    for(locus2idx.rewind(); !locus2idx.isend(); locus2idx.next()) {
        if ( !is_beginning ) {
            // check if the neighboring locus are overlapping
            if ( locus2idx.it->first.overlaps(prev_locus) ) {
                error("Overlapping loci found: %s and %s. Please use PLINK files with non-overlapping regions", prev_locus.toString(), locus2idx.it->first.toString());
                return false;
            }
        }
        prev_locus = locus2idx.it->first;
        is_beginning = false;
    }
    return true;
}

bool MultiPgenIdxReader::add_pgen(const char* chrom, int32_t beg, int32_t end, const char* pgenf, const char* pivarf, const char* psamf) {
    // add one pgen file at a time
    // make sure that the regions are non-overlapping
    // make sure that the sample ids are consistent
    // return false if any error occurs
    PgenIdxReader* p_reader = new PgenIdxReader();
    if ( ! p_reader->prep_pgen(pgenf, pivarf, psamf) ) {
        error("Failed to prepare pgen/pivar/psam files %s/%s/%s", pgenf, pivarf, psamf);
    }
    int32_t idx = p_readers.size();
    if ( idx > 0 ) {
        // make sure that the sample sizes are consistent
        if ( p_reader->get_loaded_sample_count() != p_readers[0]->get_loaded_sample_count() ) {
            error("Sample sizes do not match between pgen files");
            return false;
        }

        // make sure that the sample ids are consistent
        int32_t n = p_reader->get_loaded_sample_count();
        const std::vector<plink_samp_t>& samp_ids = p_reader->get_all_samples();
        const std::vector<plink_samp_t>& samp_ids_0 = p_readers[0]->get_all_samples();
        for(int32_t i=0; i < n; ++i) {
            if ( samp_ids[i].indID != samp_ids_0[i].indID ) {
                error("Sample ids do not match in PSAM file %s", psamf);
                return false;
            }
        }
    }
    p_readers.push_back(p_reader);
    loci.push_back(genomeLocus(chrom, beg, end));
    locus2idx.add(chrom, beg, end, idx);

    // how do I effciently ensure that the loci are non-overlapping?
    return true;
}

bool MultiPgenIdxReader::set_single_chunk_pgen(const char* pgenf, const char* pivarf, const char* psamf) {
    if ( p_readers.size() > 0 ) {
        error("Single chunk pgen files cannot be set with multiple pgen files");
        return false;
    }
    PgenIdxReader* p_reader = new PgenIdxReader();
    if ( ! p_reader->prep_pgen(pgenf, pivarf, psamf) ) {
        error("Failed to prepare pgen/pivar/psam files %s/%s/%s", pgenf, pivarf, psamf);
    }
    p_readers.push_back(p_reader);

    single_chunk_mode = true;
    idx_cur_reader = 0;
    return true;
}

void MultiPgenIdxReader::subset_sample_ids(const std::vector<std::string>& samp_ids, bool exclude) {
    for(int32_t i=0; i < p_readers.size(); ++i) {
        p_readers[i]->subset_sample_ids(samp_ids, exclude);
        if ( p_readers[i]->get_loaded_sample_count() != p_readers[0]->get_loaded_sample_count() ) {
            error("Sample sizes do not match between pgen files after subsetting");
        }
    }   
}

bool MultiPgenIdxReader::read_pos(const char* chrom, int32_t pos) { // change the current variant position to a specific CPRA
    // find the locus that contains the position
    if ( single_chunk_mode ) {
        return p_readers[0]->read_pos(chrom, pos);
    }
    else if ( locus2idx.moveTo(chrom, pos) ) {
        // overlapping region exists
        int32_t idx = locus2idx.it->second;
        idx_cur_reader = idx;
        return p_readers[idx]->read_pos(chrom, pos);
    }
    return false;
}

bool MultiPgenIdxReader::read_pivar(const char* cpra) {      // change the current variant position to a specific CPRA
    if ( cpra == NULL ) {
        if ( single_chunk_mode ) {
            return p_readers[0]->read_pivar();
        }
        else if ( idx_cur_reader >= 0 ) {
            bool ret = p_readers[idx_cur_reader]->read_pivar();
            while ( !ret ) {
                locus2idx.next();
                if ( locus2idx.isend() ) {
                    return false;
                }
                idx_cur_reader = locus2idx.it->second;
                ret = p_readers[idx_cur_reader]->read_pivar();
            }
            return true;
        }
        else {
            return false;
        }
    }
    else {
        if ( single_chunk_mode ) {
            return p_readers[0]->read_pivar(cpra);
        }
        else {
            cpra_t cpra_obj(cpra);
            if ( locus2idx.moveTo(cpra_obj.chrom.c_str(), cpra_obj.pos) ) {
                int32_t idx = locus2idx.it->second;
                idx_cur_reader = idx;
                return p_readers[idx]->read_pivar(cpra);
            }
            return false;
        }
    }
}

bool MultiPgenIdxReader::get_genos() {                               // read the genotypes at the current variant position
    if ( idx_cur_reader >= 0 ) {
        bool ret = p_readers[idx_cur_reader]->get_genos();
        if ( ret ) {
            if ( dosage_present ) {
                dbl_buf = (double*)p_readers[idx_cur_reader]->get_dbl_buf();
            }
            else {
                //std::copy(p_readers[idx_cur_reader]->get_int_buf().begin(), p_readers[idx_cur_reader]->get_int_buf().end(), int_buf.begin());
                int_buf = p_readers[idx_cur_reader]->get_int_buf();
            }
        }
        return ret;
    }
    return false;
}

bool MultiPgenIdxReader::get_genos_sparse() {                        // read genotypes, keeping the sparse representation when possible
    if ( idx_cur_reader >= 0 ) {
        PgenIdxReader* pr = p_readers[idx_cur_reader];
        bool ret = pr->get_genos_sparse();
        if ( ret ) {
            if ( dosage_present ) {
                dbl_buf = (double*)pr->get_dbl_buf();
            }
            else if ( !pr->is_sparse() ) {
                // dense hardcalls: mirror get_genos() copy behavior
                int_buf = pr->get_int_buf();
            }
            // sparse: exposed via is_sparse()/get_sparse_*() delegating accessors
        }
        return ret;
    }
    return false;
}

// bool MultiPgenIdxReader::compute_geno_stats() {
//     if ( idx_cur_reader >= 0 ) {
//         bool ret = p_readers[idx_cur_reader]->compute_geno_stats();
//         return ret;   
//     }
//     return false;
// }

void MultiPgenIdxReader::set_n_threads(int32_t n) { 
    nthreads = n; 
    for(int32_t i=0; i < p_readers.size(); ++i) {
        p_readers[i]->set_n_threads(n);
    }
}

void MultiPgenIdxReader::set_icol_pivar_idx(int32_t idx) {
    icol_pivar_idx = idx;
    for(int32_t i=0; i < p_readers.size(); ++i) {
        p_readers[i]->set_icol_pivar_idx(idx);
    }
}

const plink_var_t& MultiPgenIdxReader::get_current_variant() const {
    if ( idx_cur_reader < 0 ) {
        error("No PGEN file is currently active");
    }
    return p_readers[idx_cur_reader]->get_current_variant();
}

int32_t MultiPgenIdxReader::get_all_sample_count() const {
    if ( p_readers.size() == 0 ) {
        error("No PGEN files are added yet");
    }
    return p_readers[0]->get_all_sample_count();
}
const std::vector<plink_samp_t>& MultiPgenIdxReader::get_all_samples() {
    if ( p_readers.size() == 0 ) {
        error("No PGEN files are added yet");
    }
    return p_readers[0]->get_all_samples();
}

int32_t MultiPgenIdxReader::get_loaded_sample_count() const {
    if ( p_readers.size() == 0 ) {
        error("No PGEN files are added yet");
    }
    return p_readers[0]->get_loaded_sample_count();
}

bool MultiPgenIdxReader::sample_ids_sorted() const {
    if ( p_readers.size() == 0 ) {
        error("No PGEN files are added yet");
    }
    return p_readers[0]->sample_ids_sorted();
}

//const std::vector<int32_t>& MultiPgenIdxReader::get_loaded_sample_indices() const;
const plink_samp_t& MultiPgenIdxReader::get_loaded_sample(int32_t idx) const {
    if ( p_readers.size() == 0 ) {
        error("No PGEN files are added yet");
    }
    return p_readers[0]->get_loaded_sample(idx);
}
