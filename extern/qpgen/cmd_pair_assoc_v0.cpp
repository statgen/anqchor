#include "qgenlib/params.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_utils.h"
#include "qgenlib/phred_helper.h"
#include "qgenlib/hts_utils.h"
#include "assoc_utils.h"
#include "qpgen_utils.h"
#include "qpgen.h"
#include "Eigen/Dense"
#include <cmath>

int32_t cmd_pair_assoc_v0(int32_t argc, char **argv)
{
    std::string pgenf;      // PLINK2 genotype file
    std::string psamf;      // PLINK2 sample file
    std::string pivarf;     // PLINK2 indexed pvar file, compatible with PgenIdxReader
    std::string phef;
    std::string covf;
    std::string pairf;
    std::string outf;
    std::string samplef;
    int32_t jump_thres_bp = 1000000; 
    int32_t max_chunk_vars = 1000; // Maximum number of variants to store at once in memory
    int32_t offset_pheno = 1; 
    int32_t offset_cov = 1;
    int32_t icol_pivar_idx = 9;
    int32_t icol_pheno_id = 1;
    int32_t icol_cov_id = 1;
    std::string colname_pheno_sample("pheno"); // column name for the sample IDs in the phenotype file
    std::string colname_geno_sample("geno"); // column name for the sample IDs in the genotype file 
    std::string colname_pair_trait("trait"); // column name for the trait ID in the pair file
    std::string colname_pair_variant("variant"); // column name of the variant ID in the pair file
    std::string colname_pair_region("region"); // column name of the region string in the pair file
    bool skip_rint = false; // skip tests based on rank-based inverse normal transformation

    paramList pl;

    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Options", NULL)
    LONG_STRING_PARAM("pgen", &pgenf, "Input PLINK2 genotype file")
    LONG_STRING_PARAM("psam", &psamf, "Input PLINK2 sample file")
    LONG_STRING_PARAM("pivar", &pivarf, "Input PLINK2 index pvar file (bgzipped and tabix)")
    LONG_STRING_PARAM("pheno", &phef, "Input phenotype matrix")
    LONG_STRING_PARAM("pairs", &pairf, "Input file containing [variant_id] [pheno_ID] pairs to be tested")
    LONG_STRING_PARAM("sample", &samplef, "Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files")
    LONG_STRING_PARAM("cov", &covf, "Input covariate matrix (optional)")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output prefix")

    LONG_PARAM_GROUP("Auxiliary options", NULL)
    LONG_INT_PARAM("jump-thres-bp", &jump_thres_bp, "Jump threshold in base pairs for the variant index (default: 1000000)")
    LONG_INT_PARAM("max-chunk-vars", &max_chunk_vars, "Maximum number of variants to store at once in memory (default: 1000)")
    LONG_INT_PARAM("offset-pheno", &offset_pheno, "The number of index columns (i.e. not containing values) in the phenotype files (default: 1)")
    LONG_INT_PARAM("offset-cov", &offset_cov, "The number of index columns (i.e. not containing values) in the covariate columns (default: 1)")
    LONG_INT_PARAM("icol-pivar-idx", &icol_pivar_idx, "1-based column index for the variant ID in the pvar file (default: 9)")
    LONG_INT_PARAM("icol-pheno-id", &icol_pheno_id, "1-based column index for the phenotype ID in the phenotype file (default: 1)")
    LONG_INT_PARAM("icol-cov-id", &icol_cov_id, "1-based column index for the covariate ID in the phenotype file (default: 1)")
    LONG_STRING_PARAM("colname-pheno-sample", &colname_pheno_sample, "When --sample is provided, the column name for the sample IDs in the phenotype file (default: 'pheno')")
    LONG_STRING_PARAM("colname-geno-sample", &colname_geno_sample, "When --sample is provided, the column name for the sample IDs in the phenotype file (default: 'geno')")
    LONG_STRING_PARAM("colname-pair-trait", &colname_pair_trait, "Column name of the trait ID in the pait file")
    LONG_STRING_PARAM("colname-pair-variant", &colname_pair_variant, "Column name of the variant ID in the pair file")
    LONG_STRING_PARAM("colname-pair-region", &colname_pair_region, "Column name of the region string in the pair file")
    LONG_PARAM("skip-rint", &skip_rint, "Skip tests based on rank-based inverse normal transformation (default: false)")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    // open the Plink file
    PgenIdxReader pr;
    if (!pr.prep_pgen(pgenf.c_str(), pivarf.c_str(), psamf.c_str()))
        error("Failed to prepare pgen files with the following files:\n%s\n%s\n%s", pgenf.c_str(), pivarf.c_str(), psamf.c_str());
    pr.set_jump_thres_bp(jump_thres_bp);
    pr.set_icol_pivar_idx(icol_pivar_idx - 1); // convert to 0-based index

    // load the sample IDs first
    if (!pr.load_psam(psamf.c_str()))
        error("Failed to load the sample IDs from %s", psamf.c_str());

    const std::map<std::string, int32_t>& psam_samp2idx = pr.get_samp2idx();

    // Rules for determining the sample IDs to analyze:
    // 1. If --sample is probided, use exactly the sample IDs provided in the input file. If any sample is missing, report an error
    //    The order of samples should be the same as the input file
    // 2. If --sample is not provided, use the overlapping sample IDs among genotypes, phenotypes, and covariates

    std::vector<std::string> list_pheno_sample_ids;         // list of sample IDs in the phenotype file
    std::vector<std::string> list_geno_sample_ids;          // list of sample IDs in the genotype file
    std::map<std::string, int32_t> list_pheno_sample2idx;   
    std::map<std::string, int32_t> list_geno_sample2idx;
    std::map<std::string, std::string> list_pheno2geno; // map from phenotype sample ID to genotype sample ID
    bool sample_list_provided = !samplef.empty();  // indicates if the sample file is provided
    if ( sample_list_provided ) { 
        tsv_reader tr_sample(samplef.c_str());
        int32_t icol_pheno = -1;
        int32_t icol_geno = -1;
        while(tr_sample.read_line()) {
            if ( icol_pheno < 0 ) { 
                // read the header line
                if ( tr_sample.nfields < 2 )
                    error("The sample file %s should have at least two columns", samplef.c_str());
                for(int32_t i = 0; i < tr_sample.nfields; ++i) {
                    if ( tr_sample.str_field_at(i) == colname_pheno_sample ) {
                        icol_pheno = i;
                    }
                    else if ( tr_sample.str_field_at(i) == colname_geno_sample ) {
                        icol_geno = i;
                    }
                }
                if ( icol_pheno < 0 || icol_geno < 0 )
                    error("Cannot find %s or %s in the first line of %s, which is required", colname_pheno_sample.c_str(), colname_geno_sample.c_str(), samplef.c_str());
            }
            else {
                // check if the sample ID exists in both
                const char* pheno_id = tr_sample.str_field_at(icol_pheno);
                const char* geno_id = tr_sample.str_field_at(icol_geno);
                list_pheno_sample_ids.push_back(pheno_id);
                list_geno_sample_ids.push_back(geno_id);
                if ( list_pheno_sample2idx.find(pheno_id) != list_pheno_sample2idx.end() )
                    error("Duplicate sample ID %s found in the phenotype sample file %s", pheno_id, samplef.c_str());
                if ( list_geno_sample2idx.find(geno_id) != list_geno_sample2idx.end() )
                    error("Duplicate sample ID %s found in the genotype sample file %s", geno_id, samplef.c_str());
                list_pheno_sample2idx[pheno_id] = (int32_t)list_pheno_sample_ids.size() - 1; // store the index of the sample ID
                list_geno_sample2idx[geno_id] = (int32_t)list_geno_sample_ids.size() - 1; // store the index of the sample ID
                list_pheno2geno[pheno_id] = geno_id; // map from phenotype sample ID to genotype sample ID

                // make sure that geno_id exists in the genotype file
                if ( psam_samp2idx.find(geno_id) == psam_samp2idx.end() ) {
                    error("Sample ID %s in file %s does not exist in the genotype file %s", geno_id, samplef.c_str(), psamf.c_str());
                }
            }
        }
        tr_sample.close();       
    }

    // read the phenotype file
    tsv_reader tr_phe(phef.c_str());
    // read the header line
    if (!tr_phe.read_line())
        error("Cannot read the header line from the phenotype file %s", phef.c_str());
    std::vector<std::string> phe_sample_ids;    // all sample IDs from the phenotypen file
    std::map<std::string, int32_t> phe_id2idx;  // 0-based indices of each phenotype sample
    for(int32_t i = offset_pheno; i < tr_phe.nfields; ++i) {
        phe_sample_ids.push_back(tr_phe.str_field_at(i));
        phe_id2idx[tr_phe.str_field_at(i)] = i - offset_pheno; // store the index of the sample ID
    }

    // read the covariate file if exists
    std::vector<int32_t> cov_idxs;  // 0-based
    std::vector<int32_t> phe_idxs;  // 0-based
    std::vector<int32_t> geno_idxs; // 1-based 
    std::vector<std::string> cov_trait_ids;
    std::vector<std::string> cov_sample_ids;

    Eigen::MatrixXd cov_mat;
    if (!covf.empty())
    {
        tsv_reader tr_cov(covf.c_str());
        // read the header line
        if (!tr_cov.read_line())
            error("Cannot read the header line from the covariate file %s", covf.c_str());

        std::map<std::string, int32_t> cov_id2idx;
        for(int32_t i = offset_cov; i < tr_cov.nfields; ++i) {
            cov_sample_ids.push_back(tr_cov.str_field_at(i));
            cov_id2idx[tr_cov.str_field_at(i)] = i - offset_cov; // store the index of the sample ID
        }

        // now determine cov_idxs, phe_idxs, and geno_idx 
        if ( sample_list_provided ) { // when sample list is provided, make sure that all samples exist
            fill_list_indices(list_pheno_sample_ids, cov_id2idx, cov_idxs);
            fill_list_indices(list_pheno_sample_ids, phe_id2idx, phe_idxs);
            fill_list_indices(list_geno_sample_ids, psam_samp2idx, geno_idxs);

            // need to make geno_idxs strictly increasing, and make phe_idxs and cov_idxs correspond to the new order.
            std::map<int32_t, int32_t> unsrt2srt_idx;
            for(int32_t i=0; i < (int32_t)geno_idxs.size(); ++i) {
                unsrt2srt_idx[geno_idxs[i]] = i; // map from unsorted index to sorted index
            }
            std::vector<int32_t> reordered_geno_idxs;
            std::vector<int32_t> reordered_phe_idxs;
            std::vector<int32_t> reordered_cov_idxs;
            std::map<int32_t, int32_t>::iterator it;
            for(it = unsrt2srt_idx.begin(); it != unsrt2srt_idx.end(); ++it) {
                reordered_geno_idxs.push_back(geno_idxs[it->second]);
                reordered_phe_idxs.push_back(phe_idxs[it->second]);
                reordered_cov_idxs.push_back(cov_idxs[it->second]);
            }
            geno_idxs = reordered_geno_idxs;
            phe_idxs = reordered_phe_idxs;
            cov_idxs = reordered_cov_idxs;

            pr.subset_sample_indices(geno_idxs, false); 
            // for(int32_t i=0; i < (int32_t)geno_idxs.size(); ++i) {
            //     notice("%d\t%d\t%d\t%d", i, geno_idxs[i], phe_idxs[i], cov_idxs[i]);
            // }
        }
        else { // when sample list is not provided, we first need to identify overlapping sample IDs
            const std::vector<plink_samp_t>& geno_samps = pr.get_all_samples();
            for(int32_t i = 0; i < (int32_t)geno_samps.size(); ++i) {
                const std::string& id = geno_samps[i].indID;
                if ( phe_id2idx.find(id) != phe_id2idx.end() ) { // exists in pheno
                    if ( cov_id2idx.find(id) != cov_id2idx.end() ) { // exists in covariates
                        cov_idxs.push_back(cov_id2idx[id]); // store the index of the sample ID
                        phe_idxs.push_back(phe_id2idx[id]);
                        geno_idxs.push_back(i + 1); // 1-based index
                    }
                }
            }
            pr.subset_sample_indices(geno_idxs, false); // subset the genotype samples to only those that are in both phenotype and covariate files
        }

        notice("Identified %zu overlapping sample IDs between genotype, phenotype and covariate files", geno_idxs.size());

        // load the covariate data
        int32_t nrow_est = 100;
        cov_mat.resize(nrow_est, cov_idxs.size());
        int32_t irow = 0;
        while( tr_cov.read_line() ) {
            cov_trait_ids.push_back(tr_cov.str_field_at(icol_cov_id)); // first field is the covariate name
            if ( irow >= nrow_est ) {
                cov_mat.conservativeResize(nrow_est * 2, cov_idxs.size());
                nrow_est *= 2;
            }
            for(int32_t i = 0; i < cov_idxs.size(); ++i) {
                if ( strcmp(tr_cov.str_field_at(cov_idxs[i] + offset_cov), "NA") == 0 ) 
                    error("Missing value is not allowed in the covariate file %s at row %d, column %d", covf.c_str(), irow + 1, cov_idxs[i] + offset_cov + 1);
                cov_mat(irow, i) = tr_cov.double_field_at(cov_idxs[i]+offset_cov);
            }
            ++irow;
        }
        cov_mat.conservativeResize(irow, cov_idxs.size());

        notice("Loaded %d covariates for %zu samples", irow, (int32_t)cov_idxs.size());
    }
    else { // no covariates are provided
        // now determine cov_idxs, phe_idxs, and geno_idx 
        if ( sample_list_provided ) { // when sample list is provided, make sure that all samples exist
            fill_list_indices(list_pheno_sample_ids, phe_id2idx, phe_idxs);
            fill_list_indices(list_geno_sample_ids, psam_samp2idx, geno_idxs);

            // need to make geno_idxs strictly increasing, and make phe_idxs and cov_idxs correspond to the new order.
            std::map<int32_t, int32_t> unsrt2srt_idx;
            for(int32_t i=0; i < (int32_t)geno_idxs.size(); ++i) {
                unsrt2srt_idx[geno_idxs[i]] = i; // map from unsorted index to sorted index
            }
            std::vector<int32_t> reordered_geno_idxs;
            std::vector<int32_t> reordered_phe_idxs;
            std::map<int32_t, int32_t>::iterator it;
            for(it = unsrt2srt_idx.begin(); it != unsrt2srt_idx.end(); ++it) {
                reordered_geno_idxs.push_back(geno_idxs[it->second]);
                reordered_phe_idxs.push_back(phe_idxs[it->second]);
            }
            geno_idxs = reordered_geno_idxs;
            phe_idxs = reordered_phe_idxs;

            pr.subset_sample_indices(geno_idxs, false); 
        }
        else { // when sample list is not provided, we first need to identify overlapping sample IDs
            const std::vector<plink_samp_t>& geno_samps = pr.get_all_samples();
            for(int32_t i = 0; i < (int32_t)geno_samps.size(); ++i) {
                const std::string& id = geno_samps[i].indID;
                if ( phe_id2idx.find(id) != phe_id2idx.end() ) { // exists in pheno
                    phe_idxs.push_back(phe_id2idx[id]);
                    geno_idxs.push_back(i + 1); // 1-based index
                }
            }
            pr.subset_sample_indices(geno_idxs, false); // subset the genotype samples to only those that are in both phenotype and covariate files
        }
        notice("Identified %zu overlapping sample IDs between genotype and phenotype files", geno_idxs.size());
    }

    // notice("phe_idxs.size() = %zu, cov_idxs.size() = %zu, geno_idxs.size() = %zu", phe_idxs.size(), cov_idxs.size(), geno_idxs.size());
    // for(int32_t i = 0; i < phe_idxs.size(); ++i) {
    //     const std::vector<plink_samp_t>& geno_samps = pr.get_all_samples();
    //     notice("%d\t%d\t%s\t%d\t%s\t%d\t%s", i, geno_idxs[i], geno_samps[geno_idxs[i]-1].indID.c_str(),
    //            phe_idxs[i], phe_sample_ids[phe_idxs[i]].c_str(),
    //            cov_idxs.empty() ? -1 : cov_idxs[i], cov_idxs.empty() ? "NA" : cov_sample_ids[cov_idxs[i]].c_str());
    // }

    std::vector<std::string> phe_trait_ids;
    std::map<std::string, int32_t> phe_trait2idx; // map phenotype ID to idx
    std::map<std::string, int32_t>::iterator phe_trait2idx_it;
    int32_t nrow_est = 100;
    int32_t irow = 0;
    Eigen::MatrixXd phe_mat;
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> phe_mask;
    uint64_t n_phe_missing = 0;
    phe_mat.resize(nrow_est, phe_idxs.size());
    phe_mask.resize(nrow_est, phe_idxs.size());
    while( tr_phe.read_line() ) {
        phe_trait_ids.push_back(tr_phe.str_field_at(icol_pheno_id-1)); // the phenotype ID
        phe_trait2idx[phe_trait_ids.back()] = irow; // store the index of the phenotype ID
        if ( irow >= nrow_est ) {
            phe_mat.conservativeResize(nrow_est * 2, phe_idxs.size());
            phe_mask.conservativeResize(nrow_est * 2, phe_idxs.size());
            nrow_est *= 2;
        }
        for(int32_t i = 0; i < phe_idxs.size(); ++i) {
            if ( strcmp(tr_phe.str_field_at(phe_idxs[i] + offset_pheno), "NA") == 0 ) {
                phe_mask(irow, i) = false;
                phe_mat(irow, i) = 0;
                ++n_phe_missing; // count the number of missing values
            }
            else {
                phe_mask(irow, i) = true;
                phe_mat(irow, i) = tr_phe.double_field_at(phe_idxs[i]+offset_pheno);
            }
        }
        ++irow;
    }
    phe_mat.conservativeResize(irow, phe_idxs.size());
    phe_mask.conservativeResize(irow, phe_idxs.size());
    notice("Loaded %d phenotypes for %zu samples with %d missing values. Dimension : %d x %d", irow, (int32_t)phe_idxs.size(), n_phe_missing, phe_mat.rows(), phe_mat.cols());

    Eigen::MatrixXd phe_mat_adj;
    if ( covf.empty() ) {
        notice("No covariates provided. Centering the phenotypes by sample means");
        // center the phenotypes by sample means
        phe_mat_adj = phe_mat;
        center_rows(phe_mat_adj);
    }
    else {
        notice("Adjusting the phenotypes by covariates");
        phe_mat_adj = bulk_adjust_for_covariates(phe_mat, cov_mat);
    }
    notice("Finished adjusting the phenotype matrix");

    // read the variant-trait pair lists
    tsv_reader tr_pair(pairf.c_str());
    std::map<int32_t, std::set<cpra_t> > pair_map;
    std::map<int32_t, std::set<cpra_t> >::iterator pair_map_it;
    std::map<int32_t, std::set<cbe_t> > region_map;
    std::map<int32_t, std::set<cbe_t> >::iterator region_map_it;

    notice("Reading the trait-variant/region pairs from %s", pairf.c_str());
    int32_t n_pairs = 0;

    // open the output file gz or plain based on the extension
    htsFile* wf = hts_open(outf.c_str(), outf.substr(outf.length() - 3).compare(".gz") == 0 ? "wz" : "w");
    if ( wf == NULL ) {
        error("Cannot open output file %s for writing", outf.c_str());
    }
    // write the header line
    hprintf(wf, "#TRAIT\tCHROM\tPOS\tID\tREF\tALT\tAF\tINFO\tN\tBETA\tSE\tTSTAT\tLOG10P");
    if ( !skip_rint ) {
        hprintf(wf, "\tBETA_RINT\tSE_RINT\tTSTAT_RINT\tLOG10P_RINT");
    }
    hprintf(wf, "\n");

    int32_t icol_pair_trait = -1;
    int32_t icol_pair_variant = -1;
    int32_t icol_pair_region = -1;
    bool region_mode = false;
    while( tr_pair.read_line() ) {
        if ( icol_pair_trait < 0 ) {
            for(int32_t i = 0; i < tr_pair.nfields; ++i) {
                if ( colname_pair_trait.compare(tr_pair.str_field_at(i)) == 0 ) {
                    icol_pair_trait = i;
                }
                else if ( colname_pair_variant.compare(tr_pair.str_field_at(i)) == 0 ) {
                    icol_pair_variant = i;
                    region_mode = false;
                }
                else if ( colname_pair_region.compare(tr_pair.str_field_at(i)) == 0 ) {
                    icol_pair_region = i;
                    region_mode = true;
                }
            }
            if ( icol_pair_trait < 0 ) {
                error("Cannot find column %s in the pair file %s", colname_pair_trait.c_str(), pairf.c_str());
            }
            if ( icol_pair_variant < 0 && icol_pair_region < 0 ) {
                error("Either %s or %s must exist in the pair file %s", colname_pair_variant.c_str(), colname_pair_region.c_str(), pairf.c_str());
            }
            if ( icol_pair_variant >= 0 && icol_pair_region >= 0 ) {
                error("Both %s and %s cannot exist in the pair file %s. Please provide either variant ID or region string", colname_pair_variant.c_str(), colname_pair_region.c_str(), pairf.c_str());
            }
        }
        else {
            if ( icol_pair_variant >= 0 ) {
                const char* var_id = tr_pair.str_field_at(icol_pair_variant);
                const char* phe_id = tr_pair.str_field_at(icol_pair_trait);
                phe_trait2idx_it = phe_trait2idx.find(phe_id);
                if ( phe_trait2idx_it == phe_trait2idx.end() ) {
                    notice("Skipping phenotype %s, which is not observed in %s", phe_id, phef.c_str());
                }
                else {
                    int32_t phe_idx = phe_trait2idx_it->second;
                    cpra_t cpra(var_id);
                    if ( pair_map[phe_idx].insert(cpra).second ) {
                        ++n_pairs;
                    }
                    else {
                        notice("Skipping duplicate pair %s\t%s", var_id, phe_id);
                    }
                }
            }
            else {
                const char* phe_id = tr_pair.str_field_at(icol_pair_trait);
                const char* region_str = tr_pair.str_field_at(icol_pair_region);
                phe_trait2idx_it = phe_trait2idx.find(phe_id);
                if ( phe_trait2idx_it == phe_trait2idx.end() ) {
                    notice("Skipping phenotype %s, which is not observed in %s", phe_id, phef.c_str());
                }
                else {
                    int32_t phe_idx = phe_trait2idx_it->second;
                    cbe_t cbe(region_str);
                    if ( region_map[phe_idx].insert(cbe).second ) {
                        ++n_pairs;
                    }
                    else {
                        notice("Skipping duplicate pair %s\t%s", region_str, phe_id);
                    }
                }
                //error("region-based pair association is not supported yet. Please provide the variant ID in the pair file");
            }
        }
    }
    notice("Finished reading %d trait-%s pairs across %zu traits", n_pairs, region_mode ? "region" : "variant", pair_map.size());

    if ( region_mode ) {
        for(region_map_it = region_map.begin(); region_map_it != region_map.end(); ++region_map_it) {
            int32_t phe_idx = region_map_it->first;
            const std::set<cbe_t>& regions = region_map_it->second;
            if ( regions.empty() ) continue;

            notice("Processing %zu regions for phenotype %s", regions.size(), phe_trait_ids[phe_idx].c_str());
            int32_t n_row = (int32_t)geno_idxs.size();
            Eigen::VectorXd phe_vec(n_row);
            phe_vec = phe_mat_adj.row(phe_idx); // get the phenotype vector for the current trait

            // process one region at a time. The variant may overlap, but let's not worry about it for now
            for(std::set<cbe_t>::const_iterator regions_it = regions.begin(); regions_it != regions.end(); ++regions_it) {
                const cbe_t& region = *regions_it;
                notice("Processing region %s for phenotype %s", region.to_string().c_str(), phe_trait_ids[phe_idx].c_str());
                Eigen::Vector<bool, Eigen::Dynamic> phe_mask_vec = phe_mask.row(phe_idx); // get the missing values for the current trait

                bool phe_has_missing = !phe_mask_vec.all();
                Eigen::VectorXd phe_rint_vec;
                if ( ! skip_rint ) {
                    // perform rank-based inverse normal transformation
                    notice("Performing rank-based inverse normal transformation for phenotype %s", phe_trait_ids[phe_idx].c_str());
                    if ( phe_has_missing ) {
                        phe_rint_vec = rint_with_missing(phe_vec, phe_mask.row(phe_idx));
                    }
                    else {
                        phe_rint_vec = rint_without_missing(phe_vec);
                    }
                    if ( phe_rint_vec.size() != n_row ) {
                        error("Rank-based inverse normal transformation failed for phenotype %s", phe_trait_ids[phe_idx].c_str());
                    }
                }
                else {
                    notice("Skipping rank-based inverse normal transformation for phenotype %s", phe_trait_ids[phe_idx].c_str());
                }

                if ( pr.read_pos(region.chrom.c_str(), region.beg1) ) { // variant exists, start reading the region
                    std::vector<cpra_t> v_cpra;
                    std::vector<int32_t> ans;
                    std::vector<double> acs;
                    std::vector<double> infos;
                    const std::vector<int32_t>& int_buf = pr.get_int_buf();
                    const double* dbl_buf = pr.get_dbl_buf();
                    int32_t n_col_est = 10;
                    Eigen::MatrixXd geno_mat(n_row, n_col_est);
                    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> geno_mask(n_row, n_col_est);
                    uint32_t n_geno_missing = 0;
                    int32_t icol = 0;

                    do {
                        const plink_var_t& var = pr.get_current_variant();
                        if ( var.pos > region.end0 ) {
                            break;
                        }
                        std::string cpra_s(var.to_string());
                        pr.get_genos();

                        if ( icol >= n_col_est ) {
                            geno_mat.conservativeResize(n_row, n_col_est * 2);
                            geno_mask.conservativeResize(n_row, n_col_est * 2);
                            n_col_est *= 2;
                        }
                        int32_t an = 0;
                        double ac = 0;
                        if ( pr.is_dosage_present() ) {
                            if ( dbl_buf == NULL) {
                                dbl_buf = pr.get_dbl_buf();
                            }
                            double sumsq = 0;
                            for(int32_t i =0; i < n_row; ++i) {
                                double ds = 2.0 - dbl_buf[i];
                                geno_mat(i, icol) = ds;
                                geno_mask(i, icol) = true; // not missing
                                ac += ds;
                                an += 2;
                                sumsq += (ds * ds);
                            } 
                            if ( ac == 0 || an == ac ) {
                                // skip monomorphic variants
                                continue;
                            }
                            // E(Var(g)) = 2 * af * (1-af) = 2 * ac * ( an - ac ) / an / an
                            // Var(g) = sumsq / (an / 2) - 4 * ac * ac / an / an
                            // ratio = ( sumsq / 2 * an - 4 * ac * ac ) / ( 2 * ac * ( an - ac ) )
                            double info = ( sumsq / 2.0 * an - 4.0 * ac * ac ) / ( 2.0 * ac * ( an - ac ) );
                            infos.push_back(info);
                        }
                        else {
                            int32_t gcs[3] = {0, 0, 0};
                            for(int32_t i = 0; i < n_row; ++i) {
                                switch(int_buf[i]) { // make sure to convert 1-based index to 0-based
                                case 0:
                                    an += 2;
                                    ac += 2;
                                    ++gcs[2];
                                    break;
                                case 1:
                                    an += 2;
                                    ++ac;
                                    ++gcs[1];
                                    break;
                                case 2:
                                    an += 2;
                                    ++gcs[0];
                                    break;
                                }
                            }
                            if ( ac == 0 || an == ac ) {
                                // skip monomorphic variants
                                continue;
                            }
                            double mean = (double)ac / (double)an * 2.0;
                            // Exp(Var(g)) = af * (1-af) * 2 = mean * (2 - mean) / 2;
                            // Var(g) = EX^2 - EX^2 = (4 * n_2 + 1 * n_1)/n - mean^2
                            double info = ((4.0 * gcs[2] + gcs[1])/(an/2.0) - mean * mean) / (mean * (2.0 - mean) / 2.0);
                            infos.push_back(info);
                            //notice("mean = %.5g, an = %d, ac = %d", mean, an, ac);
                            for(int32_t i = 0; i < n_row; ++i) {
                                switch(int_buf[i]) {
                                case 0:
                                    geno_mat(i, icol) = 2.0 - mean; // homalt
                                    geno_mask(i, icol) = true; // not missing
                                    break;
                                case 1:
                                    geno_mat(i, icol) = 1.0 - mean; // het
                                    geno_mask(i, icol) = true; // not missing
                                    break;
                                case 2:
                                    geno_mat(i, icol) = 0.0 - mean; // homref
                                    geno_mask(i, icol) = true; // not missing
                                    break;
                                default:
                                    geno_mat(i, icol) = 0; // missing - mean imputation
                                    geno_mask(i, icol) = false; // missing
                                    ++n_geno_missing;
                                    break;
                                }
                            }
                        }
                        v_cpra.push_back(cpra_t(cpra_s.c_str())); 
                        acs.push_back(ac);
                        ans.push_back(an);
                        ++icol;

                        if ( icol >= max_chunk_vars ) { // process the current chunk perform association mapping
                            notice("Processing %d variants for phenotype %s", icol, phe_trait_ids[phe_idx].c_str());
                            
                            geno_mat.conservativeResize(n_row, icol);
                            geno_mask.conservativeResize(n_row, icol);

                            assoc_single_trait( 
                                wf, // output file handle
                                phe_trait_ids[phe_idx].c_str(), // phenotype ID
                                phe_vec,      // phenotype vector
                                phe_rint_vec, // rinted phenotype vector 
                                phe_mask_vec, // phenotype mask vector
                                phe_has_missing, // if the phenotype has missing values
                                skip_rint, // skip rank-based inverse normal transformation
                                geno_mat,    // genotype matrix
                                geno_mask,   // genotype mask matrix
                                n_geno_missing > 0, // if the genotype has missing values
                                v_cpra,      // variant pairs
                                ans,         // allele counts
                                acs,         // allele counts
                                infos       // information values
                            );
                            icol = 0;
                            n_col_est = 10;
                            ans.clear();
                            acs.clear();
                            infos.clear();
                            v_cpra.clear();
                            geno_mat.resize(n_row, n_col_est);
                            geno_mask.resize(n_row, n_col_est);
                            n_geno_missing = 0;
                        }
                    }
                    while ( pr.read_pivar() );

                    if ( icol > 0 ) {
                        if ( icol < n_col_est ) {
                            geno_mat.conservativeResize(n_row, icol);
                            geno_mask.conservativeResize(n_row, icol);
                        }

                        notice("Processing %d variants for phenotype %s", icol, phe_trait_ids[phe_idx].c_str());
                        
                        geno_mat.conservativeResize(n_row, icol);
                        geno_mask.conservativeResize(n_row, icol);

                        assoc_single_trait( 
                            wf, // output file handle
                            phe_trait_ids[phe_idx].c_str(), // phenotype ID
                            phe_vec,      // phenotype vector
                            phe_rint_vec, // rinted phenotype vector 
                            phe_mask_vec, // phenotype mask vector
                            phe_has_missing, // if the phenotype has missing values
                            skip_rint, // skip rank-based inverse normal transformation
                            geno_mat,    // genotype matrix
                            geno_mask,   // genotype mask matrix
                            n_geno_missing > 0, // if the genotype has missing values
                            v_cpra,      // variant pairs
                            ans,         // allele counts
                            acs,         // allele counts
                            infos       // information values
                        );
                        icol = 0;
                        ans.clear();
                        acs.clear();
                        infos.clear();
                        v_cpra.clear();
                        geno_mat.resize(n_row, n_col_est);
                        geno_mask.resize(n_row, n_col_est);
                        n_geno_missing = 0;
                    }
                }
            }
        }   
    }
    else {
        for(pair_map_it = pair_map.begin(); pair_map_it != pair_map.end(); ++pair_map_it) {
            int32_t phe_idx = pair_map_it->first;
            const std::set<cpra_t>& pairs = pair_map_it->second;
            if ( pairs.empty() ) continue;

            notice("Processing %zu pairs for phenotype %s", pairs.size(), phe_trait_ids[phe_idx].c_str());
            int32_t n_col_est = 10;
            int32_t n_row = (int32_t)geno_idxs.size();
            Eigen::MatrixXd geno_mat(n_row, n_col_est);
            Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> geno_mask(n_row, n_col_est);
            uint32_t n_geno_missing = 0;
            int32_t icol = 0;
            Eigen::VectorXd phe_vec(n_row);
            phe_vec = phe_mat_adj.row(phe_idx); // get the phenotype vector for the current trait

            Eigen::Vector<bool, Eigen::Dynamic> phe_mask_vec = phe_mask.row(phe_idx); // get the missing values for the current trait

            bool phe_has_missing = !phe_mask_vec.all();
            Eigen::VectorXd phe_rint_vec;
            if ( ! skip_rint ) {
                // perform rank-based inverse normal transformation
                notice("Performing rank-based inverse normal transformation for phenotype %s", phe_trait_ids[phe_idx].c_str());
                if ( phe_has_missing ) {
                    phe_rint_vec = rint_with_missing(phe_vec, phe_mask.row(phe_idx));
                }
                else {
                    phe_rint_vec = rint_without_missing(phe_vec);
                }
                if ( phe_rint_vec.size() != n_row ) {
                    error("Rank-based inverse normal transformation failed for phenotype %s", phe_trait_ids[phe_idx].c_str());
                }
            }
            else {
                notice("Skipping rank-based inverse normal transformation for phenotype %s", phe_trait_ids[phe_idx].c_str());
            }

            // for(int32_t i = 0; i < n_row; ++i) {
            //     //phe_vec(i) = phe_mat_adj(phe_idx, phe_idxs[i]);
            //     phe_vec(i) = phe_mat_adj(phe_idx, i);
            // }
            // read the genotypes for each pair
            std::vector<cpra_t> v_cpra;
            std::vector<int32_t> ans;
            std::vector<double> acs;
            std::vector<double> infos;
            const std::vector<int32_t>& int_buf = pr.get_int_buf();
            const double* dbl_buf = pr.get_dbl_buf();
            // notice("int_buf.size() = %zu", int_buf.size());
            for(std::set<cpra_t>::iterator it = pairs.begin(); it != pairs.end(); ++it) {
                std::string cpra_s(it->to_string());
                //notice("Reading variant %s", cpra_s.c_str());
                if ( !pr.read_pivar(cpra_s.c_str()) ) {
                    notice("Skipping pair %s\t%s, which is not found in the pvar file", it->to_string().c_str(), phe_trait_ids[phe_idx].c_str());
                    continue;
                }   
                notice("Reading genotypes for variant %s", cpra_s.c_str());
                pr.get_genos();
                notice("Finished reading genotypes for variant %s", cpra_s.c_str());
                // construct the input for association analysis
                if ( icol >= n_col_est ) {
                    geno_mat.conservativeResize(n_row, n_col_est * 2);
                    geno_mask.conservativeResize(n_row, n_col_est * 2);
                    n_col_est *= 2;
                }
                int32_t an = 0;
                double ac = 0;
                if ( pr.is_dosage_present() ) {
                    if ( dbl_buf == NULL) {
                        dbl_buf = pr.get_dbl_buf();
                    }
                    double sumsq = 0;
                    for(int32_t i =0; i < n_row; ++i) {
                        //notice("Dosage[%d] = %.5g", i, dbl_buf[i]);
                        double ds = 2.0 - dbl_buf[i];
                        geno_mat(i, icol) = ds;
                        geno_mask(i, icol) = true; // not missing
                        ac += ds;
                        an += 2;
                        sumsq += (ds * ds);
                    } 
                    if ( ac == 0 || an == ac ) {
                        // skip monomorphic variants
                        continue;
                    }

                    // E(Var(g)) = 2 * af * (1-af) = 2 * ac * ( an - ac ) / an / an
                    // Var(g) = sumsq / (an / 2) - 4 * ac * ac / an / an
                    // ratio = ( sumsq / 2 * an - 4 * ac * ac ) / ( 2 * ac * ( an - ac ) )
                    double info = ( sumsq / 2 * an - 4 * ac * ac ) / ( 2 * ac * ( an - ac ) );
                    infos.push_back(info);
                }
                else {
                    int32_t gcs[3] = {0, 0, 0};
                    for(int32_t i = 0; i < n_row; ++i) {
                        switch(int_buf[i]) { // make sure to convert 1-based index to 0-based
                        case 0:
                            an += 2;
                            ac += 2;
                            ++gcs[2];
                            break;
                        case 1:
                            an += 2;
                            ++ac;
                            ++gcs[1];
                            break;
                        case 2:
                            an += 2;
                            ++gcs[0];
                            break;
                        }
                    }

                    if ( ac == 0 || an == ac ) {
                        // skip monomorphic variants
                        continue;
                    }

                    double mean = (double)ac / (double)an * 2.0;
                    // Exp(Var(g)) = af * (1-af) * 2 = mean * (2 - mean) / 2;
                    // Var(g) = EX^2 - EX^2 = (4 * n_2 + 1 * n_1)/n - mean^2
                    double info = ((4.0 * gcs[2] + gcs[1])/(an/2.0) - mean * mean) / (mean * (2.0 - mean) / 2.0);
                    infos.push_back(info);
                    //notice("mean = %.5g, an = %d, ac = %d", mean, an, ac);
                    for(int32_t i = 0; i < n_row; ++i) {
                        switch(int_buf[i]) {
                        case 0:
                            geno_mat(i, icol) = 2.0 - mean; // homalt
                            geno_mask(i, icol) = true; // not missing
                            break;
                        case 1:
                            geno_mat(i, icol) = 1.0 - mean; // het
                            geno_mask(i, icol) = true; // not missing
                            break;
                        case 2:
                            geno_mat(i, icol) = 0.0 - mean; // homref
                            geno_mask(i, icol) = true; // not missing
                            break;
                        default:
                            geno_mat(i, icol) = 0; // missing - mean imputation
                            geno_mask(i, icol) = false; // missing
                            ++n_geno_missing;
                            break;
                        }
                    }
                }
                v_cpra.push_back(*it); 
                acs.push_back(ac);
                ans.push_back(an);
                ++icol;
            }

            if ( icol > 0 ) {
                if ( icol < n_col_est ) {
                    geno_mat.conservativeResize(n_row, icol);
                    geno_mask.conservativeResize(n_row, icol);
                }

                assoc_single_trait( 
                    wf, // output file handle
                    phe_trait_ids[phe_idx].c_str(), // phenotype ID
                    phe_vec,      // phenotype vector
                    phe_rint_vec, // rinted phenotype vector 
                    phe_mask_vec, // phenotype mask vector
                    phe_has_missing, // if the phenotype has missing values
                    skip_rint, // skip rank-based inverse normal transformation
                    geno_mat,    // genotype matrix
                    geno_mask,   // genotype mask matrix
                    n_geno_missing > 0, // if the genotype has missing values
                    v_cpra,      // variant pairs
                    ans,         // allele counts
                    acs,         // allele counts
                    infos       // information values
                );
            }

            // std::vector<slr_sumstat_t> sumstats;
            // std::vector<slr_sumstat_t> sumstats_rint;
            // if ( ( n_geno_missing == 0 ) && !phe_has_missing ) {
            //     simple_linear_regression_without_missing(phe_vec, geno_mat, sumstats);
            //     if ( ! skip_rint ) {
            //         simple_linear_regression_without_missing(phe_rint_vec, geno_mat, sumstats_rint);
            //     }
            // }
            // else {
            //     simple_linear_regression_with_missing(phe_vec, phe_mask_vec, geno_mat, geno_mask, sumstats);
            //     if ( ! skip_rint ) {
            //         simple_linear_regression_with_missing(phe_rint_vec, phe_mask_vec, geno_mat, geno_mask, sumstats_rint);
            //     }
            // }

            // // print the results
            // for(int32_t i=0; i < (int32_t)v_cpra.size(); ++i) {
            //     const slr_sumstat_t& ss = sumstats[i];
            //     hprintf(wf, "%s\t%s\t%d\t%s\t%s\t%s\t%.5g\t%.5g\t%d\tADD\t%.6g\t%.6g\t%.6g\t%.6g",
            //         phe_trait_ids[phe_idx].c_str(), // TRAIT
            //         v_cpra[i].chrom.c_str(), // CHROM
            //         v_cpra[i].pos,           // POS
            //         v_cpra[i].to_string().c_str(), // ID
            //         v_cpra[i].ref.c_str(),  // REF
            //         v_cpra[i].alts.c_str(), // ALT
            //         (double)acs[i] / (double)ans[i], // AF
            //         infos[i],  // INFO - placeholder, not calculated
            //         ss.n_obs, // N - number of samples
            //         ss.beta, // BETA
            //         ss.se,   // SE
            //         ss.tstat, // TSTAT
            //         ss.log10p); // LOG10P
            //     if ( ! skip_rint ) {
            //         const slr_sumstat_t& ss_rint = sumstats_rint[i];
            //         hprintf(wf, "\t%.6g\t%.6g\t%.6g\t%.6g", // BETA_RINT, SE_RINT, TSTAT_RINT, LOG10P_RINT
            //             ss_rint.beta, ss_rint.se, ss_rint.tstat, ss_rint.log10p);
            //     }
            //     hprintf(wf, "\n");
            // }
        }
    }
    hts_close(wf); // close the output file
    notice("Analysis finished");
    return 0;
}