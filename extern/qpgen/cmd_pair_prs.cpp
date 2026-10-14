#include "qgenlib/params.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_utils.h"
#include "qgenlib/phred_helper.h"
#include "qgenlib/hts_utils.h"
#include "qpgen.h"
#include "qpgen_utils.h"
#include "Eigen/Dense"
#include <cmath>

int32_t cmd_pair_prs(int32_t argc, char **argv)
{
    std::string pgenf;      // PLINK2 genotype file
    std::string psamf;      // PLINK2 sample file
    std::string pivarf;     // PLINK2 indexed pvar file, compatible with PgenIdxReader
    std::string listf;      // CHROM BEG END PGEN PSAM PIVAR containing the list of region-specific BED files
    std::string pairf;      // file containing TRAIT VARIANT BETA SE
    std::string outf;       // output containing PRS values
    std::string samplef;    // sample IDs to be used
    int32_t jump_thres_bp = 1000000; 
    int32_t max_chunk_vars = 1000; // Maximum number of variants to store at once in memory
    int32_t icol_pivar_idx = 9;
    std::string colname_pair_trait("trait"); // column name for the trait ID in the pair file
    std::string colname_pair_variant("variant"); // column name of the variant ID in the pair file
    std::string colname_pair_beta("beta");   // column name of the effect sizes
    std::string colname_pair_se("se");       // column name of the standard errors

    bool space_delim = false; // space- instead of tab-delimited input tables

    paramList pl;

    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Options", NULL)
    LONG_STRING_PARAM("pgen", &pgenf, "Input PLINK2 genotype file")
    LONG_STRING_PARAM("psam", &psamf, "Input PLINK2 sample file")
    LONG_STRING_PARAM("pivar", &pivarf, "Input PLINK2 index pvar file (bgzipped and tabix)")
    LONG_STRING_PARAM("pgen-list", &listf, "Input file containing CHROM BEG END PGEN PSAM PIVAR")
    LONG_STRING_PARAM("pairs", &pairf, "Input file containing TRAIT VARIANT BETA SE summary statistics")
    LONG_STRING_PARAM("sample", &samplef, "Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output prefix")

    LONG_PARAM_GROUP("Auxiliary options", NULL)
    LONG_INT_PARAM("jump-thres-bp", &jump_thres_bp, "Jump threshold in base pairs for the variant index (default: 1000000)")
    LONG_INT_PARAM("max-chunk-vars", &max_chunk_vars, "Maximum number of variants to store at once in memory (default: 1000)")
    LONG_INT_PARAM("icol-pivar-idx", &icol_pivar_idx, "1-based column index for the variant ID in the pvar file (default: 9)")
    LONG_STRING_PARAM("colname-pair-trait", &colname_pair_trait, "Column name of the trait ID in the pair file")
    LONG_STRING_PARAM("colname-pair-variant", &colname_pair_variant, "Column name of the variant ID in the pair file")
    LONG_STRING_PARAM("colname-pair-beta", &colname_pair_beta, "Column name of the effect sizes in the pair file")
    LONG_STRING_PARAM("colname-pair-se", &colname_pair_se, "Column name of the standard errors in the pair file")
    LONG_PARAM("space-delim", &space_delim, QPGEN_SPACE_DELIM_HELP)
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();
    const char table_delim = table_delimiter(space_delim);

    // open the Plink file
    MultiPgenIdxReader mpr;
    if ( !listf.empty() ) { // list is provided
        //notice("foo");
        if ( pgenf.empty() && pivarf.empty() && psamf.empty() ) {
            if ( !mpr.prep_pgen_list(listf.c_str()) ) {
                error("Failed to prepare pgen files with the following list file: %s", listf.c_str());
            }
        }
        else {
            error("When --list is provided, --pgen, --pivar, and --psam should not be provided");
        }
    }
    else {
        //notice("bar");
        if ( pgenf.empty() || pivarf.empty() || psamf.empty() ) {
            error("When --list is not provided, --pgen, --pivar, and --psam should be provided");
        }
        if ( !mpr.set_single_chunk_pgen(pgenf.c_str(), pivarf.c_str(), psamf.c_str()) ) {
            error("Failed to add pgen files with the following files:\n%s\n%s\n%s", pgenf.c_str(), pivarf.c_str(), psamf.c_str());
        }
    }
    mpr.set_jump_thres_bp(jump_thres_bp);
    mpr.set_icol_pivar_idx(icol_pivar_idx - 1); // convert to 0-based index

    // Rules for determining the sample IDs to analyze:
    // 1. If --sample is probided, use exactly the sample IDs provided in the input file. If any sample is missing, report an error
    //    The order of samples should be the same as the input file
    // 2. If --sample is not provided, use the overlapping sample IDs
    bool sample_list_provided = !samplef.empty();  // indicates if the sample file is provided
    if ( sample_list_provided ) { 
        std::vector<std::string> list_geno_sample_ids;          // list of sample IDs in the genotype file
        std::map<std::string, int32_t> list_geno_sample2idx;
        tsv_reader tr_sample(samplef.c_str());
        tr_sample.delimiter = table_delim;
        while(tr_sample.read_line()) {
            if ( tr_sample.nfields > 1 ) {
                error("The sample file %s should have at most one column, containing only IID", samplef.c_str());
            }
            const char* sample_id = tr_sample.str_field_at(0);
            if ( list_geno_sample2idx.find(sample_id) != list_geno_sample2idx.end() ) {
                error("Sample ID %s is duplicated in the sample file %s", sample_id, samplef.c_str());
            }
            list_geno_sample2idx[sample_id] = list_geno_sample_ids.size();
            list_geno_sample_ids.push_back(sample_id);
        }
        tr_sample.close();       

        // subset the samples
        mpr.subset_sample_ids(list_geno_sample_ids);
    }
    //int32_t n_geno_samples = mpr.get_loaded_sample_count();

    // read the variant-trait pair lists
    tsv_reader tr_pair(pairf.c_str());
    tr_pair.delimiter = table_delim;
    std::map<int32_t, std::map<cpra_t, std::pair<double, double> > > pair_map;
    std::map<int32_t, std::map<cpra_t, std::pair<double, double> > >::iterator pair_map_it;
    std::map<cpra_t, std::vector<int32_t> > var2traits;
    std::vector<std::string> trait_ids;
    std::map<std::string, int32_t> trait_id2idx;

    notice("Reading the trait-variant pair sumstats from %s", pairf.c_str());
    int32_t n_pairs = 0;

    int32_t icol_pair_trait = -1;
    int32_t icol_pair_variant = -1;
    int32_t icol_pair_beta = -1;
    int32_t icol_pair_se = -1;
    bool region_mode = false;
    while( tr_pair.read_line() ) {
        if ( icol_pair_trait < 0 ) {
            for(int32_t i = 0; i < tr_pair.nfields; ++i) {
                if ( colname_pair_trait.compare(tr_pair.str_field_at(i)) == 0 ) {
                    icol_pair_trait = i;
                }
                else if ( colname_pair_variant.compare(tr_pair.str_field_at(i)) == 0 ) {
                    icol_pair_variant = i;
                }
                else if ( colname_pair_beta.compare(tr_pair.str_field_at(i)) == 0 ) {
                    icol_pair_beta = i;
                }
                else if ( colname_pair_se.compare(tr_pair.str_field_at(i)) == 0 ) {
                    icol_pair_se = i;
                }
            }
            if ( icol_pair_trait < 0 ) {
                error("Cannot find column %s in the pair file %s", colname_pair_trait.c_str(), pairf.c_str());
            }
            if ( icol_pair_variant < 0 ) {
                error("Cannot find column %s in the pair file %s", colname_pair_variant.c_str(), pairf.c_str());
            }
            if ( icol_pair_beta < 0 ) {
                error("Cannot find column %s in the pair file %s", colname_pair_beta.c_str(), pairf.c_str());
            }
        }
        else {
            const char* var_id = tr_pair.str_field_at(icol_pair_variant);
            const char* phe_id = tr_pair.str_field_at(icol_pair_trait);
            double beta = tr_pair.double_field_at(icol_pair_beta);
            double se = icol_pair_se < 0 ? 0 :tr_pair.double_field_at(icol_pair_se);

            // check if trait already exists
            int32_t phe_idx = -1;
            if ( trait_id2idx.find(phe_id) == trait_id2idx.end() ) {
                trait_id2idx[phe_id] = trait_ids.size();
                trait_ids.push_back(phe_id);
                phe_idx = trait_ids.size() - 1;
            } 
            else {
                phe_idx = trait_id2idx[phe_id];
            }

            cpra_t cpra(var_id);
            var2traits[cpra].push_back(phe_idx);

            // store the effect sizes
            pair_map[phe_idx][cpra] = std::make_pair(beta, se);

            ++n_pairs;
        }
    }
    notice("Finished reading %d trait-variant pairs across %zu traits and %zu variants", n_pairs, trait_ids.size(), var2traits.size());

    // construct PRS for each variant
    int32_t n_geno_samples = mpr.get_loaded_sample_count();
    int32_t n_phe = trait_ids.size();
    uint64_t n_geno_missing = 0;
    // initialize the PRS matrix to 0
    Eigen::MatrixXd prs_mat(n_geno_samples, n_phe);
    prs_mat.setZero();
    Eigen::MatrixXd var_mat(n_geno_samples, n_phe);
    var_mat.setZero();
    Eigen::MatrixXd count_mat(n_geno_samples, n_phe);
    count_mat.setZero();
    Eigen::VectorXd geno_vec(n_geno_samples);
    geno_vec.setZero();
    Eigen::Vector<bool, Eigen::Dynamic> mask_vec(n_geno_samples);
    mask_vec.setZero();
    Eigen::Vector<bool, Eigen::Dynamic> mask_any_present(n_phe);
    mask_any_present.setZero();
    
    // iterate each variant 
    // every variant is counted once: scored (n_pass), not found or unreadable (n_skip), or
    // monomorphic among the samples (n_mono)
    int32_t n_pass = 0, n_skip = 0, n_mono = 0;
    std::map<cpra_t, std::vector<int32_t> >::iterator var2traits_it;
    for(var2traits_it = var2traits.begin(); var2traits_it != var2traits.end(); ++var2traits_it) {
        const int32_t n_done = n_pass + n_skip + n_mono;
        if ( n_done > 0 && n_done % 100 == 0 ) {
            notice("Processed %d of %zu variants: %d scored, %d not found, %d monomorphic", n_done, var2traits.size(), n_pass, n_skip, n_mono);
        }
        cpra_t cpra = var2traits_it->first;
        const std::vector<int32_t>& phe_idxs = var2traits_it->second;
        // extract the genotypes
        // notice("int_buf.size() = %zu", int_buf.size());
        std::string cpra_s(cpra.to_string());
        if ( !mpr.read_pivar(cpra_s.c_str()) ) {
            notice("Skipping variant %s which is not found in the pvar file", cpra_s.c_str());
            ++n_skip;
            continue;
        }   
        //notice("Reading genotypes for variant %s", cpra_s.c_str());
        if ( !mpr.get_genos() ) {
            notice("Skipping variant %s, which failed to read genotypes", cpra_s.c_str());
            ++n_skip;
            continue;
        }
        const std::vector<int32_t>& int_buf = mpr.get_int_buf();
        const double* dbl_buf = mpr.get_dbl_buf();
        //notice("Finished reading genotypes for variant %s", cpra_s.c_str());
        int32_t an = 0;
        double ac = 0;
        double info = 0;
        if ( mpr.is_dosage_present() ) {
            if ( dbl_buf == NULL) {
                dbl_buf = mpr.get_dbl_buf();
            }
            double sumsq = 0;
            for(int32_t i =0; i < n_geno_samples; ++i) {
                //notice("Dosage[%d] = %.5g", i, dbl_buf[i]);
                double ds = 2.0 - dbl_buf[i];
                geno_vec(i) = ds;
                mask_vec(i) = true; // not missing
                ac += ds;
                an += 2;
                sumsq += (ds * ds);
            } 
            if ( ac == 0 || an == ac ) {
                // skip monomorphic variants
                ++n_mono;
                continue;
            }

            // E(Var(g)) = 2 * af * (1-af) = 2 * ac * ( an - ac ) / an / an
            // Var(g) = sumsq / (an / 2) - 4 * ac * ac / an / an
            // ratio = ( sumsq / 2 * an - 4 * ac * ac ) / ( 2 * ac * ( an - ac ) )
            info = ( sumsq / 2 * an - 4 * ac * ac ) / ( 2 * ac * ( an - ac ) );
        }
        else {
            int32_t gcs[3] = {0, 0, 0};
            for(int32_t i = 0; i < n_geno_samples; ++i) {
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
                ++n_mono;
                continue;
            }

            double mean = (double)ac / (double)an * 2.0;
            // Exp(Var(g)) = af * (1-af) * 2 = mean * (2 - mean) / 2;
            // Var(g) = EX^2 - EX^2 = (4 * n_2 + 1 * n_1)/n - mean^2
            info = ((4.0 * gcs[2] + gcs[1])/(an/2.0) - mean * mean) / (mean * (2.0 - mean) / 2.0);
            //notice("mean = %.5g, an = %d, ac = %d", mean, an, ac);
            for(int32_t i = 0; i < n_geno_samples; ++i) {
                switch(int_buf[i]) {
                case 0:
                    geno_vec(i) = 2.0 - mean; // homalt
                    mask_vec(i) = true; // not missing
                    break;
                case 1:
                    geno_vec(i) = 1.0 - mean; // het
                    mask_vec(i) = true; // not missing
                    break;
                case 2:
                    geno_vec(i) = 0.0 - mean; // homref
                    mask_vec(i) = true; // not missing
                    break;
                default:
                    geno_vec(i) = 0; // missing - mean imputation
                    mask_vec(i) = false; // missing
                    ++n_geno_missing;
                    break;
                }
            }
        }
        ++n_pass;
        // update each trait
        for(int32_t i=0; i < phe_idxs.size(); ++i) {
            int32_t phe_idx = phe_idxs[i];
            // add the effect size to the PRS
            std::pair<double, double> beta_se = pair_map[phe_idx][cpra];
            // update the PRS
            prs_mat.col(phe_idx) += beta_se.first * geno_vec;
            // update the variance
            var_mat.col(phe_idx).array() += (beta_se.second * geno_vec).array().square();
            // update the count
            count_mat.col(phe_idx) += mask_vec.cast<double>();
            mask_any_present(phe_idx) = true;
        }        
    }
    notice("Processed all %zu variants: %d scored, %d not found, %d monomorphic", var2traits.size(), n_pass, n_skip, n_mono);

    std::string out_prs = outf + ".prs.tsv.gz";
    htsFile* wf = hts_open(out_prs.c_str(), "wz");
    if ( wf == NULL ) {
        error("Cannot open output file %s for writing", outf.c_str());
    }
    // write the predicted PRS matrix in the format of Regenie phenotype files
    hprintf(wf, "FID\tIID");
    for(int32_t i=0; i < trait_ids.size(); ++i) {
        if ( mask_any_present(i) ) {
            hprintf(wf, "\t%s", trait_ids[i].c_str());
        }
    }
    hprintf(wf, "\n");
    // write down the PRS matrix
    for(int32_t i=0; i < n_geno_samples; ++i) {
        const plink_samp_t& samp = mpr.get_loaded_sample(i);
        hprintf(wf, "%s\t%s", samp.famID.empty() ? samp.indID.c_str() : samp.famID.c_str(), samp.indID.c_str());
        for(int32_t j=0; j < trait_ids.size(); ++j) {
            if ( count_mat(i, j) > 0 ) {
                hprintf(wf, "\t%.5g", prs_mat(i, j));
            }
            else if ( mask_any_present(j) ) {
                hprintf(wf, "\tNA");
            }
        }
        hprintf(wf, "\n");
    }
    hts_close(wf); // close the output file

    if ( icol_pair_se >= 0 ) {
        std::string out_se = outf + ".se.tsv.gz";
        htsFile* wf = hts_open(out_se.c_str(), "wz");
        if ( wf == NULL ) {
            error("Cannot open output file %s for writing", outf.c_str());
        }
        // write the predicted PRS matrix in the format of Regenie phenotype files
        hprintf(wf, "FID\tIID");
        for(int32_t i=0; i < trait_ids.size(); ++i) {
            if ( mask_any_present(i) ) {
                hprintf(wf, "\t%s", trait_ids[i].c_str());
            }
        }
        hprintf(wf, "\n");
        // write down the SE matrix
        for(int32_t i=0; i < n_geno_samples; ++i) {
            const plink_samp_t& samp = mpr.get_loaded_sample(i);
            hprintf(wf, "%s\t%s", samp.famID.empty() ? samp.indID.c_str() : samp.famID.c_str(), samp.indID.c_str());
            for(int32_t j=0; j < trait_ids.size(); ++j) {
                if ( count_mat(i, j) > 0 ) {
                    hprintf(wf, "\t%.5g", sqrt(var_mat(i, j)));
                }
                else if ( mask_any_present(j) ) {
                    hprintf(wf, "\tNA");
                }
            }
            hprintf(wf, "\n");
        }
        hts_close(wf); // close the output file
    }
    notice("Analysis finished");
    return 0;
}