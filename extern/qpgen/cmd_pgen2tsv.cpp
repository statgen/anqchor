#include "qgenlib/params.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_utils.h"
#include "qgenlib/hts_utils.h"
#include "qpgen.h"
#include "qpgen_utils.h"
#include <cmath>
#include <vector>
#include <string>

// pgen2tsv: extract PLINK 2.0 (.pgen) genotypes into TSV, using the modern
// PgenIdxReader. Unlike geno2tsv, this command:
//  - is PGEN-only (no .bed support)
//  - emits genotypes in ALT-dosage coding: 0 = REF/REF, 1 = het, 2 = ALT/ALT, NA = missing
//  - supports selecting variants either by explicit index (--var) or by genomic region (--region)
//  - supports allele-frequency / allele-count filters (--min-af/--max-af/--min-ac/--max-ac)
//  - opportunistically uses the .pgen's sparse (difflist) representation to read
//    and (in --sparse output mode) emit only the non-REF/REF carriers.
//
// Interpreting the sparse buffer (when reader.is_sparse() == true):
//   The reader returns the genotypes in REF-count coding (hom-ref = 2). The bulk
//   of the samples share a single "common" genotype, get_sparse_common_geno(),
//   and only the listed exceptions differ:
//     - get_sparse_sample_idxs()[k] : 0-based index (within the loaded sample
//                                     subset) of the k-th non-common sample
//     - get_sparse_genos()[k]       : that sample's genotype (REF-count, -9=missing)
//   Every sample NOT listed has genotype get_sparse_common_geno().
//   Here we convert REF-count -> ALT-count for output via alt = 2 - ref.

// Materialize the full per-sample ALT-coded genotype vector (0/1/2, -9=missing)
// from whichever representation the reader produced.
static void materialize_alt(PgenIdxReader& r, int nsamps, std::vector<int>& alt) {
    alt.assign(nsamps, 0);
    if (r.is_sparse()) {
        int32_t c_ref = r.get_sparse_common_geno();
        int c_alt = (c_ref < 0) ? -9 : (2 - c_ref);
        std::fill(alt.begin(), alt.end(), c_alt);
        const std::vector<int32_t>& ids = r.get_sparse_sample_idxs();
        const std::vector<int32_t>& gs  = r.get_sparse_genos();
        for (size_t k = 0; k < ids.size(); ++k) {
            int g = gs[k];
            alt[ids[k]] = (g < 0) ? -9 : (2 - g);
        }
    } else {
        const std::vector<int>& ib = r.get_int_buf(); // REF-count dense
        for (int i = 0; i < nsamps; ++i) alt[i] = (ib[i] < 0) ? -9 : (2 - ib[i]);
    }
}

// Compute ALT allele count and non-missing allele number directly from the
// reader's representation (O(#carriers) when sparse).
static void compute_ac_an(PgenIdxReader& r, int nsamps, long& ac, long& an) {
    ac = 0; an = 0;
    if (r.is_sparse()) {
        int32_t c_ref = r.get_sparse_common_geno();
        const std::vector<int32_t>& ids = r.get_sparse_sample_idxs();
        const std::vector<int32_t>& gs  = r.get_sparse_genos();
        int n_common = nsamps - (int)ids.size();
        if (c_ref >= 0) { ac += (long)(2 - c_ref) * n_common; an += 2L * n_common; }
        for (size_t k = 0; k < gs.size(); ++k) {
            if (gs[k] >= 0) { ac += (2 - gs[k]); an += 2; }
        }
    } else {
        const std::vector<int>& ib = r.get_int_buf();
        for (int i = 0; i < nsamps; ++i) {
            if (ib[i] >= 0) { ac += (2 - ib[i]); an += 2; }
        }
    }
}

// Write the AC/AN columns followed by the genotype columns for the current variant.
// info_prefix is the already-tab-joined leading columns (no trailing tab, may be empty).
static void write_variant(htsFile* wh, PgenIdxReader& r, int nsamps, bool sparse_out,
                          long ac, long an, const std::string& info_prefix,
                          std::vector<int>& altbuf) {
    if (!info_prefix.empty()) hprintf(wh, "%s\t", info_prefix.c_str());
    hprintf(wh, "%ld\t%ld", ac, an);

    if (sparse_out && r.is_sparse() && r.get_sparse_common_geno() == 2) {
        // Common genotype is REF/REF (ALT count 0); emit only the carriers.
        // pgenlib returns difflist sample ids in ascending order.
        const std::vector<int32_t>& ids = r.get_sparse_sample_idxs();
        const std::vector<int32_t>& gs  = r.get_sparse_genos();
        for (size_t k = 0; k < ids.size(); ++k) {
            int g = gs[k];
            if (g < 0)            hprintf(wh, "\t%d:NA", ids[k]);
            else if ((2 - g) != 0) hprintf(wh, "\t%d:%d", ids[k], 2 - g);
        }
    } else if (sparse_out) {
        materialize_alt(r, nsamps, altbuf);
        for (int i = 0; i < nsamps; ++i) {
            if (altbuf[i] < 0)      hprintf(wh, "\t%d:NA", i);
            else if (altbuf[i] != 0) hprintf(wh, "\t%d:%d", i, altbuf[i]);
        }
    } else {
        materialize_alt(r, nsamps, altbuf);
        for (int i = 0; i < nsamps; ++i) {
            if (altbuf[i] < 0) hprintf(wh, "\tNA");
            else               hprintf(wh, "\t%d", altbuf[i]);
        }
    }
    hprintf(wh, "\n");
}

int32_t cmd_pgen2tsv(int32_t argc, char **argv)
{
    std::string pfile;   // PLINK2 prefix (derives .pgen/.psam/.pvar.idx.gz)
    std::string pgenf;
    std::string psamf;
    std::string pivarf;
    std::string varf;    // explicit variant-index file
    std::string region;  // CHROM:BEG-END
    std::string samplist;
    std::string sampf;
    std::string outf;
    int32_t ibegin = 0;
    int32_t icol_pivar_idx = 9; // 1-based column of the variant index in the pvar file (this GTEx example uses 8)
    double min_af = 0.0, max_af = 1.0;
    double min_ac = 0.0, max_ac = 1e18;
    bool sparse = false;
    bool space_delim = false; // space- instead of tab-delimited input tables

    paramList pl;

    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input genotypes", NULL)
    LONG_STRING_PARAM("pfile", &pfile, "Input PLINK 2.0 file prefix (expects .pgen/.psam/.pvar.idx.gz)")
    LONG_STRING_PARAM("pgen", &pgenf, "Input PLINK 2.0 genotype file (.pgen); overrides --pfile")
    LONG_STRING_PARAM("psam", &psamf, "Input PLINK 2.0 sample file (.psam); overrides --pfile")
    LONG_STRING_PARAM("pivar", &pivarf, "Input tabixed/indexed pvar file (required for --region)")
    LONG_STRING_PARAM("sample-list", &samplist, "Comma-separated list of samples to include (optional)")
    LONG_STRING_PARAM("sample-file", &sampf, "File containing the list of samples to include (optional)")

    LONG_PARAM_GROUP("Variant selection (use one of --var or --region)", NULL)
    LONG_STRING_PARAM("var", &varf, "Input variant file (TSV) with [variant index] [extra columns]")
    LONG_INT_PARAM("idx-begin", &ibegin, "Starting index for --var (default: 0). Use 1 for 1-based indices")
    LONG_STRING_PARAM("region", &region, "Genomic region CHROM:BEG-END to stream (requires --pivar)")
    LONG_INT_PARAM("icol-pivar-idx", &icol_pivar_idx, "1-based column index for the variant ID in the pvar file (default: 9)")

    LONG_PARAM_GROUP("Variant filters", NULL)
    LONG_DOUBLE_PARAM("min-af", &min_af, "Minimum ALT allele frequency (default: 0.0)")
    LONG_DOUBLE_PARAM("max-af", &max_af, "Maximum ALT allele frequency (default: 1.0)")
    LONG_DOUBLE_PARAM("min-ac", &min_ac, "Minimum ALT allele count (default: 0)")
    LONG_DOUBLE_PARAM("max-ac", &max_ac, "Maximum ALT allele count (default: 1e18)")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output file (use .gz suffix for bgzipped output)")
    LONG_PARAM("sparse", &sparse, "Emit sparse [IDX]:[GENO] for non-REF/REF samples instead of full genotypes")
    LONG_PARAM("space-delim", &space_delim, QPGEN_SPACE_DELIM_HELP)
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();
    const char table_delim = table_delimiter(space_delim);

    // resolve inputs from --pfile prefix when explicit paths are not given
    if (!pfile.empty()) {
        if (pgenf.empty())  pgenf  = pfile + ".pgen";
        if (psamf.empty())  psamf  = pfile + ".psam";
        if (pivarf.empty()) pivarf = pfile + ".pvar.idx.gz";
    }
    if (pgenf.empty() || psamf.empty())
        error("Genotype (--pgen/--pfile) and sample (--psam/--pfile) inputs are required");
    if (varf.empty() + region.empty() != 1)
        error("Exactly one of --var or --region must be specified");
    if (!region.empty() && pivarf.empty())
        error("--region requires an indexed pvar file (--pivar/--pfile)");
    if (samplist.empty() + sampf.empty() == 0)
        error("Only one of --sample-file or --sample-list must be specified");
    if (outf.empty())
        error("Output file (--out) must be specified");

    // open the reader and load all samples
    PgenIdxReader reader;
    reader.set_icol_pivar_idx(icol_pivar_idx - 1); // convert to 0-based
    if (!reader.prep_pgen(pgenf.c_str(), pivarf.empty() ? "" : pivarf.c_str(), psamf.c_str()))
        error("Failed to prepare pgen/psam files %s / %s", pgenf.c_str(), psamf.c_str());

    // subset samples if requested
    if (!samplist.empty()) {
        std::vector<std::string> samp_ids;
        split(samp_ids, ",", samplist);
        reader.subset_sample_ids(samp_ids);
    }
    else if (!sampf.empty()) {
        tsv_reader tr_samp(sampf.c_str());
        tr_samp.delimiter = table_delim;
        std::vector<std::string> samp_ids;
        while (tr_samp.read_line()) samp_ids.push_back(tr_samp.str_field_at(0));
        reader.subset_sample_ids(samp_ids);
    }
    int32_t nsamps = reader.get_loaded_sample_count();
    notice("Number of samples to be included = %d", nsamps);

    htsFile* wh = hts_open(outf.c_str(), outf.size() >= 3 && outf.compare(outf.size() - 3, 3, ".gz") == 0 ? "wz" : "w");
    if (wh == NULL) error("Cannot open output file %s", outf.c_str());

    std::vector<int> altbuf;
    int64_t n_in = 0, n_out = 0;

    if (!region.empty()) {
        // ---- stream by genomic region ----
        cbe_t cbe(region.c_str());
        // header
        hprintf(wh, "#CHROM\tPOS\tID\tREF\tALT\tAC\tAN");
        for (int32_t i = 0; i < nsamps; ++i)
            hprintf(wh, "\t%s", reader.get_loaded_sample(i).indID.c_str());
        hprintf(wh, "\n");

        bool have = reader.read_pos(cbe.chrom.c_str(), cbe.beg1);
        while (have) {
            const plink_var_t& v = reader.get_current_variant();
            if (v.pos > cbe.end0) break;
            ++n_in;
            reader.get_genos_sparse();

            long ac, an;
            compute_ac_an(reader, nsamps, ac, an);
            double af = (an > 0) ? (double)ac / (double)an : 0.0;
            bool keep = (an > 0) && (ac >= min_ac) && (ac <= max_ac) && (af >= min_af) && (af <= max_af);
            if (keep) {
                std::string prefix = v.to_string(); // CHROM:POS:REF:ALT not desired; build columns
                // build CHROM\tPOS\tID\tREF\tALT
                std::string alts;
                for (size_t a = 0; a < v.alts.size(); ++a) { if (a) alts += ","; alts += v.alts[a]; }
                char posbuf[32]; snprintf(posbuf, sizeof(posbuf), "%d", v.pos);
                std::string info = v.schrom + "\t" + posbuf + "\t" + v.vid + "\t" + v.ref + "\t" + alts;
                write_variant(wh, reader, nsamps, sparse, ac, an, info, altbuf);
                ++n_out;
            }
            have = reader.read_pivar(NULL); // advance to the next variant
        }
    }
    else {
        // ---- explicit variant indices from --var ----
        tsv_reader tr_var(varf.c_str());
        tr_var.delimiter = table_delim;
        bool header_written = false;
        uint32_t nlines = 0;
        while (tr_var.read_line()) {
            const char* s = tr_var.str_field_at(0);
            if (s[0] == '#') { // header line(s)
                if (s[1] == '#') continue; // meta line
                // echo the index-file columns (skip the leading index column), then AC/AN/samples
                std::string hdr;
                for (int32_t i = 1; i < tr_var.nfields; ++i) { if (i > 1) hdr += "\t"; hdr += tr_var.str_field_at(i); }
                if (!hdr.empty()) hprintf(wh, "%s\t", hdr.c_str());
                hprintf(wh, "AC\tAN");
                for (int32_t i = 0; i < nsamps; ++i)
                    hprintf(wh, "\t%s", reader.get_loaded_sample(i).indID.c_str());
                hprintf(wh, "\n");
                header_written = true;
                continue;
            }
            if (!header_written && nlines == 0) {
                // headerless: synthesize V1,V2,... for the echoed columns
                std::string hdr;
                for (int32_t i = 1; i < tr_var.nfields; ++i) { if (i > 1) hdr += "\t"; hdr += "V" + std::to_string(i); }
                if (!hdr.empty()) hprintf(wh, "%s\t", hdr.c_str());
                hprintf(wh, "AC\tAN");
                for (int32_t i = 0; i < nsamps; ++i)
                    hprintf(wh, "\t%s", reader.get_loaded_sample(i).indID.c_str());
                hprintf(wh, "\n");
                header_written = true;
            }
            // data line
            int32_t idx = tr_var.int_field_at(0) - ibegin;
            ++n_in;
            if (!reader.get_genos_sparse(idx))
                error("Failed to get genotypes for variant index %d", idx);

            long ac, an;
            compute_ac_an(reader, nsamps, ac, an);
            double af = (an > 0) ? (double)ac / (double)an : 0.0;
            bool keep = (an > 0) && (ac >= min_ac) && (ac <= max_ac) && (af >= min_af) && (af <= max_af);
            if (keep) {
                std::string info;
                for (int32_t i = 1; i < tr_var.nfields; ++i) { if (i > 1) info += "\t"; info += tr_var.str_field_at(i); }
                write_variant(wh, reader, nsamps, sparse, ac, an, info, altbuf);
                ++n_out;
            }
            ++nlines;
        }
        tr_var.close();
    }

    hts_close(wh);
    notice("Processed %lld variants, wrote %lld after filtering", (long long)n_in, (long long)n_out);
    notice("Analysis finished");
    return 0;
}
