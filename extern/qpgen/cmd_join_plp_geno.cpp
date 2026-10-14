#include "qgenlib/params.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_utils.h"
#include "qgenlib/phred_helper.h"
#include "qgenlib/hts_utils.h"
#include "qpgen.h"
#include <cmath>

int32_t join_plp_geno(int32_t argc, char **argv)
{
    std::string pfile;
    std::string bfile;
    std::string plpf;
    std::string varf;
    std::string sampID;
    std::string outf;
    bool hetOnly = false;
    paramList pl;

    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input pileup", NULL)
    LONG_STRING_PARAM("plp", &plpf, "Input pileup file (TSV format)")
    LONG_STRING_PARAM("var", &varf, "Input variant file associated with pileups")

    LONG_PARAM_GROUP("Input genotypes", NULL)
    LONG_STRING_PARAM("pfile", &pfile, "Input PLINK 2.0 file prefix (PGEN format)")
    LONG_STRING_PARAM("bfile", &bfile, "Input PLINK 1.9 file prefix (BED format)")
    LONG_STRING_PARAM("id", &sampID, "Sample ID that should be used for joining")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output file that joins the pileup and genotypes")

    LONG_PARAM_GROUP("Other options", NULL)
    LONG_PARAM("het-only", &hetOnly, "Output only heterozygous genotypes")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if (pfile.empty() + bfile.empty() != 1)
        error("Only one of --pfile or --bfile must be specified");

    if (sampID.empty())
        error("--sampleID must be specified");

    if (outf.empty())
        error("Output file name must be specified");

    // open the Plink file
    PlinkReader pr;
    if (!pfile.empty())
    {
        std::string pgenf = pfile + ".pgen";
        std::string psamf = pfile + ".psam";
        std::string pvarf = pfile + ".pvar";
        if (!pr.prep_pgen(pgenf.c_str(), pvarf.c_str(), psamf.c_str()))
            error("Failed to prepare the prefix %s", pfile.c_str());
    }
    else
    {
        std::string bedf = bfile + ".bed";
        std::string bimf = bfile + ".bim";
        std::string famf = bfile + ".fam";
        if (!pr.prep_bed(bedf.c_str(), bimf.c_str(), famf.c_str()))
            error("Failed to prepare the prefix %s", bfile.c_str());
    }

    std::vector<std::string> samp_ids;
    samp_ids.push_back(sampID);
    pr.set_filter_sample_id(samp_ids);

    // Read the pileup variants and identify the variant IDs
    notice("Loading the variant list from %s", varf.c_str());
    tsv_reader tr_var(varf.c_str());
    std::vector<std::string> cpras;
    tr_var.read_line(); // skip the header
    char buf[256];
    while (tr_var.read_line())
    {
        const char *s = strncmp(tr_var.str_field_at(1), "chr", 3) == 0 ? tr_var.str_field_at(1) + 3 : tr_var.str_field_at(1); // remove chr prefix
        snprintf(buf, 256, "%s:%d:%s:%s", s, tr_var.int_field_at(2), tr_var.str_field_at(3), tr_var.str_field_at(4));
        cpras.push_back(buf);
    }
    notice("Finished loading %zu variants", cpras.size());
    // notice("cpras[0] = %s", cpras[0].c_str());

    std::vector<int32_t> variant_idx;
    int32_t nmatches = pr.get_variant_idx_from_cpra(cpras, variant_idx);

    notice("Identified %d variants matching with pileup", nmatches);

    // read the pileup file
    notice("Loading the pileup file from %s", plpf.c_str());
    tsv_reader tr_plp(plpf.c_str());
    tr_plp.read_line(); // skip the header

    int32_t nsamps = pr.samp_idx.size();

    if ( nsamps == 0 )
        error("Sample ID %s could not be found in the input file", sampID.c_str());

    // print the summary file
    htsFile *wh = hts_open(outf.c_str(), outf.compare(outf.size()-3, 3, ".gz") == 0 ? "wz" : "w");
    if (wh == NULL)
        error("Failed to open %s", outf.c_str());
    hprintf(wh, "#Barcode\tVariant\tGenotype\tDepth\tnR\tnA\tnO\tAllele\tBQ\n");

    while (tr_plp.read_line())
    {
        int32_t vidx = tr_plp.int_field_at(1);
        if (variant_idx[vidx] < 0)
            continue; // skip if the variant is not in the list

        // get the genotypes
        pr.get_genos_at(variant_idx[vidx]);

        int32_t geno = pr.int_buf[0]; // AA:0, RA:1, RR:2
        if (geno > 2) {
            geno = -1;
        }

        int32_t depth = tr_plp.int_field_at(2);
        const char *als = tr_plp.str_field_at(3);
        const char *bqs = tr_plp.str_field_at(4);

        int32_t nA = 0, nR = 0, nO = 0;
        for(int32_t i=0; i < depth; ++i) {
            if (als[i] == '0') {
                ++nR;
            } else if (als[i] == '1') {
                ++nA;
            } else {
                ++nO;
            }
        }

        if ( hetOnly == false || geno == 1 ) {
            hprintf(wh, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n",
                        tr_plp.str_field_at(0), vidx, geno, depth,
                        nA, nR, nO, als, bqs);
        }

        if (vidx % 10000 == 0)
        {
            notice("Processing variant index %d at %s", vidx);
        }
    }

    hts_close(wh);

    notice("Analysis Finished");

    return 0;
}
