#include "qgenlib/params.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_utils.h"
#include "qgenlib/phred_helper.h"
#include "qgenlib/hts_utils.h"
#include "qpgen.h"
#include <cmath>

int32_t cmd_geno2tsv(int32_t argc, char **argv)
{
    std::string pfile;
    std::string bfile;
    std::string varf;
    std::string samplist;
    std::string sampf;
    std::string outf;
    int32_t ibegin = 0;
    //bool dosage = false;
    bool sparse = false;
    paramList pl;

    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input variant info", NULL)
    LONG_STRING_PARAM("var", &varf, "Input variant file (TSV) containing [variant index] [extra columns]")
    LONG_INT_PARAM("idx-begin", &ibegin, "Starting index for the variant (default: 0). Use 1 for 1-based index")

    LONG_PARAM_GROUP("Input genotypes", NULL)
    LONG_STRING_PARAM("pfile", &pfile, "Input PLINK 2.0 file prefix (PGEN format)")
    LONG_STRING_PARAM("bfile", &bfile, "Input PLINK 1.9 file prefix (BED format)")
    LONG_STRING_PARAM("sample-list", &samplist, "Comma-separated list of samples to be included (optional)")
    LONG_STRING_PARAM("sample-file", &sampf, "File containing the list of samples to be included (optional)")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output prefix")
    LONG_PARAM("sparse", &sparse, "Output sparse representation [IDX]:[GENO] instead of full genotypes")
    //LONG_PARAM("dosage", &dosage, "Output dosage instead of genotypes (default: genotypes). Only available for PGEN format")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if (pfile.empty() + bfile.empty() != 1)
        error("Only one of --pfile or --bfile must be specified");

    if (samplist.empty() + sampf.empty() == 0)
        error("Only one of --sample-file or --sample-list must be specified");

    if (outf.empty())
        error("Output prefix (--out) must be specified");

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

    // read the sample ids to subset
    if (!samplist.empty())
    { // sample list was specified
        std::vector<std::string> samp_ids;
        split(samp_ids, ",", samplist);
        pr.set_filter_sample_id(samp_ids);
    }
    else if (!sampf.empty())
    {
        tsv_reader tr_samp(sampf.c_str());
        std::vector<std::string> samp_ids;
        while (tr_samp.read_line())
        {
            samp_ids.push_back(tr_samp.str_field_at(0));
        }
        pr.set_filter_sample_id(samp_ids);
    }
    notice("Number of samples to be included = %zu", pr.samp_idx.size());

    int32_t nsamps = (int32_t)pr.samp_idx.size();


    // Read the pileup variants and identify the variant IDs
    notice("Loading the variant list from %s", varf.c_str());
    tsv_reader tr_var(varf.c_str());
    htsFile *wh = hts_open(outf.c_str(), outf.compare(outf.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");
    bool has_header = false;
    uint32_t nlines = 0;
    while (tr_var.read_line())
    {
        const char *s = tr_var.str_field_at(0);
        if (s[0] == '#')
        { // header line
            if (s[1] == '#')
            { // meta line - skip}
                continue;
            }
            else
            { // output the header, except for the index column
                for (int32_t i = 1; i < tr_var.nfields; ++i)
                {
                    if (i > 1)
                        hprintf(wh, "\t");
                    hprintf(wh, "%s", tr_var.str_field_at(i));
                }
                if (tr_var.nfields > 1)
                    hprintf(wh, "\t");
                hprintf(wh, "AC\tAN");
                for (int32_t i = 0; i < nsamps; ++i)
                {
                    hprintf(wh, "\t%s", pr.samps[pr.samp_idx[i]-1].indID.c_str());
                }
                hprintf(wh, "\n");
                has_header = true;
            }
        }
        else if (nlines == 0)
        { // first data line without headers
            for (int32_t i = 1; i < tr_var.nfields; ++i)
            {
                if (i > 1)
                    hprintf(wh, "\t");
                hprintf(wh, "V%d", i);
            }
            if (tr_var.nfields > 1)
                hprintf(wh, "\t");
            hprintf(wh, "AC\tAN");
            for (int32_t i = 0; i < nsamps; ++i)
            {
                hprintf(wh, "\t%s", pr.samps[pr.samp_idx[i]-1].indID.c_str());
            }
            hprintf(wh, "\n");
        }
        else
        { // data line
            // get the genotypes
            int32_t idx = tr_var.int_field_at(0);
            if (!pr.get_genos_at(idx - ibegin))
                error("Failed to get genotypes for variant %d", idx);

            int32_t ac = 0;
            int32_t an = 0;
            for (int32_t i = 0; i < nsamps; ++i)
            {
                switch (pr.int_buf[i])
                {
                case 0:
                    an += 2;
                    ac += 2;
                    break;
                case 1:
                    an += 2;
                    ++ac;
                    break;
                case 2:
                    an += 2;
                    break;
                }
            }

            // print the variant info
            for (int32_t i = 1; i < tr_var.nfields; ++i)
            {
                if (i > 1)
                    hprintf(wh, "\t");
                hprintf(wh, "%s", tr_var.str_field_at(i));
            }
            if (tr_var.nfields > 1)
                hprintf(wh, "\t");
            hprintf(wh, "%d\t%d", ac, an);
            if ( sparse ) {
                for (int32_t i = 0; i < nsamps; ++i)
                {
                    switch (pr.int_buf[i])
                    {
                    case 0:
                        hprintf(wh, "\t%d:2", i);
                        break;
                    case 1:
                        hprintf(wh, "\t%d:1", i);
                        break;
                    case 2:
                        //hprintf(wh, "0");
                        break;
                    default: // missing
                        hprintf(wh, "\t%d:NA", i);
                        break;
                    }
                }
                hprintf(wh, "\n");
            }
            else {
                for (int32_t i = 0; i < nsamps; ++i)
                {
                    switch (pr.int_buf[i])
                    {
                    case 0:
                        hprintf(wh, "\t2");
                        break;
                    case 1:
                        hprintf(wh, "\t1");
                        break;
                    case 2:
                        hprintf(wh, "\t0");
                        break;
                    default: // missing
                        hprintf(wh, "\tNA");
                        break;
                    }
                }
                hprintf(wh, "\n");
            }
        }
        ++nlines;
    }
    notice("Finished reading %llu lines", nlines);
    tr_var.close();
    hts_close(wh);

    notice("Analysis finished");
    return 0;
}
