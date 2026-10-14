#include "qgenlib/params.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_utils.h"
#include "qgenlib/phred_helper.h"
#include "qgenlib/hts_utils.h"
#include "qpgen.h"
#include <cmath>

int32_t match_plp_geno(int32_t argc, char** argv) {
  std::string pfile;
  std::string bfile;
  std::string plpf;
  std::string varf;
  std::string samplist;
  std::string sampf;
  int32_t maxdepth = 100;
  int32_t topk = 10;
  bool skipfull = false;
  std::string outf;
  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input pileup", NULL)
    LONG_STRING_PARAM("plp", &plpf, "Input pileup file (TSV format)")
    LONG_STRING_PARAM("var", &varf, "Input variant file associated with pileups")

    LONG_PARAM_GROUP("Input genotypes", NULL)
    LONG_STRING_PARAM("pfile", &pfile, "Input PLINK 2.0 file prefix (PGEN format)")
    LONG_STRING_PARAM("bfile", &bfile, "Input PLINK 1.9 file prefix (BED format)")
    LONG_STRING_PARAM("sample-list", &samplist, "Comma-separated list of samples to be included (optional)")
    LONG_STRING_PARAM("sample-file", &sampf, "File containing the list of samples to be included (optional)")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output prefix")

    LONG_PARAM_GROUP("Other options", NULL)
    LONG_INT_PARAM("max-depth", &maxdepth, "Maximum depth per site to calculate the likelihood (default: 100)")
    LONG_INT_PARAM("top-k", &topk, "Number of top matching samples to report in the summary file")
    LONG_PARAM("skip-full", &skipfull, "Skip writing the full report containing all samples")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( pfile.empty() + bfile.empty() != 1 )
    error("Only one of --pfile or --bfile must be specified");

  if ( samplist.empty() + sampf.empty() == 0 )
    error("Only one of --sample-file or --sample-list must be specified");

  if ( outf.empty() ) 
    error("Output prefix (--out) must be specified");

  // open the Plink file
  PlinkReader pr;
  if ( !pfile.empty() ) {
    std::string pgenf = pfile + ".pgen";
    std::string psamf = pfile + ".psam";
    std::string pvarf = pfile + ".pvar";
    if ( !pr.prep_pgen(pgenf.c_str(), pvarf.c_str(), psamf.c_str()) ) 
      error("Failed to prepare the prefix %s", pfile.c_str());
  }
  else {
    std::string bedf = bfile + ".bed";
    std::string bimf = bfile + ".bim";
    std::string famf = bfile + ".fam";
    if ( !pr.prep_bed(bedf.c_str(), bimf.c_str(), famf.c_str()) )
      error("Failed to prepare the prefix %s", bfile.c_str());
  }

  // read the sample ids to subset
  if ( !samplist.empty() ) { // sample list was specified
    std::vector<std::string> samp_ids;
    split(samp_ids, ",", samplist);
    pr.set_filter_sample_id(samp_ids);
  }
  else if ( !sampf.empty() ) {
    tsv_reader tr_samp(sampf.c_str());
    std::vector<std::string> samp_ids;
    while( tr_samp.read_line() ) {
      samp_ids.push_back(tr_samp.str_field_at(0));
    }
    pr.set_filter_sample_id(samp_ids);
  }
  notice("Number of samples to be included = %zu", pr.samp_idx.size());

  // Read the pileup variants and identify the variant IDs
  notice("Loading the variant list from %s", varf.c_str());
  tsv_reader tr_var(varf.c_str());
  std::vector<std::string> cpras;
  tr_var.read_line(); // skip the header
  char buf[256];
  while( tr_var.read_line() ) {
    const char* s = strncmp(tr_var.str_field_at(1), "chr", 3) == 0 ? tr_var.str_field_at(1) + 3 : tr_var.str_field_at(1); // remove chr prefix
    snprintf(buf, 256, "%s:%d:%s:%s", s, tr_var.int_field_at(2), tr_var.str_field_at(3), tr_var.str_field_at(4));
    cpras.push_back(buf);
  }
  notice("Finished loading %zu variants", cpras.size());
  //notice("cpras[0] = %s", cpras[0].c_str());

  std::vector<int32_t> variant_idx;
  int32_t nmatches = pr.get_variant_idx_from_cpra(cpras, variant_idx);

  notice("Identified %d variants matching with pileup", nmatches);

  // read the pileup file
  notice("Loading the pileup file from %s", plpf.c_str());
  tsv_reader tr_plp(plpf.c_str());
  tr_plp.read_line(); // skip the header

  int32_t nsamps = pr.samp_idx.size();

  std::vector<double> llks(nsamps, 0);
  double logten = log(10);
  double log3 = log(3)/logten;
  uint64_t nreads = 0;

  while( tr_plp.read_line() ) {
    int32_t vidx = tr_plp.int_field_at(1);
    if ( variant_idx[vidx] < 0 ) continue; // skip if the variant is not in the list

    // get the genotypes
    pr.get_genos_at(variant_idx[vidx]);

    // calculate per-site likelihood for each possible genotype
    double lkRR = 0, lkRA = 0, lkAA = 0;
    int32_t depth = tr_plp.int_field_at(2);
    if ( depth > maxdepth ) depth = maxdepth;
    const char* als = tr_plp.str_field_at(3);
    const char* bqs = tr_plp.str_field_at(4);
    nreads += depth;
    for(int32_t i=0; i < depth; ++i) {
      int32_t phred = bqs[i] - 33;
      if ( als[i] == '0' ) {
        lkRR += phredConv.phred2LogMat[ phred ];
        lkRA += phredConv.phred2HalfLogMat3[ phred ];
        lkAA += ( -0.1 * phred - log3 );
      }
      else if ( als[i] == '1' ) {
        lkAA += phredConv.phred2LogMat[ phred ];
        lkRA += phredConv.phred2HalfLogMat3[ phred ];
        lkRR += ( -0.1 * phred - log3 );
      }
      else {
        lkRR += ( -0.1 * phred - log3 );
        lkRA += ( -0.1 * phred - log3 );
        lkAA += ( -0.1 * phred - log3 );
      }
    }

    // calculate allele frequencies
    int32_t an = 0, ac = 0;
    for(int32_t i=0; i < nsamps; ++i) {
      an += 2;
      ac += pr.int_buf[i];
    }
    double af = an == 0 ? 0.5 : (double)ac / (double)an; 
    double gfs[3] = { (1-af)*(1-af), 2*af*(1-af), af*af }; // AA, RA, RR
    double lkmax = lkAA;
    if ( lkRA > lkmax ) lkmax = lkRA;
    if ( lkRR > lkmax ) lkmax = lkRR;
    double lkmiss = lkmax + log( pow(10.0, lkmax - lkAA) * gfs[0] + pow(10.0, lkmax - lkRA) * gfs[1] + pow(10.0, lkmax - lkRR) * gfs[2] )/logten;

    // add log-likelihoods per sample
    for(int32_t i=0; i < nsamps; ++i) {
      if ( pr.int_buf[i] == 0 ) {
        llks[i] += lkAA;
      }
      else if ( pr.int_buf[i] == 1 ) {
        llks[i] += lkRA;
      }
      else if ( pr.int_buf[i] == 2 ) {
        llks[i] += lkRR;
      }
      else {  // if missing, assume genotype frequency for imputing missing genotypes
        llks[i] += lkmiss;
      }
    }

    // intermediate checking and print output
    if ( vidx % 10000 == 0 ) {
      // identify the best and the second best samples 
      int32_t best_idx = -1, second_idx = -1;
      double best_lk = -1e100, second_lk = -1e100;
      for(int32_t i=0; i < nsamps; ++i) {
        if ( llks[i] > best_lk ) {
          second_lk = best_lk;
          second_idx = best_idx;
          best_lk = llks[i];
          best_idx = i;
        }
        else if ( llks[i] > second_lk ) {
          second_lk = llks[i];
          second_idx = i;
        }
      }
      notice("Processing variant index %d at %s. Best : %s (%.5lf), Next : %s (%.5lf)", 
                vidx, cpras[vidx].c_str(), pr.samps[pr.samp_idx[best_idx]-1].indID.c_str(), best_lk,
                pr.samps[pr.samp_idx[second_idx]-1].indID.c_str(), second_lk);
    }
  }

  // identify top k matching samples
  std::vector<int32_t> best_idx(topk, -1);
  std::vector<double> best_llks(topk, -1e100);
  for(int32_t i=0; i < nsamps; ++i) {
    for(int32_t j=0; j < topk; ++j) {
      if ( llks[i] > best_llks[j] ) {
        for(int32_t k=topk-1; k > j; --k) {
          best_llks[k] = best_llks[k-1];
          best_idx[k] = best_idx[k-1];
        }
        best_llks[j] = llks[i];
        best_idx[j] = i;
        break;
      }
    }
  }

  // print the summary file
  htsFile* fp = hts_open((outf + ".summary").c_str(), "w");
  if ( fp == NULL ) 
    error("Failed to open %s", (outf + ".summary").c_str());
  
  hprintf(fp, "plp_file:\t%s\n", plpf.c_str());
  hprintf(fp, "var_file:\t%s\n", varf.c_str());
  hprintf(fp, "n_samples:\t%zu\n", pr.samp_idx.size());
  hprintf(fp, "n_variants:\t%d\n", nmatches);
  hprintf(fp, "n_reads:\t%llu\n", nreads);
  hprintf(fp, "max_depth:\t%d\n", maxdepth);
  for(int32_t i=0; i < topk; ++i) {
    if ( best_idx[i] >= 0 ) {
      hprintf(fp, "top%d:\t%s\t%.5lf\n", i+1, 
                pr.samps[pr.samp_idx[best_idx[i]]-1].indID.c_str(), 
                best_llks[i]);
    }
    else {
      hprintf(fp, "top%d:\tNA\tNA\n", i+1);
    }
  }
  hts_close(fp);

  fp = hts_open((outf + ".best").c_str(), "w");
  if ( fp == NULL ) 
    error("Failed to open %s", (outf + ".best").c_str());

  if ( best_idx[0] >= 0 ) {
    hprintf(fp, "%s\t%.5lf\n", 
                pr.samps[pr.samp_idx[best_idx[0]]-1].indID.c_str(), 
                best_llks[0]);
  }
  else {
    hprintf(fp, "NA\tNA\n");
  }
  hts_close(fp);

  if ( !skipfull ) {
    fp = hts_open((outf + ".full.tsv.gz").c_str(), "wz");
    if ( fp == NULL ) 
      error("Failed to open %s", (outf + ".best").c_str());

    hprintf(fp, "#SAMPLE_ID\tLOG10LLK\n");
    for(int32_t i=0; i < nsamps; ++i) {
      hprintf(fp, "%s\t%.5lf\n", 
                  pr.samps[pr.samp_idx[i]-1].indID.c_str(), 
                  llks[i]);
    }
    hts_close(fp);
  }

  return 0;
}
