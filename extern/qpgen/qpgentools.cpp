#include "qpgentools.h"
#include "qgenlib/commands.h"
#include "qgenlib/qgen_utils.h"
#include "qgenlib/qgen_error.h"
#include "qpgen_version.h" // generated at build time (cmake/qpgen_version.cmake)

int32_t test_qpgen(int32_t argc, char** argv);
int32_t match_plp_geno(int32_t argc, char** argv);
int32_t join_plp_geno(int32_t argc, char** argv);
int32_t cmd_geno2tsv(int32_t argc, char** argv);
int32_t cmd_pgen2tsv(int32_t argc, char** argv);
int32_t cmd_pair_assoc(int32_t argc, char** argv);
int32_t cmd_pair_assoc_v0(int32_t argc, char** argv);
int32_t cmd_rect_assoc(int32_t argc, char** argv);
int32_t cmd_pair_prs(int32_t argc, char** argv);
int32_t cmd_match_prs_pheno(int32_t argc, char** argv);
int32_t cmd_region_assoc(int32_t argc, char** argv);
int32_t main(int32_t argc, char** argv) {
  commandList cl;

  BEGIN_LONG_COMMANDS(longCommandlines)
    LONG_COMMAND_GROUP("Utilities for xqtls", NULL)
    LONG_COMMAND("pair-assoc", &cmd_pair_assoc, "Perform pairwise association analysis for specific pairs of phenotype variant pairs")
    LONG_COMMAND("pair-assoc-v0", &cmd_pair_assoc_v0, "Perform pairwise association analysis for specific pairs of phenotype variant pairs (old version)")
    LONG_COMMAND("pair-prs", &cmd_pair_prs, "Perform pairwise PRS generation based on summary statistics")
    LONG_COMMAND("rect-assoc", &cmd_rect_assoc, "Perform rectangular association analysis for multiple phenotypes and variants")
    LONG_COMMAND("region-assoc", &cmd_region_assoc, "Perform association analysis for a specific region and multiple phenotypes")
    LONG_COMMAND("match-prs-pheno", &cmd_match_prs_pheno, "Match PRS and phenotype matrices based on overlapping samples and compute weights for each phenotype")

    LONG_COMMAND_GROUP("Utilities for sequence data", NULL)
    LONG_COMMAND("match-plp-geno", &match_plp_geno, "Check the concordance between pileup and genotypes")
    LONG_COMMAND("join-plp-geno",  &join_plp_geno, "Join the pileup and individual genotypes")

    LONG_COMMAND_GROUP("Utilities for genotype extraction", NULL)
    LONG_COMMAND("geno2tsv", &cmd_geno2tsv, "Extract genotypes into TSV format")
    LONG_COMMAND("pgen2tsv", &cmd_pgen2tsv, "Extract PGEN genotypes (0/1/2 ALT coding) into TSV with AF/AC filters and sparse support")

    LONG_COMMAND_GROUP("Other Utilities", NULL)
    LONG_COMMAND("test-qpgen", &test_qpgen, "Test software to check the qpgen library")
  END_LONG_COMMANDS();

  cl.Add(new longCommands("Available Commands", longCommandlines));
  
  if ( argc >= 2 && strcmp(argv[1], "--version") == 0 ) {
    printf("qpgentools %s (git %s)\n", QPGEN_VERSION, QPGEN_GIT);
    return 0;
  }

  if ( argc < 2 ) {
    printf("[qpgentools] -- Another command line toolkit for PLINK files\n");
    printf("Version %s (git %s)\n\n", QPGEN_VERSION, QPGEN_GIT);
    fprintf(stderr, " Copyright (c) 2023-2025 by Hyun Min Kang\n");
    fprintf(stderr, " Licensed under the Apache License v2.0 and LGPLv3 (see LICENSE for details)\n");    
    fprintf(stderr, "To run a specific command      : %s [command] [options]\n",argv[0]);
    fprintf(stderr, "For detailed instructions, run : %s --help\n",argv[0]);
    fprintf(stderr, "To print the version, run      : %s --version\n",argv[0]);
    cl.Status();
    return 1;
  }
  else {
    if ( strcmp(argv[1],"--help") == 0 ) {
      cl.HelpMessage();
    }
    else {
      // the version goes first in every command's log, for provenance
      notice("qpgentools %s (git %s)", QPGEN_VERSION, QPGEN_GIT);
      return cl.Read(argc, argv);
    }
  }
  return 0;
}
