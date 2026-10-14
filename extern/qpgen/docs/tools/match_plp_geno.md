# qpgentools match-plp-geno

## Summary

`qpgentools match-plp-geno` matches a sample with low-pass sequencing data (represented as a pileup file) to a set of reference genotypes. It calculates the likelihood of the pileup data given each sample's genotype in the reference panel and identifies the best matching sample. This is useful for sample authentication or identifying sample swaps.

A typical running example command is given below:

```bash
qpgentools match-plp-geno --plp [pileup_file] --var [variant_file] --pfile [reference_pgen] --out [out_prefix]
```

## Input File Formats

### Variant File (`--var`)

The variant file is a TSV file that defines the variants included in the pileup. It must contain a header line (which is skipped).

Columns expected:
1.  **Variant ID**: Variant ID 
2.  **CHROM**: Chromosome (e.g., `chr1`, `1`)
3.  **POS**: Position (1-based)
4.  **REF**: Reference allele
5.  **ALT**: Alternative allele

### Pileup File (`--plp`)

The pileup file is a TSV file containing the sequencing summary at each variant site. It must contain a header line (which is skipped).

Columns expected (0-based indices used in parsing):
0.  **Barcode**: Cell barcode of sample ID
1.  **Variant Index**: 0-based index of the variant in the `--var` file.
2.  **Depth**: Read depth at this site.
3.  **Alleles**: String representing alleles observed in reads (`0` for REF, `1` for ALT).
4.  **Base Qualities**: String representing base qualities for each read (Phred+33 ASCII).

**Example Pileup Row:**
```
barcode1    0    5    00010    IIIII
```
This indicates:
-   At Variant Index 0 (first variant in `.var` file),
-   Depth is 5,
-   Reads observed: Ref, Ref, Ref, Alt, Ref (`00010`),
-   Qualities are all 'I' (Phred score 40).

## Required options

*   **Input Pileup**:
    *   `--plp`: Input pileup file (TSV format).
    *   `--var`: Input variant file associated with pileups.
*   **Input Genotypes** (One of the following):
    *   `--pfile`: Input PLINK 2.0 file prefix (PGEN format).
    *   `--bfile`: Input PLINK 1.9 file prefix (BED format).
*   **Sample Filtering** (One of the following is required):
    *   `--sample-list`: Comma-separated list of sample IDs to include.
    *   `--sample-file`: File with one sample ID per line.
*   `--out`: Output prefix.

## Additional Options

*   `--max-depth` : Maximum depth per site to use for likelihood calculation (default: 100).
*   `--top-k` : Number of top matching samples to report (default: 10).
*   `--skip-full` : Do not write the full report of all samples (default: false).

## Expected Output

*   `[out_prefix].summary`: A summary file containing run parameters and the top-k matching samples with their log-likelihoods.
*   `[out_prefix].best`: A simple file containing `[SAMPLE_ID] [LOG_LIKELIHOOD]` of the single best match.
*   `[out_prefix].full.tsv.gz`: (Unless skipped) A gzip-compressed TSV containing `[SAMPLE_ID] [LOG10LLK]` for all samples tested.

## Full Usage

The full usage of `qpgentools match-plp-geno` can be viewed with the `--help` option:

```
$ ./qpgentools match-plp-geno --help
[bin/qpgentools match-plp-geno] -- Check the concordance between pileup and genotypes

 Copyright (c) 2009-2024 by Hyun Min Kang and Adrian Tan
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input pileup ==
   --plp         [STR: ]             : Input pileup file (TSV format)
   --var         [STR: ]             : Input variant file associated with pileups

== Input genotypes ==
   --pfile       [STR: ]             : Input PLINK 2.0 file prefix (PGEN format)
   --bfile       [STR: ]             : Input PLINK 1.9 file prefix (BED format)
   --sample-list [STR: ]             : Comma-separated list of samples to be included (optional)
   --sample-file [STR: ]             : File containing the list of samples to be included (optional)

== Output options ==
   --out         [STR: ]             : Output prefix

== Other options ==
   --max-depth   [INT: 100]          : Maximum depth per site to calculate the likelihood (default: 100)
   --top-k       [INT: 10]           : Number of top matching samples to report in the summary file
   --skip-full   [FLG: OFF]          : Skip writing the full report containing all samples


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
