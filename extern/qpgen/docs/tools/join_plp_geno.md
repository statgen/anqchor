# qpgentools join-plp-geno

## Summary

`qpgentools join-plp-geno` joins the variant file and pileup file together with the genotypes of a specific sample. The primary purpose of this tool is to perform sanity checking of the input files and to verify whether the pileup was generated reasonably well by comparing the observed reads with the ground truth genotypes.

A typical running example command is given below:

```bash
qpgentools join-plp-geno --plp [pileup_file] --var [variant_file] --pfile [reference_pgen] --id [sample_id] --out [output_file]
```

## Input File Formats

The input file formats for `--var` and `--plp` are identical to those used in [`match-plp-geno`](match_plp_geno.md).

### Variant File (`--var`)

Expected columns (TSV with header):
1.  **Variant ID**: Variant ID
2.  **CHROM**: Chromosome
3.  **POS**: Position (1-based)
4.  **REF**: Reference allele
5.  **ALT**: Alternative allele

### Pileup File (`--plp`)

Expected columns (TSV with header):
1.  **Barcode**: Cell barcode or sample ID
2.  **Variant Index**: 0-based index of the variant in the `--var` file
3.  **Depth**: Read depth at this site
4.  **Alleles**: String representing alleles observed in reads (`0` for REF, `1` for ALT)
5.  **Base Qualities**: String representing base qualities for each read

## Required options

*   **Input Pileup**:
    *   `--plp`: Input pileup file (TSV format).
    *   `--var`: Input variant file associated with pileups.
*   **Input Genotypes** (One of the following):
    *   `--pfile`: Input PLINK 2.0 file prefix (PGEN format).
    *   `--bfile`: Input PLINK 1.9 file prefix (BED format).
*   **Sample Selection**:
    *   `--id`: Sample ID from the genotype file that should be used for joining.
*   `--out`: Output file name.

## Additional Options

*   `--het-only` : Output only heterozygous genotypes (genotype 1). Useful for checking reference bias.

## Expected Output

The output is a TSV file containing the joined information:

*   `#Barcode`: Barcode/Sample ID from the pileup file.
*   `Variant`: 0-based variant index.
*   `Genotype`: Genotype of the selected sample (`--id`) at this variant.
    *   `0`: Homozygous Alternative (AA) - *Note: qpgentools internal encoding*
    *   `1`: Heterozygous (RA)
    *   `2`: Homozygous Reference (RR)
    *   `-1`: Missing or other
*   `Depth`: Read depth from the pileup.
*   `nR`: Number of reference alleles observed in the reads.
*   `nA`: Number of alternative alleles observed in the reads.
*   `nO`: Number of other alleles observed in the reads.
*   `Allele`: String of alleles observed (from pileup).
*   `BQ`: String of base qualities (from pileup).

## Full Usage

The full usage of `qpgentools join-plp-geno` can be viewed with the `--help` option:

```
$ ./qpgentools join-plp-geno --help
[bin/qpgentools join-plp-geno] -- Join the pileup and individual genotypes

 Copyright (c) 2009-2024 by Hyun Min Kang and Adrian Tan
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input pileup ==
   --plp      [STR: ]             : Input pileup file (TSV format)
   --var      [STR: ]             : Input variant file associated with pileups

== Input genotypes ==
   --pfile    [STR: ]             : Input PLINK 2.0 file prefix (PGEN format)
   --bfile    [STR: ]             : Input PLINK 1.9 file prefix (BED format)
   --id       [STR: ]             : Sample ID that should be used for joining

== Output options ==
   --out      [STR: ]             : Output file that joins the pileup and genotypes

== Other options ==
   --het-only [FLG: OFF]          : Output only heterozygous genotypes


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
