# qpgentools geno2tsv

## Summary

`qpgentools geno2tsv` extracts genotypes from a PLINK/PLINK2 file for a specific list of variants and samples, and outputs them in a TSV format. This is useful for exporting a subset of genotype data for downstream analysis or manual inspection.

A typical running example command is given below:

```bash
qpgentools geno2tsv --pfile [pgen_prefix] --var [variant_file] --out [output_prefix]
```

## Input File Formats

### Variant File (`--var`)

The variant file is a TSV file that defines the variants to be extracted.
-   The first column must contain the **variant index** in the genotype file.
-   Additional columns are optional and will be included in the output "as is" to describe the variant.
-   If the file has a header (starts with `#`), the header will be preserved in the output.

## Required options

*   **Input Genotypes** (One of the following):
    *   `--pfile`: Input PLINK 2.0 file prefix (PGEN format).
    *   `--bfile`: Input PLINK 1.9 file prefix (BED format).
*   **Input Variant Info**:
    *   `--var`: Input variant file (TSV) containing `[variant index] [extra columns]`.
*   **Sample Selection**:
    *   `--sample-list`: Comma-separated list of sample IDs to be included.
    *   `--sample-file`: File containing the list of samples to be included.
    *   (Note: While the code checks for sample selection, if neither is provided, it may error or default depending on implementation logic, but the help message implies they are optional for filtering, yet the code enforces one of them to be present).
*   `--out`: Output prefix.

## Additional Options

*   `--idx-begin` : Starting index for the variant in the input `--var` file (default: 0). Use 1 if your variant file uses 1-based indexing.
*   `--sparse` : Output sparse representation `[IDX]:[GENO]` instead of full genotypes. This is useful for very large datasets where most genotypes are reference (0).

## Expected Output

The output is a TSV file containing:

1.  **Variant Columns**: All columns from the input `--var` file except the first (index) column.
2.  **AC**: Allele Count of the alternative allele within the selected samples.
3.  **AN**: Allele Number (total alleles called, usually 2 * N_samples) within the selected samples.
4.  **Sample Columns**: One column per selected sample.
    *   **Default Format**: `0` (HomRef), `1` (Het), `2` (HomAlt), `NA` (Missing). **Note**: The code outputs allele counts for the *alternative* allele, so `0`=Ref/Ref, `1`=Ref/Alt, `2`=Alt/Alt.
    *   **Sparse Format** (`--sparse`): `[SampleIndex]:[Genotype]`. Only non-reference genotypes are printed.
        *   `[i]:1` : Heterozygous
        *   `[i]:2` : Homozygous Alternative
        *   `[i]:NA` : Missing
        *   Homozygous Reference genotypes are omitted.

## Full Usage

The full usage of `qpgentools geno2tsv` can be viewed with the `--help` option:

```
$ ./qpgentools geno2tsv --help
[bin/qpgentools geno2tsv] -- Extract genotypes into TSV format

 Copyright (c) 2009-2024 by Hyun Min Kang and Adrian Tan
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input variant info ==
   --var         [STR: ]             : Input variant file (TSV) containing [variant index] [extra columns]
   --idx-begin   [INT: 0]            : Starting index for the variant (default: 0). Use 1 for 1-based index

== Input genotypes ==
   --pfile       [STR: ]             : Input PLINK 2.0 file prefix (PGEN format)
   --bfile       [STR: ]             : Input PLINK 1.9 file prefix (BED format)
   --sample-list [STR: ]             : Comma-separated list of samples to be included (optional)
   --sample-file [STR: ]             : File containing the list of samples to be included (optional)

== Output options ==
   --out         [STR: ]             : Output prefix
   --sparse      [FLG: OFF]          : Output sparse representation [IDX]:[GENO] instead of full genotypes


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
