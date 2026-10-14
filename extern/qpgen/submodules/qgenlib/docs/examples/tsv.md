# Handling TSV Files

## Introduction

`qgenlib` provides functions for reading tab-limited files to handle tidy data types easily through `tsv_reader` and `dataframe_t` classes.

## Streamlined reading of tab-limited files

The `tsv_reader` class offers a streamlined reading of tab-delimited file, in plaintext format or gzip- or bgzip-compressed formats. For bgzipped-compressed files, it also offers the ability to read from a specific region of the file through `tabix` index files.

Below are examples of using `tsv_reader`:

```cpp
// streamline reading of tab-delimited files
// this example read the second column of the file and sum up the values
tsv_reader tr1("input.tsv.gz");
double sum = 0.0;
while( tr1.read_line() ) {
    sum += tr1.get_double_at(1); // sum up the values in the second column
}
printf("Sum of the second column: %lf\n", sum);

// read a specific region of a bgzipped-compressed file
tsv_reader tr2("input.tsv.gz");
tr2.jump_to("chr1:1000000-2000000"); // read a specific region of the file
sum = 0.0;
while( tr2.read_line() ) {
    sum += tr2.get_double_at(1); // sum up the values in the second column
}
printf("Sum of the second column: %lf\n", sum);
``` 

The class `tsv_reader` would be particularly useful when writing a big tsv file that may not fit into the memory.

## Loading a tab-limited files into a `dataframe_t` object

`dataframe_t` class is intended to handle  tab-limited files as a dataframe, which is a collection of columns. It reads the entire file into memory and provides a convenient way to access the data, so it is recommended to use it only for small files.

Below are examples of using `dataframe_t`:

```cpp
// Read a tsv file into a dataframe_t object
dataframe_t df("input.tsv.gz");
if ( df.has_column("value") ) { // check if the column 'value' exists
    int idx = df.get_colidx("value");
    double sum = 0.0;
    for(int i = 0; i < df.nrows; i++) {
        sum += df.get_double_elem(i, idx);
    }
    // print the sum of the value column
    printf("Sum of the value column: %lf\n", sum);
}
else { // if the column 'value' does not exist, print an error message
    error("The column 'value' does not exist in the file");
}
```