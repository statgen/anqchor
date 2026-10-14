# Example Usage of qgenlib

## Coding with qgenlib

`qgenlib` provides C++ libraries for following functionalities.

* [Command Line Arguments](examples/params.md) A macro-based command line parsing library implemented by Gonçalo Abecasis and used in many [tools](https://github.com/statgen).
* [I/O Functions](examples/io.md) A collection of helper functions to output messages, and handle basic formats supported by [htslib](https://www.htslib.org/).
* [Handling TSV Files](examples/tsv.md) Classes and libraries to aid basic data wrangling using tab-limited files.
* [Genomics](examples/genomics.md) Simple classes to handle simple genomic data types.
* [Miscellaneous](examples/misc.md) Other utility functions to handle strings.

Note that these examples contain only a subset of functions that `qgenlib` provides. A complete list of API reference may be available later. 

## Building a C++ program using qgenlib

A recommended way to build a new C++ program is to place your new program in parallel to the `qgenlib` and `htslib` directories. The following is an example of the directory structure.

```sh
/path/to/install
|-- htslib
|   `-- libhts.a
|-- qgenlib
|   `-- lib
|       `-- libqgen.a
`-- your_new_program
```

An example program that can be tested this way can be found at [qbcf](https://github.com/hyunminkang/qbcf)


