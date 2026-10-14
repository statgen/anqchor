# qgenlib
Libraries for quantitative genomics

## Overview

`qgenlib` is a collection of C++ libraries (and test binaries) to help quantitative genomic analysis.

## Installing qgenlib

Before installing `qgenlib`, you need to install
[htslib](https://github.com/samtools/htslib) in the same directory you
want to install `qgenlib' (i.e. `qgenlib` and `htslib` should be
siblings directories). You also need [cmake](https://cmake.org/) installed in your system.

After installing `htslib`, you can clone the current snapshot of this repository to install as well

<pre>
$ mkdir build

$ cd build

$ cmake ..
</pre>

In case any required libraries is missing, you may specify customized installing path by replacing "cmake .." with:

<pre>
For libhts:
  - $ cmake -DHTS_INCLUDE_DIRS=/htslib_absolute_path/include/  -DHTS_LIBRARIES=/htslib_absolute_path/lib/libhts.a ..

For bzip2:
  - $ cmake -DBZIP2_INCLUDE_DIRS=/bzip2_absolute_path/include/ -DBZIP2_LIBRARIES=/bzip2_absolute_path/lib/libbz2.a ..

For lzma:
  - $ cmake -DLZMA_INCLUDE_DIRS=/lzma_absolute_path/include/ -DLZMA_LIBRARIES=/lzma_absolute_path/lib/liblzma.a ..
</pre>

Finally, to build the binary, run

<pre>
$ make
</pre>

### List of test applications contained in `qgenlib`

`qgenlib` contains many in-house C++ tools that are currently under
the hood development phase. To list the available commands of tools, type:

<pre>
testqgen --help
</pre>

To see the usage of individual commands, type:

<pre>
testqgen [command] --help
</pre>
