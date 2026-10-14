# Quickstart for qgenlib

## Installing qgenlib

Before installing `qgenlib`, you need to install
[htslib](https://github.com/samtools/htslib) in the **same directory** you
want to install `qgenlib` (i.e. `qgenlib` and `htslib` should be
siblings directories). You also need [cmake](https://cmake.org/) installed in your system.

After installing `htslib`, you can clone the current snapshot of this repository to install as well

```sh
$ git clone https://github.com/hyunminkang/qgenlib.git
$ cd qgenlib
$ mkdir build
$ cd build
$ cmake ..
$ make

## test the example binary to check if the build was successful
$ ../bin/testqgen --tsv ../test/dummy.dat
```

If you encounter any difficulties, see [Install](install.md) for more details.

## Writing a C++ program using `qgenlib`

The [test program](https://github.com/hyunminkang/qgenlib/blob/main/test_qgen.cpp) accompanied by this library contains a minimal example of how to use `qgenlib` to make your own C++ program. 

A more sophisticated examples can be found in the following example tools:

* [qbcf](https://github.com/hyunminkang/qbcf): A collection of utilities to handle BCF/VCF files