# Installing qpgen

## Requirements

* A C++14-capable compiler with OpenMP support (`cmake` configures with `find_package(OpenMP REQUIRED)`).
    * On macOS, the Apple Clang shipped with Xcode does not include OpenMP; install it with `brew install libomp`.
* [cmake](https://cmake.org/) (2.8 or later)
* Standard compression libraries used by `htslib` : `zlib`, `bzip2`, `lzma`, and (optionally) `libcurl`, `libdeflate`, `libcrypto`
* `bgzip` and `tabix` (from [htslib](https://github.com/samtools/htslib)) are needed to build the
  indexed pvar files described in [Genotype file formats](formats/genotypes.md).

`qpgentools` bundles [htslib](https://github.com/samtools/htslib),
[qgenlib](https://github.com/hyunminkang/qgenlib), and [Eigen](https://gitlab.com/libeigen/eigen)
as git submodules, so all three are built from within the repository.

## Installing qpgen

Because `htslib`, `qgenlib`, and `eigen` are submodules, you need to clone the repository recursively, and build the submodules before building the `qpgen` package. An example instruction is given below.

```sh
## STEP 1 : CLONE THE REPOSITORY
## clone the current snapshot of this repository
git clone --recursive https://github.com/hyunminkang/qpgen.git

## move to the qpgen directory
cd qpgen

## STEP 2 : BUILD THE SUBMODULES
## move to the submodules directory
cd submodules

## build the submodules using build.sh script
sh -x build.sh

## move back to the qpgen directory
cd ..

## STEP 3 : BUILD QPGENTOOLS
## create a build directory
mkdir build
cd build

## Run cmake to configure the build
cmake ..

## Build the qpgentools package
make
```

If you cloned the repository without `--recursive`, run `git submodule update --init --recursive`
before STEP 2. `cmake` fails with `Eigen submodule missing` when the `submodules/eigen` directory
has not been populated.

If `cmake` is not found, you need to install [cmake](https://cmake.org/) in your system.

The build produces two artifacts:

* `bin/qpgentools` : the command line executable
* `lib/libqpgen.a` : the static library

To check which build you are running:

```
$ ./bin/qpgentools --version
qpgentools 0.2.0 (git 938bf19)
```

The version is set in `CMakeLists.txt`, and the git commit is read at every build. A `-dirty`
suffix means the build includes uncommitted changes to tracked files; `unknown` means the source
is not a git checkout. Every command also writes this line at the start of its log.

## (Optional) Customized specification of the library path

In case any required libraries is missing in `cmake`, you may specify customized installing path by replacing "cmake .." with the following options:

```sh
## If qgenlib is missing or installed in a different directory
$ cmake -DQGEN_INCLUDE_DIRS=/qgenlib_absolute_path/include \
        -DQGEN_LIBRARIES=/qgenlib_absolute_path/lib/libqgen.a ..

## If htslib is missing or installed in a different directory
$ cmake -DHTS_INCLUDE_DIRS=/htslib_absolute_path/include/  \
        -DHTS_LIBRARIES=/htslib_absolute_path/libhts.a ..

## If bzip2 is missing or installed in a different directory
$ cmake -DBZIP2_INCLUDE_DIRS=/bzip2_absolute_path/include/ \
        -DBZIP2_LIBRARIES=/bzip2_absolute_path/lib/libbz2.a ..

## If lzma is missing or installed in a different directory
$ cmake -DLZMA_INCLUDE_DIRS=/lzma_absolute_path/include/ \
        -DLZMA_LIBRARIES=/lzma_absolute_path/lib/liblzma.a ..

## You may combine the multiples options above if needed.
## Other missing libraries can be handled in a similar way.
```

## Testing qpgen

To test whether build was successful, you can run the following command:

```sh
## Current directory: /path/to/install/qpgen/build
$ ../bin/qpgentools --help
```
