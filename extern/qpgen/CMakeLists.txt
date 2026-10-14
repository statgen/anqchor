cmake_minimum_required(VERSION 2.8...3.23)
project(qpgen VERSION 0.2.0)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread -std=c++14")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")

find_path(HTS_INCLUDE_DIRS htslib/bgzf.h ../htslib/bgzf.h HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ submodules/htslib ../htslib)
if(NOT HTS_INCLUDE_DIRS )
    message(FATAL_ERROR "libhts HTS_INCLUDE_DIRS not found")
endif()

message(STATUS ${HTS_INCLUDE_DIRS})

find_library(HTS_LIBRARIES libhts.a HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ submodules/htslib ../htslib/)
if(NOT HTS_LIBRARIES)
    message(FATAL_ERROR "libhts HTS_LIBRARIES not found")
endif()

message(STATUS ${HTS_LIBRARIES})

find_path(QGEN_INCLUDE_DIRS qgenlib/qgen_error.h  HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ submodules/qgenlib ../qgenlib)
if(NOT QGEN_INCLUDE_DIRS )
    message(FATAL_ERROR "libqgen QGEN_INCLUDE_DIRS not found")
endif()

message(STATUS ${QGEN_INCLUDE_DIRS})

find_library(QGEN_LIBRARIES libqgen.a  HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ submodules/qgenlib/lib/ ../qgenlib/lib/)
if(NOT QGEN_LIBRARIES)
    message(FATAL_ERROR "libqgen QGEN_LIBRARIES not found")
endif()

message(STATUS ${QGEN_LIBRARIES})

# Eigen
set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/submodules/eigen")
if (NOT EXISTS "${EIGEN3_INCLUDE_DIR}/Eigen/Dense")
    message(FATAL_ERROR "Eigen submodule missing. Did you run 'git submodule update --init'?")
endif()
add_library(Eigen3::Eigen INTERFACE IMPORTED)
set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")


include_directories(. simde/ ${HTS_INCLUDE_DIRS} ${QGEN_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})

if(APPLE)
    message(STATUS "Configuring for macOS")
    # Prefer Homebrew LLVM: the unversioned keg first, then the newest versioned keg (llvm@NN).
    file(GLOB _LLVM_DIRS /opt/homebrew/opt/llvm /opt/homebrew/opt/llvm@*)
    list(SORT _LLVM_DIRS ORDER DESCENDING)
    set(_LLVM_FOUND FALSE)
    foreach(_dir IN LISTS _LLVM_DIRS)
        if(EXISTS "${_dir}/bin/clang++" AND NOT _LLVM_FOUND)
            set(CMAKE_C_COMPILER "${_dir}/bin/clang")
            set(CMAKE_CXX_COMPILER "${_dir}/bin/clang++")
            message(STATUS "Using Homebrew LLVM at ${_dir}")
            set(_LLVM_FOUND TRUE)
        endif()
    endforeach()
    if(NOT _LLVM_FOUND)
        message(WARNING "Homebrew LLVM not found under /opt/homebrew/opt; using default compiler ${CMAKE_CXX_COMPILER}")
    endif()
endif()


set(APP_FILES
    qpgentools.cpp
    qpgentools.h
    cmd_test_qpgen.cpp
    cmd_match_plp_geno.cpp
    cmd_join_plp_geno.cpp
    cmd_geno2tsv.cpp
    cmd_pgen2tsv.cpp
    cmd_pair_assoc.cpp
    cmd_pair_assoc_v0.cpp
    cmd_pair_prs.cpp
    cmd_rect_assoc.cpp
    cmd_match_prs_pheno.cpp
    cmd_region_assoc.cpp
)

set(SOURCE_FILES
    qpgen.cpp
    qpgen.h
    pheno.cpp
    pheno.h
    qpgen_utils.cpp
    qpgen_utils.h
    assoc_utils.h
    assoc_utils.cpp
    susie_utils.h
    susie_utils.cpp
    susie_export.h
    susie_export.cpp
    blas_threads.h
    blas_threads.cpp
    include/pgenlib_misc.cc
    include/pgenlib_misc.h
    include/pgenlib_read.cc
    include/pgenlib_read.h
    include/plink2_base.cc
    include/plink2_base.h
    include/plink2_bits.cc
    include/plink2_bits.h
    include/plink2_simd.cc
    include/plink2_simd.h
    include/pgenlib_ffi_support.cc
    include/pgenlib_ffi_support.h
    pgenlibr.cpp
    pgenlibr.h
    pvar_ffi_support.cc
    pvar_ffi_support.h
    # SIMDe headers are header-only and resolved via include_directories(simde/);
    # not enumerated here so the build is robust to SIMDe version updates.
)

ADD_LIBRARY(qpgen STATIC ${SOURCE_FILES})

SET(APP_EXE qpgentools)

ADD_EXECUTABLE(${APP_EXE} ${SOURCE_FILES} ${APP_FILES})

# Version string (qpgentools --version, and the first log line of every command): the project
# version above plus the git commit, regenerated at every build so it is never stale.
set(QPGEN_VERSION_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated")
add_custom_target(qpgen_version
    COMMAND ${CMAKE_COMMAND} -DSRC_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -DOUT_FILE=${QPGEN_VERSION_DIR}/qpgen_version.h
            -DQPGEN_VERSION=${PROJECT_VERSION}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/qpgen_version.cmake
    BYPRODUCTS ${QPGEN_VERSION_DIR}/qpgen_version.h
    COMMENT "Updating qpgen version information")
add_dependencies(${APP_EXE} qpgen_version)
target_include_directories(${APP_EXE} PRIVATE ${QPGEN_VERSION_DIR})

find_library(ZLIB z HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT ZLIB)
    message(FATAL_ERROR "libz library not found")
endif()

find_library(CURLLIB curl HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CURLLIB)
    message(FATAL_ERROR "libcurl library not found")
endif()

find_package (BZip2)
if (NOT BZIP2_FOUND)
    message(FATAL_ERROR "libbz2 library not found")
else()
    include_directories(${BZIP2_INCLUDE_DIRS})
    set(BZIP2 ${BZIP2_LIBRARIES})
endif()

find_library(LZMA lzma HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT LZMA)
    message(WARNING "liblzma library not found, if you specified --disable-lzma when compiling libhts, please ignore this warning")
    set(LZMA "")
else()
    set(LZMA ${LZMA_LIBRARIES})
endif()

find_library(DEFLATELIB deflate HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT DEFLATELIB)
    message(FATAL_ERROR "libdeflate library not found")
endif()

find_library(CRYPTOLIB crypto HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CRYPTOLIB)
    message(FATAL_ERROR "libcrypto library not found")
endif()

if(APPLE)
    # Check for Homebrew install path on Apple Silicon
    if(EXISTS "/opt/homebrew/opt/libomp")
        set(OpenMP_ROOT "/opt/homebrew/opt/libomp")
    # Check for Homebrew install path on Intel Mac
    elseif(EXISTS "/usr/local/opt/libomp")
        set(OpenMP_ROOT "/usr/local/opt/libomp")
    endif()
endif()
find_package(OpenMP REQUIRED)

# Optional LAPACK/BLAS (Accelerate on macOS; OpenBLAS/MKL/reference on Linux,
# selectable with -DBLA_VENDOR=...). Used by SuSiE-inf/ash for the genotype
# cross-product and its eigendecomposition, which are multithreaded and several
# times faster than Eigen's built-in single-threaded routines. Without it the
# code falls back to Eigen. The backend in use is detected at run time and
# logged; region-assoc --threads sets its thread count (see blas_threads.h).
option(QPGEN_USE_LAPACK "Use system LAPACK/BLAS for SuSiE linear algebra" ON)
if(QPGEN_USE_LAPACK)
    find_package(LAPACK)
    if(LAPACK_FOUND)
        message(STATUS "LAPACK found: ${LAPACK_LIBRARIES}")
        # the linked library list is embedded so region-assoc can report it
        string(REPLACE ";" " " QPGEN_LAPACK_LIBS_STR "${LAPACK_LIBRARIES}")
        foreach(tgt qpgen ${APP_EXE})
            target_compile_definitions(${tgt} PRIVATE QPGEN_USE_LAPACK
                "QPGEN_LAPACK_LIBS=\"${QPGEN_LAPACK_LIBS_STR}\"")
        endforeach()
        target_link_libraries(${APP_EXE} PRIVATE ${LAPACK_LIBRARIES} ${CMAKE_DL_LIBS})
    else()
        message(WARNING "LAPACK not found: SuSiE-inf/ash will use single-threaded Eigen linear algebra (slower for large regions)")
    endif()
endif()

#TARGET_LINK_LIBRARIES(${APP_EXE} PRIVATE ${QGEN_LIBRARIES} ${HTS_LIBRARIES} ${ZLIB} ${LZMA} ${BZIP2} ${CURLLIB} ${DEFLATELIB} ${CRYPTOLIB} "${OpenMP_CXX_FLAGS}")
#TARGET_LINK_LIBRARIES(${APP_EXE} PRIVATE Eigen3::Eigen ${QGEN_LIBRARIES} ${HTS_LIBRARIES} ${ZLIB} ${LZMA} ${BZIP2} ${CURLLIB} ${DEFLATELIB} ${CRYPTOLIB} "${OpenMP_CXX_FLAGS}")
TARGET_LINK_LIBRARIES(${APP_EXE} PRIVATE Eigen3::Eigen ${QGEN_LIBRARIES} ${HTS_LIBRARIES} ${ZLIB} ${LZMA} ${BZIP2} ${CURLLIB} ${DEFLATELIB} ${CRYPTOLIB} "${OpenMP_CXX_FLAGS}")
#target_compile_options(${APP_EXE} PRIVATE "${OpenMP_CXX_FLAGS}")


install(TARGETS ${APP_EXE} RUNTIME DESTINATION bin)
