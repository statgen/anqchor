cmake_minimum_required(VERSION 3.10) # 2.8 is very old; 3.10+ is recommended for modern targets
project(qgenlib)

set(CMAKE_CXX_STANDARD 14)
# Don't use -O3 or -std in flags globally; use target_compile_options if needed.
# This prevents overriding the parent project's flags.

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")

# --- 1. Find HTSLib (Keep this if it's custom, but standard paths usually work better) ---
find_path(HTS_INCLUDE_DIRS htslib/bgzf.h HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ ../htslib/)
find_library(HTS_LIBRARIES NAMES libhts.a hts HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ ../htslib/)

if(NOT HTS_INCLUDE_DIRS OR NOT HTS_LIBRARIES)
    message(FATAL_ERROR "libhts not found")
endif()
include_directories(${HTS_INCLUDE_DIRS})

# --- 2. Define Sources ---
set(SOURCE_FILES
    genome_interval.cpp gtf.cpp hts_utils.cpp params.cpp phred_helper.cpp
    qgen_error.cpp qgen_utils.cpp tsv_reader.cpp commands.cpp dataframe.cpp
    # Headers usually don't need to be listed in SOURCE_FILES for compilation, 
    # but it helps IDEs display them.
)

ADD_LIBRARY(qgen STATIC ${SOURCE_FILES})

# --- 3. Use Standard Packages (The Fix) ---

# Replace manual find_library with find_package for ZLIB
find_package(ZLIB REQUIRED)

# Replace manual find_library with find_package for CURL
find_package(CURL REQUIRED)

# Replace manual find_library with find_package for OpenSSL (Crypto)
find_package(OpenSSL REQUIRED)

# BZip2 (You already had this correct!)
find_package(BZip2 REQUIRED)
if(BZIP2_FOUND)
    include_directories(${BZIP2_INCLUDE_DIRS})
endif()

# LZMA (Still manual because it's less standard, but okay)
find_library(LZMA lzma HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)

# --- 4. Link Libraries ---
SET(APP_EXE testqgen)
ADD_EXECUTABLE(${APP_EXE} test_qgen.cpp)

# --- Find libdeflate dynamically ---
find_library(DEFLATE_LIBRARY NAMES deflate libdeflate HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ /opt/homebrew/lib)

if(DEFLATE_LIBRARY)
    message(STATUS "Found libdeflate: ${DEFLATE_LIBRARY}")
else()
    message(STATUS "libdeflate not found. Skipping.")
    # Set to empty so target_link_libraries doesn't try to link a "NOTFOUND" variable
    set(DEFLATE_LIBRARY "") 
endif()

# Use the variables provided by find_package (ZLIB::ZLIB, CURL::libcurl, OpenSSL::Crypto)
# These automatically handle include dirs and linking flags.
target_link_libraries(${APP_EXE} 
    qgen 
    ${HTS_LIBRARIES} 
    ZLIB::ZLIB        # Better than ${ZLIB}
    CURL::libcurl     # Better than ${CURLLIB}
    OpenSSL::Crypto   # Better than ${CRYPTOLIB}
    ${BZIP2_LIBRARIES} 
    ${LZMA}
    ${DEFLATE_LIBRARY} # Now conditionally linked!
)

install(TARGETS ${APP_EXE} RUNTIME DESTINATION bin)
