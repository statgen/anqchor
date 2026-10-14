#ifndef __BLAS_THREADS_H
#define __BLAS_THREADS_H

// Identify and control the system LAPACK/BLAS that qpgen was linked against
// (used by SuSiE-inf/ash for the genotype cross-product and eigendecomposition).
//
// The backend is detected at RUN time from the symbols actually loaded, so it
// is correct even when a generic libblas/liblapack is a symlink to OpenBLAS,
// MKL, BLIS or FlexiBLAS (e.g. Debian alternatives, EasyBuild modules).

#include <string>

namespace blas_threads {

// Human-readable backend, e.g. "Apple Accelerate", "OpenBLAS 0.3.26 ...",
// "Intel MKL 2024.1 ...", or "none (Eigen built-in, single-threaded)".
std::string backend();

// Current thread setting as a string ("4", "auto", "1", or "unknown").
std::string threads_str();

// Set the number of BLAS/LAPACK threads (n >= 1). Must be called before the
// first BLAS call (Accelerate reads its setting once). Returns false, leaving
// the library default, if the backend offers no way to control it.
bool set_threads(int n);

} // namespace blas_threads

#endif
