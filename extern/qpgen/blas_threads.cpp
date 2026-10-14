#ifndef _GNU_SOURCE
#define _GNU_SOURCE // RTLD_DEFAULT on glibc
#endif
#include "blas_threads.h"
#include <cstdlib>
#include <cstring>
#include <string>

#ifdef QPGEN_USE_LAPACK
#include <dlfcn.h>
#endif

namespace blas_threads {

#ifdef QPGEN_USE_LAPACK

// Library-specific entry points, looked up among the loaded symbols. Only the
// ones belonging to the linked backend resolve.
template <typename F>
static F sym(const char* name) { return reinterpret_cast<F>(dlsym(RTLD_DEFAULT, name)); }

enum Vendor { FLEXIBLAS, OPENBLAS, MKL, BLIS, ACCELERATE, UNKNOWN };

static Vendor vendor() {
    // FlexiBLAS first: it wraps (and may itself load) another backend
    if (sym<void*>("flexiblas_set_num_threads"))  return FLEXIBLAS;
    if (sym<void*>("openblas_set_num_threads"))   return OPENBLAS;
    if (sym<void*>("MKL_Set_Num_Threads"))        return MKL;
    if (sym<void*>("bli_thread_set_num_threads")) return BLIS;
#ifdef __APPLE__
    return ACCELERATE;
#else
    return UNKNOWN;
#endif
}

std::string backend() {
    switch (vendor()) {
    case FLEXIBLAS: {
        std::string s = "FlexiBLAS";
        typedef int (*F)(char*, size_t);
        if (F f = sym<F>("flexiblas_current_backend")) {
            char buf[256] = {0};
            f(buf, sizeof(buf) - 1);
            if (buf[0]) s += std::string(" (backend: ") + buf + ")";
        }
        return s;
    }
    case OPENBLAS: {
        typedef char* (*F)();
        F f = sym<F>("openblas_get_config");
        return std::string("OpenBLAS") + (f ? std::string(" ") + f() : "");
    }
    case MKL: {
        typedef void (*F)(char*, int);
        if (F f = sym<F>("MKL_Get_Version_String")) {
            char buf[256] = {0};
            f(buf, (int)sizeof(buf) - 1);
            return std::string(buf);
        }
        return "Intel MKL";
    }
    case BLIS: {
        typedef const char* (*F)();
        F f = sym<F>("bli_info_get_version_str");
        return std::string("BLIS") + (f ? std::string(" ") + f() : "");
    }
    case ACCELERATE:
        return "Apple Accelerate";
    default:
#ifdef QPGEN_LAPACK_LIBS
        return std::string("unrecognized LAPACK/BLAS (linked: ") + QPGEN_LAPACK_LIBS + ")";
#else
        return "unrecognized LAPACK/BLAS";
#endif
    }
}

std::string threads_str() {
    switch (vendor()) {
    case FLEXIBLAS: {
        typedef int (*F)();
        F f = sym<F>("flexiblas_get_num_threads");
        return f ? std::to_string(f()) : "unknown";
    }
    case OPENBLAS: {
        typedef int (*F)();
        F f = sym<F>("openblas_get_num_threads");
        return f ? std::to_string(f()) : "unknown";
    }
    case MKL: {
        typedef int (*F)();
        F f = sym<F>("MKL_Get_Max_Threads");
        return f ? std::to_string(f()) : "unknown";
    }
    case BLIS: {
        typedef long (*F)();
        F f = sym<F>("bli_thread_get_num_threads");
        return f ? std::to_string(f()) : "unknown";
    }
    case ACCELERATE: {
        const char* v = getenv("VECLIB_MAXIMUM_THREADS");
        return (v && *v) ? std::string(v) : "auto";
    }
    default: {
        const char* v = getenv("OMP_NUM_THREADS");
        return (v && *v) ? std::string(v) : "unknown";
    }
    }
}

bool set_threads(int n) {
    if (n < 1) return false;
    switch (vendor()) {
    case FLEXIBLAS: {
        typedef void (*F)(int);
        if (F f = sym<F>("flexiblas_set_num_threads")) { f(n); return true; }
        return false;
    }
    case OPENBLAS: {
        typedef void (*F)(int);
        if (F f = sym<F>("openblas_set_num_threads")) { f(n); return true; }
        return false;
    }
    case MKL: {
        typedef void (*F)(int);
        if (F f = sym<F>("MKL_Set_Num_Threads")) { f(n); return true; }
        return false;
    }
    case BLIS: {
        typedef void (*F)(long);
        if (F f = sym<F>("bli_thread_set_num_threads")) { f((long)n); return true; }
        return false;
    }
    case ACCELERATE: {
        // Accelerate reads VECLIB_MAXIMUM_THREADS at its first BLAS call, so
        // setting it in-process works as long as no BLAS call has happened yet.
        // BLASSetThreading (macOS 15+) additionally forces single-threading.
        setenv("VECLIB_MAXIMUM_THREADS", std::to_string(n).c_str(), 1);
        typedef int (*F)(unsigned int);
        if (F f = sym<F>("BLASSetThreading"))
            f(n == 1 ? 1u /*SINGLE_THREADED*/ : 0u /*MULTI_THREADED*/);
        return true;
    }
    default: {
        // unknown backend: try OpenMP, which most threaded BLAS builds use
        setenv("OMP_NUM_THREADS", std::to_string(n).c_str(), 1);
        typedef void (*F)(int);
        if (F f = sym<F>("omp_set_num_threads")) { f(n); return true; }
        return false;
    }
    }
}

#else // !QPGEN_USE_LAPACK

std::string backend() { return "none (Eigen built-in, single-threaded)"; }
std::string threads_str() { return "1"; }
bool set_threads(int) { return false; }

#endif

} // namespace blas_threads
