#include "susie_export.h"
#include "qgenlib/qgen_error.h"
#include "htslib/bgzf.h"
#include <cmath>
#include <cstring>

namespace susie_export {

namespace {

// Thin BGZF writer with the primitive types used by the format.
class BgzfWriter {
public:
    explicit BgzfWriter(const char* path) : path_(path) {
        fp_ = bgzf_open(path, "w");
        if ( fp_ == NULL ) error("Cannot open %s for writing", path);
    }
    ~BgzfWriter() { if ( fp_ ) bgzf_close(fp_); }

    void raw(const void* data, size_t nbytes) {
        // bgzf_write splits large buffers into blocks internally; write in
        // moderately sized pieces anyway so ssize_t never overflows.
        const char* p = static_cast<const char*>(data);
        const size_t step = (size_t)1 << 30;
        while ( nbytes > 0 ) {
            size_t len = nbytes < step ? nbytes : step;
            if ( bgzf_write(fp_, p, len) < 0 ) error("Failed to write to %s", path_.c_str());
            p += len; nbytes -= len;
        }
    }
    void magic(const char* m8) { raw(m8, 8); }
    void i32(int32_t v) { raw(&v, sizeof(v)); }
    void str(const std::string& s) { raw(s.c_str(), s.size() + 1); } // NUL-terminated
    void dbl(double v) { raw(&v, sizeof(v)); }
    void dbls(const Eigen::MatrixXd& m) { // column-major, matches R's matrix()
        raw(m.data(), sizeof(double) * (size_t)m.rows() * (size_t)m.cols());
    }
    void dbls(const Eigen::VectorXd& v) { raw(v.data(), sizeof(double) * (size_t)v.size()); }
    void i32s(const std::vector<int32_t>& v) { if ( !v.empty() ) raw(v.data(), sizeof(int32_t) * v.size()); }

    void close() { if ( fp_ ) { if ( bgzf_close(fp_) < 0 ) error("Failed to close %s", path_.c_str()); fp_ = NULL; } }

private:
    std::string path_;
    BGZF* fp_;
};

void write_header(BgzfWriter& w, const char* magic, const RegionMeta& meta, int32_t p, int32_t K) {
    if ( (int32_t)meta.variant_ids.size() != p || (int32_t)meta.positions.size() != p || (int32_t)meta.af.size() != p )
        error("susie_export: variant metadata size does not match p=%d", p);
    if ( (int32_t)meta.trait_ids.size() != K )
        error("susie_export: %d trait IDs but K=%d", (int32_t)meta.trait_ids.size(), K);
    w.magic(magic);
    w.i32(2);                 // format version
    w.i32(meta.n_samples);
    w.i32(p);
    w.i32(K);
    w.i32(meta.n_cov);
    w.i32(meta.test);
    w.str(meta.region);
    for(int32_t k = 0; k < K; ++k) w.str(meta.trait_ids[k]);
    for(int32_t j = 0; j < p; ++j) w.str(meta.variant_ids[j]);
    w.i32s(meta.positions);
    Eigen::VectorXd af = Eigen::Map<const Eigen::VectorXd>(meta.af.data(), p);
    w.dbls(af);
}

} // anonymous namespace

void write_suff_stats(const char* path, const RegionMeta& meta,
                      const Eigen::MatrixXd& XtX, const Eigen::MatrixXd& Xty,
                      const Eigen::VectorXd& yty)
{
    const int32_t p = (int32_t)XtX.rows();
    const int32_t K = (int32_t)Xty.cols();
    if ( XtX.cols() != p || Xty.rows() != p || yty.size() != K )
        error("susie_export: inconsistent dimensions for sufficient statistics");
    BgzfWriter w(path);
    write_header(w, "QPGNSUFF", meta, p, K);
    w.dbls(XtX);
    w.dbls(Xty);
    w.dbls(yty);
    w.close();
}

void write_rss_stats(const char* path, const RegionMeta& meta,
                     const Eigen::MatrixXd& R, const Eigen::MatrixXd& z,
                     const Eigen::MatrixXd& bhat, const Eigen::MatrixXd& shat,
                     const Eigen::VectorXd& var_y)
{
    const int32_t p = (int32_t)R.rows();
    const int32_t K = (int32_t)z.cols();
    if ( R.cols() != p || z.rows() != p || bhat.rows() != p || bhat.cols() != K ||
         shat.rows() != p || shat.cols() != K || var_y.size() != K )
        error("susie_export: inconsistent dimensions for RSS statistics");
    BgzfWriter w(path);
    write_header(w, "QPGN_RSS", meta, p, K);
    w.dbls(R);
    w.dbls(z);
    w.dbls(bhat);
    w.dbls(shat);
    w.dbls(var_y);
    w.close();
}

Eigen::MatrixXd xtx_to_corr(const Eigen::MatrixXd& XtX)
{
    const Eigen::Index p = XtX.rows();
    Eigen::VectorXd inv_sd(p);
    for(Eigen::Index j = 0; j < p; ++j) {
        const double d = XtX(j, j);
        inv_sd(j) = ( std::isfinite(d) && d > 0.0 ) ? 1.0 / std::sqrt(d) : 0.0;
    }
    // Fill the lower triangle and mirror it so R is bit-exactly symmetric
    // (susieR checks symmetry and warns otherwise).
    Eigen::MatrixXd R(p, p);
    for(Eigen::Index j = 0; j < p; ++j) {
        R(j, j) = 1.0; // exact unit diagonal, also for zero-variance columns
        for(Eigen::Index i = j + 1; i < p; ++i) {
            const double r = XtX(i, j) * inv_sd(i) * inv_sd(j);
            R(i, j) = r;
            R(j, i) = r;
        }
    }
    return R;
}

} // namespace susie_export
