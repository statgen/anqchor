## qpgen_susie_io.R -- read the binary summary-data files written by
##   qpgentools region-assoc --out-suff / --out-rss
## and hand them to susieR (https://github.com/stephenslab/susieR).
##
## Usage:
##   source("scripts/qpgen_susie_io.R")
##   ss  <- read_qpgen_suff("out.suff.bin.gz")   # X'X, X'y, y'y, n
##   rss <- read_qpgen_rss("out.rss.bin.gz")     # z, R, n, bhat, shat, var_y, test
##   fit <- susie_from_qpgen(ss,  trait = "ENSG00000164308.17", L = 10)  # susie_ss(XtX, Xty, yty, n)
##   fit <- susie_from_qpgen(rss, trait = "ENSG00000164308.17", L = 10)  # susie_rss(bhat, shat, R, n, var_y)
##   fit$sets$cs; susieR::susie_get_pip(fit)
##
## Both files are gzip (BGZF) streams readable with gzfile(); no extra
## packages are needed to read them. Matrices are stored column-major in
## little-endian float64, so readBin() + matrix() reconstructs them directly.

.qpgen_read_header <- function(con, expected_magic) {
  magic <- readChar(con, 8, useBytes = TRUE)
  if (!identical(magic, expected_magic))
    stop(sprintf("unexpected magic '%s' (expected '%s')", magic, expected_magic))
  version <- readBin(con, "integer", 1, size = 4, endian = "little")
  if (!(version %in% c(1L, 2L))) stop("unsupported format version ", version)
  dims <- readBin(con, "integer", 4, size = 4, endian = "little")
  n <- dims[1]; p <- dims[2]; K <- dims[3]; n_cov <- dims[4]
  ## version 2 records the test of the .assoc statistics; version 1 files do not
  test <- NA_character_
  if (version >= 2L) {
    code <- readBin(con, "integer", 1, size = 4, endian = "little")
    test <- c("wald", "score")[code + 1L]
    if (is.na(test)) stop("unknown test code ", code)
  }
  region   <- readBin(con, "character", 1)
  traits   <- readBin(con, "character", K)
  variants <- readBin(con, "character", p)
  pos      <- readBin(con, "integer", p, size = 4, endian = "little")
  af       <- readBin(con, "double", p, size = 8, endian = "little")
  list(version = version, n = n, p = p, K = K, n_cov = n_cov, test = test, region = region,
       traits = traits, variants = variants, pos = pos, af = af)
}

.qpgen_read_matrix <- function(con, nrow, ncol, rownames = NULL, colnames = NULL) {
  m <- matrix(readBin(con, "double", as.numeric(nrow) * as.numeric(ncol),
                      size = 8, endian = "little"), nrow = nrow, ncol = ncol)
  dimnames(m) <- list(rownames, colnames)
  m
}

## Sufficient statistics: list(XtX [p x p], Xty [p x K], yty [K], n, ...)
read_qpgen_suff <- function(path) {
  con <- gzfile(path, "rb"); on.exit(close(con))
  h <- .qpgen_read_header(con, "QPGNSUFF")
  h$XtX <- .qpgen_read_matrix(con, h$p, h$p, h$variants, h$variants)
  h$Xty <- .qpgen_read_matrix(con, h$p, h$K, h$variants, h$traits)
  h$yty <- setNames(readBin(con, "double", h$K, size = 8, endian = "little"), h$traits)
  h
}

## RSS summary statistics: list(R [p x p], z/bhat/shat [p x K], var_y [K], n, ...)
read_qpgen_rss <- function(path) {
  con <- gzfile(path, "rb"); on.exit(close(con))
  h <- .qpgen_read_header(con, "QPGN_RSS")
  h$R     <- .qpgen_read_matrix(con, h$p, h$p, h$variants, h$variants)
  h$z     <- .qpgen_read_matrix(con, h$p, h$K, h$variants, h$traits)
  h$bhat  <- .qpgen_read_matrix(con, h$p, h$K, h$variants, h$traits)
  h$shat  <- .qpgen_read_matrix(con, h$p, h$K, h$variants, h$traits)
  h$var_y <- setNames(readBin(con, "double", h$K, size = 8, endian = "little"), h$traits)
  h
}

## Fit susieR directly on the sufficient statistics of one trait:
##   susie_ss(XtX, Xty, yty, n)        in susieR >= 0.14 (exported as susie_ss)
##   susie_suff_stat(XtX, Xty, yty, n) in susieR <= 0.12 (CRAN)
## Extra arguments (L, coverage, min_abs_corr, estimate_residual_variance, ...)
## pass through. Coefficients are on the per-ALT-allele scale of the input.
susie_from_qpgen_suff <- function(ss, trait = 1, ...) {
  if (is.character(trait)) trait <- match(trait, ss$traits)
  if (is.na(trait)) stop("trait not found in file")
  Xty <- ss$Xty[, trait]; yty <- unname(ss$yty[trait])
  if (exists("susie_ss", envir = asNamespace("susieR")))
    susieR::susie_ss(XtX = ss$XtX, Xty = Xty, yty = yty, n = ss$n, ...)
  else
    susieR::susie_suff_stat(XtX = ss$XtX, Xty = Xty, yty = yty, n = ss$n, ...)
}

## Fit susieR on the RSS statistics of one trait via susie_rss(). Prefer the
## sufficient-statistics file and susie_from_qpgen_suff() when available.
## Supplying bhat/shat/var_y (rather than z alone) lets susie_rss rebuild X'X,
## X'y and y'y from R, n and var_y, but susie_rss() fixes the residual
## variance by default; pass estimate_residual_variance = TRUE to reproduce
## the sufficient-statistics fit (matches susie_ss() to ~1e-9 on in-sample R).
##
## susie_rss() rebuilds X'X, X'y and y'y from bhat, shat, var_y and n assuming
## shat is the Wald SE with df = n - 2. Unless rescale_shat = FALSE, shat in the
## file is first converted to that SE (n_cov = covariate degrees of freedom and
## test from the header, z = bhat / shat):
##   Wald (default):  shat * sqrt((n - 2 - n_cov) / (n - 2))
##   score (--score): shat * sqrt((n - 1 - n_cov - z^2) / (n - 2))
## Both are exact, except after --rint-after-adj or LOCO, where the Wald SE uses
## the residual sum of squares of the transformed phenotype rather than its y'y
## (relative error ~ the share of y'y the covariates explain, which region-assoc
## prints; typically < 1e-3). Version-1 files do not record the test; pass
## test = "score" for those written with --score (the default assumes Wald).
susie_from_qpgen_rss <- function(rss, trait = 1, rescale_shat = TRUE, test = NULL, ...) {
  if (is.character(trait)) trait <- match(trait, rss$traits)
  if (is.na(trait)) stop("trait not found in file")
  if (is.null(test)) test <- rss$test
  if (is.na(test)) {
    if (rescale_shat) warning("format version 1 does not record the test; assuming Wald. ",
                              "Pass test = \"score\" if the file was written with --score")
    test <- "wald"
  }
  test <- match.arg(test, c("wald", "score"))
  bhat <- rss$bhat[, trait]; shat <- rss$shat[, trait]
  if (rescale_shat) {
    if (test == "wald") {
      shat <- shat * sqrt((rss$n - 2 - rss$n_cov) / (rss$n - 2))
    } else {
      z2 <- (bhat / shat)^2  # score z^2 = r^2 (n - 1 - n_cov) <= n - 1 - n_cov
      shat <- shat * sqrt(pmax(rss$n - 1 - rss$n_cov - z2, 0) / (rss$n - 2))
    }
  }
  susieR::susie_rss(bhat = bhat, shat = shat, R = rss$R,
                    n = rss$n, var_y = unname(rss$var_y[trait]), ...)
}

## Dispatch on the file type: sufficient statistics -> susie_ss(); RSS -> susie_rss().
susie_from_qpgen <- function(obj, trait = 1, ...) {
  if (!is.null(obj$XtX)) susie_from_qpgen_suff(obj, trait, ...)
  else susie_from_qpgen_rss(obj, trait, ...)
}

## Derive the RSS quantities (bhat, shat, z, R, var_y) from sufficient
## statistics, using the same formulas as region-assoc's marginal test
## (centered X and y, no intercept, df = n - 2 - n_cov). Provided for convenience
## (e.g. to feed other RSS-based tools); not needed to run susieR.
suff_to_rss <- function(ss, trait = 1) {
  if (is.character(trait)) trait <- match(trait, ss$traits)
  d    <- diag(ss$XtX)
  Xty  <- ss$Xty[, trait]
  bhat <- Xty / d
  sse  <- pmax(ss$yty[trait] - Xty^2 / d, 0)
  shat <- sqrt(sse / (ss$n - 2 - ss$n_cov) / d)
  R    <- cov2cor(ss$XtX); R <- (R + t(R)) / 2   # bit-exact symmetry
  list(bhat = bhat, shat = shat, z = bhat / shat, R = R, n = ss$n,
       var_y = unname(ss$yty[trait]) / (ss$n - 1))
}
