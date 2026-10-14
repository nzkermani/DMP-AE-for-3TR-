## Checks for the internal sPLS-DA engine in R/lomics_spls.R
##   1. keepX = all variables reproduces pls::plsr(method = "oscorespls")
##      predictions on the scaled class-indicator matrix (PLS2 / NIPALS).
##   2. keepX selects exactly keepX variables per component.
##   3. With planted signal in 5 of 200 variables, sPLS-DA (keepX = 5)
##      recovers them.
## Run from the project root:  Rscript tests/test_spls_internal.R
source("R/lomics_spls.R")
set.seed(1)
n <- 40; p <- 30
X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("v", 1:p)))
y <- factor(rep(c("A", "B"), each = n / 2))
X[y == "B", 1:3] <- X[y == "B", 1:3] + 1

ok <- TRUE
f <- .lomics_spls_internal(X, y, ncomp = 3)
P1 <- lomics_spls_predict(f, X)
if (requireNamespace("pls", quietly = TRUE)) {
  Y <- stats::model.matrix(~ y - 1)
  Ys <- scale(Y); Xs <- scale(X)
  m <- pls::plsr(Ys ~ Xs, ncomp = 3, method = "oscorespls", scale = FALSE)
  P2 <- drop(stats::predict(m, ncomp = 3))
  P2 <- sweep(sweep(P2, 2, attr(Ys, "scaled:scale"), "*"), 2, attr(Ys, "scaled:center"), "+")
  d <- max(abs(P1 - P2))
  message(sprintf("1. max |internal - pls::plsr| prediction difference: %.2e", d))
  ok <- ok && d < 1e-6
} else message("1. pls not installed - skipped")

f5 <- .lomics_spls_internal(X, y, ncomp = 2, keepX = c(4, 6))
nz <- colSums(f5$loadings != 0)
message("2. selected per component: ", paste(nz, collapse = ", "))
ok <- ok && all(nz == c(4, 6))

p <- 200
X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("v", 1:p)))
X[y == "B", 1:5] <- X[y == "B", 1:5] + 1.2
sel <- lomics_spls_selected(.lomics_spls_internal(X, y, ncomp = 1, keepX = 5))$feature
message("3. selected: ", paste(sel, collapse = ", "))
ok <- ok && mean(sel %in% paste0("v", 1:5)) >= 0.8

message(if (ok) "PASS" else "FAIL")
if (!ok) quit(status = 1)
