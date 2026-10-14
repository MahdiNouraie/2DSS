# Affymetrix Rat Genome 230 2.0 (GSE5680), response = TRIM32 probe.
# Part A: full-data selection (as in the paper).
# Part B: repeated train/test splits: own rule and matched-size (top-k) comparison.
# Outputs: real_data_summary.csv, real_data_selected_probes.csv,
#          real_split_raw.csv, real_split_summary.csv, real_split_probe_freq.csv,
#          plot_real_split.pdf

options(warn = -1)
for (pk in c("GEOquery", "glmnet")) if (!requireNamespace(pk)) install.packages(pk)
library(GEOquery); library(glmnet); library(parallel)

## ---- settings ----
N_SPLITS  <- 50      # number of random train/test splits
TRAIN_FRAC <- 2/3
B_FULL    <- 100     # resamples per level, Part A
B_SPLIT   <- 50      # resamples per level, Part B (speed)
N_CORES   <- max(1, detectCores() - 1)   # set 1 on Windows
delta_grid <- seq(0.05, 2.5, length.out = 10)

set.seed(26)

## ---- data ----
gset <- getGEO("GSE5680", GSEMatrix = TRUE, getGPL = FALSE)
if (length(gset) > 1) idx <- grep("GPL1355", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
ex <- exprs(gset)
response_var <- ex["1389163_at", ]
ex <- ex[rownames(ex) != "1389163_at", ]
ex <- data.frame(t(ex)); rm(gset, idx); gc()
max_expr <- as.numeric(apply(ex, 2, max))
filtered_ex <- ex[, max_expr > quantile(max_expr, 0.25)]
range_expr <- as.numeric(apply(filtered_ex, 2, function(x) diff(range(x))))
filtered_ex <- filtered_ex[, range_expr >= 2]
x <- scale(as.matrix(filtered_ex), center = TRUE, scale = TRUE)   # unsupervised
y <- as.numeric(scale(response_var, center = TRUE, scale = FALSE))
rm(ex, filtered_ex, max_expr, range_expr, response_var); gc()
n <- nrow(x); p <- ncol(x)
cat("n =", n, " p =", p, "\n")

## ---- methods ----
aggregate_frequencies <- function(X, y, grid, n1_use, B_use, lam) {
  n_obs <- nrow(X); p_full <- ncol(X)
  freq_mat <- matrix(0, length(grid), p_full)
  for (i in seq_along(grid)) {
    delta <- grid[i]
    cnt <- numeric(p_full)
    for (b in seq_len(B_use)) {
      ib <- sample.int(n_obs, n1_use)
      Wb <- if (delta > 0) matrix(rnorm(n1_use * p_full, sd = delta), n1_use, p_full) else 0
      fit <- glmnet(X[ib, , drop = FALSE] + Wb, y[ib], alpha = 1, lambda = lam, standardize = FALSE)
      cnt <- cnt + (as.numeric(coef(fit)[-1]) != 0)
    }
    freq_mat[i, ] <- cnt / B_use
  }
  colMeans(freq_mat)      # equal weights
}
largest_gap_select <- function(fbar) {
  ord <- order(fbar, decreasing = TRUE); sf <- fbar[ord]
  gaps <- sf[-length(sf)] - sf[-1]
  k <- which.max(gaps); tau <- 0.5 * (sf[k] + sf[k + 1])
  list(Shat = which(fbar >= tau), tau_hat = tau)
}
run_methods <- function(X, y, B_use, lam) {
  n_obs <- nrow(X); n1 <- floor(n_obs / 2)
  fit <- glmnet(X, y, alpha = 1, lambda = lam, standardize = FALSE)
  S_lasso <- which(as.numeric(coef(fit)[-1]) != 0)
  f_ss  <- aggregate_frequencies(X, y, 0, n1, B_use, lam)
  f_jit <- aggregate_frequencies(X, y, delta_grid, n_obs, B_use, lam)
  f_2d  <- aggregate_frequencies(X, y, delta_grid, n1, B_use, lam)
  list(sets = list(Lasso = S_lasso,
                   SS = largest_gap_select(f_ss)$Shat,
                   SS_tau06 = which(f_ss >= 0.6),
                   Jitter = largest_gap_select(f_jit)$Shat,
                   TwoDSS = largest_gap_select(f_2d)$Shat),
       tau = c(SS = largest_gap_select(f_ss)$tau_hat,
               Jitter = largest_gap_select(f_jit)$tau_hat,
               TwoDSS = largest_gap_select(f_2d)$tau_hat),
       lasso_fit = fit, freqs = list(SS = f_ss, Jitter = f_jit, TwoDSS = f_2d))
}

## =========================== Part A: full data ===========================
foldid <- sample(rep(1:10, length.out = n))
lam_full <- cv.glmnet(x, y, alpha = 1, standardize = FALSE, foldid = foldid)$lambda.1se
rA <- run_methods(x, y, B_FULL, lam_full)
nm <- colnames(x)
sel <- lapply(rA$sets, function(s) nm[s])
cat("\nFull data: sizes\n"); print(sapply(sel, length))
summary_df <- data.frame(
  method = c("Lasso", "SS (largest-gap)", "SS (tau=0.6)", "Jitter-only", "2D-SS"),
  n_selected = sapply(sel, length),
  tau_hat = c(NA, rA$tau["SS"], 0.6, rA$tau["Jitter"], rA$tau["TwoDSS"]))
write.csv(summary_df, "real_data_summary.csv", row.names = FALSE)
pr <- data.frame(probe = unique(unlist(sel)))
for (k in names(sel)) pr[[k]] <- pr$probe %in% sel[[k]]
write.csv(pr, "real_data_selected_probes.csv", row.names = FALSE)

## ====================== Part B: repeated splits ==========================
## (i) each method with its own rule (largest-gap; Lasso at CV lambda), and
## (ii) matched size: top-k probes, k in KS, ranked by selection frequency
##      (SS, Jitter-only, 2D-SS) or by order of entry on the Lasso path.
KS <- c(3, 5, 10)
meths <- c("Lasso", "SS", "Jitter", "TwoDSS")
rank_ties <- function(f) order(-f, runif(length(f)))     # random tie-breaking
lasso_topk <- function(X, y, k) {
  fp <- glmnet(X, y, alpha = 1, standardize = FALSE, nlambda = 100)
  cf <- as.matrix(fp$beta) != 0
  entry <- apply(cf, 1, function(r) if (any(r)) which(r)[1] else Inf)
  last <- abs(as.numeric(fp$beta[, ncol(fp$beta)]))
  o <- order(entry, -last, runif(length(entry)))
  o[seq_len(k)]
}
refit_mse <- function(S, Xtr, ytr, Xte, yte) {
  if (length(S) == 0) return(mean((yte - mean(ytr))^2))
  f <- lm(y ~ ., data = data.frame(y = ytr, Xtr[, S, drop = FALSE]))
  pr <- as.numeric(predict(f, data.frame(Xte[, S, drop = FALSE])))
  pr[is.na(pr)] <- mean(ytr)
  mean((yte - pr)^2)
}
one_split <- function(r) {
  set.seed(1000 + r)
  tr <- sample.int(n, floor(TRAIN_FRAC * n)); te <- setdiff(seq_len(n), tr)
  Xtr <- x[tr, ]; ytr <- y[tr]; Xte <- x[te, ]; yte <- y[te]
  lam <- cv.glmnet(Xtr, ytr, alpha = 1, standardize = FALSE, nfolds = 10)$lambda.1se
  R <- run_methods(Xtr, ytr, B_SPLIT, lam)
  base <- mean((yte - mean(ytr))^2)
  rows <- list()
  add <- function(m, rule, k, S) {
    rows[[length(rows) + 1]] <<- data.frame(split = r, method = m, rule = rule, k = k,
                                            size = length(S), relmse = refit_mse(S, Xtr, ytr, Xte, yte) / base,
                                            sel = I(paste(S, collapse = ",")))
  }
  for (m in meths) add(m, "gap", NA, R$sets[[m]])
  for (k in KS) {
    add("Lasso", "topk", k, lasso_topk(Xtr, ytr, k))
    for (m in c("SS", "Jitter", "TwoDSS")) add(m, "topk", k, rank_ties(R$freqs[[m]])[seq_len(k)])
  }
  do.call(rbind, rows)
}
res <- mclapply(seq_len(N_SPLITS), one_split, mc.cores = N_CORES)
raw <- do.call(rbind, res)
write.csv(raw, "real_split_raw.csv", row.names = FALSE)

nog <- function(Z) {
  M <- nrow(Z); pf <- colMeans(Z); kbar <- mean(rowSums(Z))
  if (kbar == 0 || kbar == ncol(Z)) return(NA)
  1 - (M / (M - 1)) * mean(pf * (1 - pf)) / ((kbar / ncol(Z)) * (1 - kbar / ncol(Z)))
}
jac <- function(Z) {
  M <- nrow(Z); v <- c()
  for (i in 1:(M - 1)) for (j in (i + 1):M) {
    u <- sum(Z[i, ] | Z[j, ]); v <- c(v, if (u == 0) NA else sum(Z[i, ] & Z[j, ]) / u)
  }
  mean(v, na.rm = TRUE)
}
se <- function(v) sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v)))
cfgs <- unique(raw[, c("method", "rule", "k")])
summ <- do.call(rbind, lapply(seq_len(nrow(cfgs)), function(i) {
  m <- cfgs$method[i]; ru <- cfgs$rule[i]; k <- cfgs$k[i]
  d <- raw[raw$method == m & raw$rule == ru & (if (is.na(k)) is.na(raw$k) else raw$k %in% k), ]
  Z <- matrix(0L, nrow(d), p)
  for (j in seq_len(nrow(d))) { s <- as.integer(strsplit(d$sel[j], ",")[[1]]); if (length(s)) Z[j, s] <- 1L }
  data.frame(method = m, rule = ru, k = k, size_mean = mean(d$size),
             relmse_mean = mean(d$relmse), relmse_se = se(d$relmse),
             relmse_median = median(d$relmse), jaccard = jac(Z), nogueira = nog(Z))
}))
write.csv(summ, "real_split_summary.csv", row.names = FALSE)

## paired differences (relative test MSE): other minus 2D-SS, per configuration
pd <- do.call(rbind, lapply(c("gap", paste0("topk", KS)), function(cf) {
  ru <- if (cf == "gap") "gap" else "topk"; k <- if (cf == "gap") NA else as.numeric(sub("topk", "", cf))
  d <- raw[raw$rule == ru & (if (is.na(k)) is.na(raw$k) else raw$k %in% k), ]
  w <- reshape(d[, c("split", "method", "relmse")], idvar = "split", timevar = "method", direction = "wide")
  do.call(rbind, lapply(c("Lasso", "SS", "Jitter"), function(m) {
    dd <- w[[paste0("relmse.", m)]] - w$relmse.TwoDSS
    data.frame(config = cf, comparison = paste(m, "minus 2D-SS"), mean = mean(dd), se = sd(dd) / sqrt(length(dd)))
  }))
}))
write.csv(pd, "real_split_paired.csv", row.names = FALSE)

## probe frequency across splits (top-5 rule)
fr <- do.call(rbind, lapply(meths, function(m) {
  d <- raw[raw$method == m & raw$rule == "topk" & raw$k %in% 5, ]
  tab <- table(unlist(strsplit(d$sel, ",")))
  data.frame(method = m, probe = nm[as.integer(names(tab))], freq = as.numeric(tab) / nrow(d))
}))
write.csv(fr, "real_split_probe_freq.csv", row.names = FALSE)
print(summ, digits = 3); print(pd, digits = 3)

## ---- plot (Okabe-Ito) ----
cols <- c(Lasso = "#999999", SS = "#E69F00", Jitter = "#009E73", TwoDSS = "#0072B2")
lab  <- c(Lasso = "Lasso", SS = "SS", Jitter = "Jitter-only", TwoDSS = "2D-SS")
pch  <- c(Lasso = 15, SS = 17, Jitter = 16, TwoDSS = 18)
pdf("plot_real_split.pdf", width = 11, height = 3.8)
par(mfrow = c(1, 3), mar = c(4, 4.2, 2.5, 1))
tk <- summ[summ$rule == "topk", ]
mk <- function(col, ylab, main, se_col = NULL) {
  yl <- range(c(tk[[col]] - if (!is.null(se_col)) tk[[se_col]] else 0,
                tk[[col]] + if (!is.null(se_col)) tk[[se_col]] else 0), na.rm = TRUE)
  plot(NA, xlim = range(KS), ylim = yl, xlab = "Model size k", ylab = ylab, main = main, las = 1, xaxt = "n")
  axis(1, KS)
  for (m in meths) {
    d <- tk[tk$method == m, ]; d <- d[order(d$k), ]
    lines(d$k, d[[col]], col = cols[m], lwd = 2); points(d$k, d[[col]], col = cols[m], pch = pch[m], cex = 1.2)
    if (!is.null(se_col)) arrows(d$k, d[[col]] - d[[se_col]], d$k, d[[col]] + d[[se_col]],
                                 angle = 90, code = 3, length = .03, col = cols[m])
  }
}
mk("relmse_mean", "Test MSE / null MSE", "Prediction at matched size", "relmse_se")
legend("topright", lab[meths], col = cols[meths], pch = pch[meths], lwd = 2, bty = "n", cex = .8)
mk("nogueira", "Nogueira stability", "Stability at matched size")
g <- summ[summ$rule == "gap", ]; g <- g[match(meths, g$method), ]
barplot(g$relmse_mean, col = cols[meths], names.arg = lab[meths], las = 1, ylim = c(0, 1),
        ylab = "Test MSE / null MSE", main = "Own rule (largest-gap)")
text(seq(0.7, by = 1.2, length.out = 4), g$relmse_mean, paste0("size ", round(g$size_mean, 1)), pos = 3, cex = .8)
dev.off()
cat("\nDone.\n")