# Study 4: comparison with measurement-error baselines. Runs everything and draws
# the figures.
#
# Baselines: corrected Lasso (Loh & Wainwright 2012) and GMUS (Rosenbaum &
# Tsybakov 2010; Sorensen et al. 2018), both from hdme; CoCoLasso (Datta & Zou
# 2017) from BDcocolasso. Each needs an input noise level, so it is run with
#   oracle  : the true delta_obs (tuned by each package's own cross-validation);
#             GMUS takes no noise level (lambda from cv.glmnet, delta by the elbow rule),
#   naive   : delta = 0 (told there is no measurement error),
#   half    : 0.5 x delta_obs (under-estimated),   [K = 1 design only]
#   double  : 2   x delta_obs (over-estimated).    [K = 1 design only]
# 2D-SS needs none of these.
#
# Designs: K = 1, K = 5, and a control with no confounder. Methods: Lasso, SS,
# Jitter-only, 2D-SS, plus the baselines with a misspecified noise level. Reports F1,
# TPR, FDR and exact recovery over 40 replicates, with paired differences.
#
#   Fig 1: F1 vs delta_obs (ours + baselines package-tuned and best-tuned ceiling), by design.
#   Fig 2: confounder inclusion vs delta_obs, by design with confounders.
#   ceiling : tuning-path point with the highest F1 at the true noise level (uses the
#             truth; an optimistic upper bound on what tuning could achieve).
#   Fig 3: sensitivity of corrected Lasso / CoCoLasso to the assumed noise level (F1).
#   Fig 4: min rule vs 1-SE rule vs ceiling for every baseline (F1), K = 1.
#
# Colours: Okabe-Ito (colour-blind safe); line types and shapes also differ.
# Set quick_test <- TRUE to time a 3-replicate run first.

suppressWarnings(options(warn = -1))

required_pkgs <- c("glmnet", "parallel", "ggplot2", "latex2exp")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs, dependencies = TRUE)
invisible(lapply(required_pkgs, library, character.only = TRUE))

if (!requireNamespace("hdme", quietly = TRUE)) install.packages("hdme")
library(hdme)

cocolasso_available <- TRUE
if (!requireNamespace("BDcocolasso", quietly = TRUE)) {
  ok <- tryCatch({ install.packages("BDcocolasso"); TRUE }, error = function(e) FALSE)
  if (!ok || !requireNamespace("BDcocolasso", quietly = TRUE)) {
    ok2 <- tryCatch({
      if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
      devtools::install_github("celiaescribe/BDcocolasso")
      TRUE
    }, error = function(e) FALSE)
    cocolasso_available <- ok2 && requireNamespace("BDcocolasso", quietly = TRUE)
  }
}
if (cocolasso_available) library(BDcocolasso)
cat("CoCoLasso available:", cocolasso_available, "\n")

RNGkind("L'Ecuyer-CMRG")
set.seed(2026)
n_cores <- max(1, min(7, parallel::detectCores() - 1))

quick_test <- FALSE

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
n <- 100
p <- 200                       # reduced for tractability of the baseline solvers
s <- 5
beta_S  <- c(-4, -2, 0.5, 0.6, 0.7)
a_conf  <- c(1/8, 1/4, 1/2, 1/2, 1/2)
idio_sd <- sqrt(1 - sum(a_conf^2))
true_support <- 1:s

K_grid         <- c(1, 5, 0)   # 0 = control without a confounder
delta_obs_grid <- c(0, 0.5, 1, 1.5, 2)
delta_grid     <- seq(0, 1, length.out = 5)
m <- length(delta_grid)
B <- 50
B_SS <- m * B

n1 <- floor(n / 2)
lambda_half <- sqrt(log(p) / n1)
lambda_full <- sqrt(log(p) / n)

cond_scale <- c(oracle = 1, naive = 0, half = 0.5, double = 2)

n_rep <- if (quick_test) 3 else 40

# Okabe-Ito
cols <- c("Lasso" = "#999999", "SS" = "#D55E00", "Jitter-only" = "#CC79A7",
          "2D-SS" = "#0072B2", "CorrectedLasso" = "#009E73", "GMUS" = "#E69F00",
          "CoCoLasso" = "#56B4E9")
shp  <- c("Lasso" = 16, "SS" = 18, "Jitter-only" = 4, "2D-SS" = 8,
          "CorrectedLasso" = 15, "GMUS" = 17, "CoCoLasso" = 3)
ours <- c("Lasso", "SS", "Jitter-only", "2D-SS")

# ---------------------------------------------------------------------------
# Data-generating process (K = 0 gives the control without a confounder)
# ---------------------------------------------------------------------------
generate_data <- function(K, delta_obs) {
  p_plain <- p - s - K
  X_S <- matrix(rnorm(n * s), n, s)
  parts <- list(X_S)
  if (K > 0) {
    base <- as.numeric(X_S %*% a_conf)
    parts[[2]] <- sapply(seq_len(K), function(k) base + idio_sd * rnorm(n))
  }
  parts[[length(parts) + 1]] <- matrix(rnorm(n * p_plain), n, p_plain)
  X0 <- do.call(cbind, parts)
  y <- as.numeric(X_S %*% beta_S + rnorm(n))
  X <- if (delta_obs > 0) X0 + matrix(rnorm(n * p, sd = delta_obs), n, p) else X0
  list(X = X, y = y)
}

# ---------------------------------------------------------------------------
# Our methods
# ---------------------------------------------------------------------------
fit_coef <- function(X, y, lam, alpha) {
  fit <- glmnet(X, y, alpha = alpha, lambda = lam, standardize = TRUE)
  as.numeric(coef(fit))[-1]
}

largest_gap_select <- function(fbar) {
  pp <- length(fbar)
  ord <- order(fbar, decreasing = TRUE)
  sf <- fbar[ord]
  khat <- which.max(sf[-pp] - sf[-1])
  which(fbar >= 0.5 * (sf[khat] + sf[khat + 1]))
}

aggregate_frequencies <- function(X, y, grid, n1_use, B_use, lam, alpha) {
  n_obs <- nrow(X); p_full <- ncol(X)
  freq_mat <- matrix(0, length(grid), p_full)
  for (i in seq_along(grid)) {
    delta <- grid[i]
    cnt <- numeric(p_full)
    for (b in seq_len(B_use)) {
      idx <- if (n1_use < n_obs) sample(n_obs, n1_use) else seq_len(n_obs)
      Xj <- X[idx, , drop = FALSE]
      if (delta > 0) Xj <- Xj + matrix(rnorm(n1_use * p_full, sd = delta), n1_use, p_full)
      cnt <- cnt + (fit_coef(scale(Xj), y[idx], lam, alpha) != 0)
    }
    freq_mat[i, ] <- cnt / B_use
  }
  colMeans(freq_mat)
}

# ---------------------------------------------------------------------------
# Baseline wrappers (following the package sources, the hdme vignette and the
# original papers)
#
# Corrected Lasso (Loh & Wainwright 2012; hdme):  needs Sigma_uu = sigma^2 I.
#   cv_corrected_lasso() picks the l1 radius by K-fold CV of the corrected loss
#   ($radius_min; $radius_1se is the 1-SE rule). Fit on the radius grid via
#   corrected_lasso(W, y, sigmaUU, radii)$betaCorr.
# CoCoLasso (Datta & Zou 2017; BDcocolasso::coco): needs tau = SD of the additive
#   error (Sigma_hat = Z'Z/n - tau^2 I). tau is on the scale of Z, so scale.Z = FALSE
#   (all columns have unit variance here). Tuning: calibrated K-fold CV over a
#   100-point lambda path: $beta.opt (CV minimum) and $beta.sd (1-SE rule).
# GMUS (Sorensen et al. 2018; hdme::gmus): does NOT take a noise level. It needs
#   (lambda, delta); the authors fix lambda at the cross-validated lasso value
#   (lambda_min or lambda_se, from cv.glmnet) and choose delta by the ELBOW RULE on
#   the number-of-nonzeros-versus-delta curve (Rosenbaum & Tsybakov 2010). We
#   automate the elbow as the point farthest below the chord of the normalised curve,
#   on delta in {0, 0.05, ..., 0.5} (the range of the authors' own elbow plots).
#
# Selections returned per baseline:
#   cv / cv1se : the packages' data-driven tuning (min rule / 1-SE rule),
#   best       : the tuning-path point with the highest F1. This uses the truth, so
#                it is an optimistic CEILING for any tuning rule.
# nonzero = |coef| > 1e-10 (hdme's threshold)
# ---------------------------------------------------------------------------
best_col <- function(B) {
  f1 <- apply(abs(B) > 1e-10, 2, function(v) selection_metrics(v, 0)[["F1"]])
  which.max(f1)
}

run_corrected_lasso <- function(X, y, sigma) {
  tryCatch({
    sigmaUU <- diag(sigma^2, p)
    cvfit <- cv_corrected_lasso(X, y, sigmaUU, n_folds = 5, no_radii = 20)
    radii <- cvfit$cv$radii
    B <- corrected_lasso(X, y, sigmaUU, radii = radii)$betaCorr      # p x 20
    pick <- function(r) abs(B[, which.min(abs(radii - r))]) > 1e-10
    list(cv = pick(cvfit$radius_min), cv1se = pick(cvfit$radius_1se),
         best = abs(B[, best_col(B)]) > 1e-10, err = NULL)
  }, error = function(e) list(cv = NULL, cv1se = NULL, best = NULL, err = conditionMessage(e)))
}

run_cocolasso <- function(X, y, sigma) {
  if (!cocolasso_available) return(list(cv = NULL, cv1se = NULL, best = NULL, err = "BDcocolasso not available"))
  tryCatch({
    invisible(utils::capture.output(
      fit <- BDcocolasso::coco(Z = X, y = matrix(y, ncol = 1), n = nrow(X), p = ncol(X),
                               center.Z = TRUE, scale.Z = FALSE, center.y = TRUE, scale.y = TRUE,
                               noise = "additive", block = FALSE, penalty = "lasso", tau = sigma)))
    B <- t(as.matrix(fit$data_beta[, -1, drop = FALSE]))              # p x path length
    list(cv = abs(as.numeric(fit$beta.opt)) > 1e-10,
         cv1se = abs(as.numeric(fit$beta.sd)) > 1e-10,
         best = abs(B[, best_col(B)]) > 1e-10, err = NULL)
  }, error = function(e) list(cv = NULL, cv1se = NULL, best = NULL, err = conditionMessage(e)))
}

# Automated elbow: point farthest below the chord of the normalised (delta, nonzeros) curve
elbow_index <- function(d, nz) {
  if (length(unique(nz)) == 1) return(1)
  x <- (d - min(d)) / (max(d) - min(d))
  y <- (nz - min(nz)) / (max(nz) - min(nz))
  which.max((y[1] + (y[length(y)] - y[1]) * x) - y)
}

run_gmus <- function(X, y) {
  tryCatch({
    d_grid <- seq(0, 0.5, by = 0.05)
    cvl <- glmnet::cv.glmnet(X, y)
    sel <- list(); allB <- NULL
    for (rule in c("min", "1se")) {
      lam <- if (rule == "min") cvl$lambda.min else cvl$lambda.1se
      B <- gmus(X, y, lambda = lam, delta = d_grid)$beta             # p x length(d_grid)
      B[is.na(B)] <- 0
      sel[[rule]] <- abs(B[, elbow_index(d_grid, colSums(abs(B) > 1e-10))]) > 1e-10
      allB <- cbind(allB, B)
    }
    list(cv = sel[["min"]], cv1se = sel[["1se"]],
         best = abs(allB[, best_col(allB)]) > 1e-10, err = NULL)
  }, error = function(e) list(cv = NULL, cv1se = NULL, best = NULL, err = conditionMessage(e)))
}

# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
selection_metrics <- function(sel_logical, K) {
  Shat <- which(sel_logical)
  tp <- length(intersect(Shat, true_support))
  fp <- length(setdiff(Shat, true_support))
  fn <- s - tp
  conf_idx <- if (K > 0) (s + 1):(s + K) else integer(0)
  c(TPR = tp / s,
    FDR = if (length(Shat) > 0) fp / length(Shat) else 0,
    F1 = if (tp > 0) 2 * tp / (2 * tp + fp + fn) else 0,
    exact = as.numeric(tp == s && fp == 0),
    conf_frac = if (K > 0) mean(conf_idx %in% Shat) else NA_real_,
    conf_any  = if (K > 0) as.numeric(any(conf_idx %in% Shat)) else NA_real_)
}

mk_row <- function(sel, K, delta_obs, rep_id, method, condition) {
  data.frame(t(selection_metrics(sel, K)), K = K, delta_obs = delta_obs, rep = rep_id,
             method = method, condition = condition, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# One replicate
# ---------------------------------------------------------------------------
run_one_replicate <- function(K, delta_obs, rep_id) {
  dat <- generate_data(K, delta_obs)
  X <- dat$X; y <- dat$y
  rows <- list(); errs <- list()
  
  rows[[1]] <- mk_row(fit_coef(scale(X), y, lambda_full, 1) != 0, K, delta_obs, rep_id, "Lasso", "-")
  f_ss <- aggregate_frequencies(X, y, 0, n1, B_SS, lambda_half, 1)
  rows[[2]] <- mk_row(seq_len(p) %in% largest_gap_select(f_ss), K, delta_obs, rep_id, "SS", "-")
  f_jo <- aggregate_frequencies(X, y, delta_grid, n, B, lambda_full, 1)
  rows[[3]] <- mk_row(seq_len(p) %in% largest_gap_select(f_jo), K, delta_obs, rep_id, "Jitter-only", "-")
  f_2d <- aggregate_frequencies(X, y, delta_grid, n1, B, lambda_half, 1)
  rows[[4]] <- mk_row(seq_len(p) %in% largest_gap_select(f_2d), K, delta_obs, rep_id, "2D-SS", "-")
  
  conds <- if (K == 1) names(cond_scale) else c("oracle", "naive")
  add_err <- function(bname, cond, msg)
    errs[[length(errs) + 1]] <<- data.frame(K = K, delta_obs = delta_obs, rep = rep_id,
                                            method = bname, condition = cond, error = msg, stringsAsFactors = FALSE)
  
  # baselines that take a noise level: corrected Lasso, CoCoLasso
  for (bname in c("CorrectedLasso", "CoCoLasso")) {
    runner <- if (bname == "CorrectedLasso") run_corrected_lasso else run_cocolasso
    cache <- list()                       # reuse fits with the same assumed noise level
    for (cond in conds) {
      sg <- cond_scale[[cond]] * delta_obs
      key <- as.character(sg)
      if (is.null(cache[[key]])) cache[[key]] <- runner(X, y, sg)
      res <- cache[[key]]
      if (is.null(res$cv)) { add_err(bname, cond, res$err); next }
      rows[[length(rows) + 1]] <- mk_row(res$cv, K, delta_obs, rep_id, bname, cond)
      if (cond == "oracle") {
        rows[[length(rows) + 1]] <- mk_row(res$cv1se, K, delta_obs, rep_id, bname, "oracle1se")
        rows[[length(rows) + 1]] <- mk_row(res$best,  K, delta_obs, rep_id, bname, "ceiling")
      }
    }
  }
  # GMUS: no noise input; lambda from cv.glmnet, delta by the elbow rule
  res <- run_gmus(X, y)
  if (is.null(res$cv)) {
    add_err("GMUS", "elbow", res$err)
  } else {
    rows[[length(rows) + 1]] <- mk_row(res$cv,    K, delta_obs, rep_id, "GMUS", "elbow")
    rows[[length(rows) + 1]] <- mk_row(res$cv1se, K, delta_obs, rep_id, "GMUS", "elbow1se")
    rows[[length(rows) + 1]] <- mk_row(res$best,  K, delta_obs, rep_id, "GMUS", "ceiling")
  }
  list(metrics = do.call(rbind, rows), errors = if (length(errs)) do.call(rbind, errs) else NULL)
}

tasks <- expand.grid(K = K_grid, delta_obs = delta_obs_grid, rep_id = seq_len(n_rep))
cat("Running", nrow(tasks), "replicates on", n_cores, "cores...\n")
res_list <- mclapply(seq_len(nrow(tasks)), function(i)
  run_one_replicate(tasks$K[i], tasks$delta_obs[i], tasks$rep_id[i]), mc.cores = n_cores)

ok <- !sapply(res_list, function(r) inherits(r, "try-error") || is.null(r$metrics))
if (any(!ok)) cat("Warning:", sum(!ok), "replicates failed and were dropped.\n")
results_df <- do.call(rbind, lapply(res_list[ok], function(r) r$metrics))
rownames(results_df) <- NULL
results_df$label <- ifelse(results_df$condition == "-", results_df$method,
                           paste0(results_df$method, " (", results_df$condition, ")"))
write.csv(results_df, "section3_study4_raw_results.csv", row.names = FALSE)

err_list <- lapply(res_list[ok], function(r) r$errors)
err_list <- err_list[!sapply(err_list, is.null)]
if (length(err_list) > 0) {
  errors_df <- do.call(rbind, err_list)
  write.csv(errors_df, "section3_study4_errors.csv", row.names = FALSE)
  cat("\n=== Baseline call failures (section3_study4_errors.csv) ===\n")
  print(unique(errors_df[, c("method", "condition", "error")]), row.names = FALSE)
} else cat("\nNo baseline call failures recorded.\n")

# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------
se <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA else sd(x) / sqrt(length(x)) }
mn <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA else mean(x) }
metrics <- c("TPR", "FDR", "F1", "exact", "conf_frac", "conf_any")
grp <- split(results_df, list(results_df$K, results_df$delta_obs, results_df$label), drop = TRUE)
summary_df <- do.call(rbind, lapply(grp, function(d) {
  out <- data.frame(K = d$K[1], delta_obs = d$delta_obs[1], method = d$method[1],
                    condition = d$condition[1], label = d$label[1], n = nrow(d))
  for (mm in metrics) { out[[paste0(mm, "_mean")]] <- mn(d[[mm]]); out[[paste0(mm, "_se")]] <- se(d[[mm]]) }
  out
}))
rownames(summary_df) <- NULL
write.csv(summary_df, "section3_study4_summary.csv", row.names = FALSE)

# paired F1 difference: 2D-SS minus each other method, same replicate
ref <- results_df[results_df$label == "2D-SS", c("K", "delta_obs", "rep", "F1")]
names(ref)[4] <- "F1_2DSS"
mg <- merge(results_df[results_df$label != "2D-SS", ], ref, by = c("K", "delta_obs", "rep"))
mg$diff <- mg$F1_2DSS - mg$F1
paired_df <- do.call(rbind, lapply(split(mg, list(mg$K, mg$delta_obs, mg$label), drop = TRUE), function(d)
  data.frame(K = d$K[1], delta_obs = d$delta_obs[1], label = d$label[1],
             diff_mean = mean(d$diff), diff_se = sd(d$diff) / sqrt(nrow(d)), n = nrow(d))))
rownames(paired_df) <- NULL
write.csv(paired_df, "section3_study4_paired.csv", row.names = FALSE)
print(subset(summary_df, label %in% c("Lasso", "SS", "Jitter-only", "2D-SS"))[, c("K", "delta_obs", "label", "F1_mean", "exact_mean", "conf_any_mean")], digits = 3, row.names = FALSE)

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
lab_K <- function(x) factor(ifelse(x == 0, "No confounder (control)",
                                   ifelse(x == 1, "1 confounder", paste0(x, " confounders"))),
                            levels = c("1 confounder", "5 confounders", "No confounder (control)"))
th <- theme_bw(base_size = 13) + theme(legend.position = "bottom", strip.text = element_text(size = 10))

summary_df$K_f <- lab_K(summary_df$K)
summary_df$method <- factor(summary_df$method, levels = names(cols))

# Fig 1 / 2: ours + baselines (package-tuned, and the best-tuned ceiling)
d <- subset(summary_df, condition %in% c("-", "oracle", "elbow", "ceiling"))
d$cond_f <- factor(ifelse(d$condition == "-", "no noise input needed",
                          ifelse(d$condition %in% c("oracle", "elbow"), "baseline: package-tuned (CV / elbow)",
                                 "baseline: best-tuned (ceiling)")),
                   levels = c("no noise input needed", "baseline: package-tuned (CV / elbow)",
                              "baseline: best-tuned (ceiling)"))
lt <- c("no noise input needed" = "solid", "baseline: package-tuned (CV / elbow)" = "dashed",
        "baseline: best-tuned (ceiling)" = "dotted")
d$lw <- ifelse(d$method == "2D-SS", 1.3, 0.7)

f1 <- ggplot(d, aes(delta_obs, F1_mean, color = method, shape = method, linetype = cond_f,
                    group = label)) +
  geom_line(aes(linewidth = lw)) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(F1_mean - F1_se, 0), ymax = pmin(F1_mean + F1_se, 1)),
                width = 0, alpha = 0.4, linetype = "solid") +
  facet_wrap(~ K_f) +
  scale_color_manual(values = cols) + scale_shape_manual(values = shp) +
  scale_linetype_manual(values = lt) + scale_linewidth_identity() +
  labs(x = TeX("$\\delta_{obs}$"), y = "F1 (support recovery)", color = NULL, shape = NULL, linetype = NULL) +
  th + guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2), linetype = guide_legend(nrow = 3))
ggsave("plot_study4_f1.pdf", f1, width = 12, height = 5.5)

dc <- subset(d, K > 0)
f2 <- ggplot(dc, aes(delta_obs, conf_any_mean, color = method, shape = method, linetype = cond_f,
                     group = label)) +
  geom_line(aes(linewidth = lw)) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(conf_any_mean - conf_any_se, 0), ymax = pmin(conf_any_mean + conf_any_se, 1)),
                width = 0, alpha = 0.4, linetype = "solid") +
  facet_wrap(~ K_f) +
  scale_color_manual(values = cols) + scale_shape_manual(values = shp) +
  scale_linetype_manual(values = lt) + scale_linewidth_identity() +
  labs(x = TeX("$\\delta_{obs}$"), y = "P(any confounder included)", color = NULL, shape = NULL, linetype = NULL) +
  th + guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2), linetype = guide_legend(nrow = 3))
ggsave("plot_study4_confounder.pdf", f2, width = 10, height = 5.5)

# Fig 3: sensitivity to the assumed noise level (baselines that take one), K = 1
d3 <- subset(summary_df, method %in% c("CorrectedLasso", "CoCoLasso") & K == 1 &
               condition %in% c("oracle", "half", "double", "naive"))
d3$cond_f <- factor(d3$condition, levels = c("oracle", "half", "double", "naive"),
                    labels = c("true level", "0.5 x true", "2 x true", "zero"))
ref3 <- subset(summary_df, method == "2D-SS" & K == 1)
ref3 <- do.call(rbind, lapply(c("CorrectedLasso", "CoCoLasso"), function(b) transform(ref3, panel = b)))
d3$panel <- as.character(d3$method)
f3 <- ggplot(d3, aes(delta_obs, F1_mean)) +
  geom_line(aes(linetype = cond_f, color = cond_f, group = cond_f)) +
  geom_point(aes(color = cond_f, shape = cond_f), size = 2) +
  geom_line(data = ref3, aes(delta_obs, F1_mean), color = "#0072B2", linewidth = 1.3, inherit.aes = FALSE) +
  facet_wrap(~ panel) +
  scale_color_manual(values = c("true level" = "#000000", "0.5 x true" = "#E69F00",
                                "2 x true" = "#CC79A7", "zero" = "#999999")) +
  scale_linetype_manual(values = c("true level" = "solid", "0.5 x true" = "dashed",
                                   "2 x true" = "dotdash", "zero" = "dotted")) +
  scale_shape_manual(values = c("true level" = 16, "0.5 x true" = 17, "2 x true" = 15, "zero" = 4)) +
  labs(x = TeX("$\\delta_{obs}$"), y = "F1 (thick blue line: 2D-SS)",
       color = "Assumed noise level", linetype = "Assumed noise level", shape = "Assumed noise level") + th
ggsave("plot_study4_sensitivity.pdf", f3, width = 9, height = 5)

# Fig 4: tuning rule (min vs 1-SE vs ceiling) for every baseline, K = 1
d4 <- subset(summary_df, method %in% c("CorrectedLasso", "GMUS", "CoCoLasso") & K == 1 &
               condition %in% c("oracle", "elbow", "oracle1se", "elbow1se", "ceiling"))
d4$rule <- factor(ifelse(d4$condition %in% c("oracle", "elbow"), "package-tuned, min rule",
                         ifelse(d4$condition %in% c("oracle1se", "elbow1se"), "package-tuned, 1-SE rule",
                                "best-tuned (ceiling)")),
                  levels = c("package-tuned, min rule", "package-tuned, 1-SE rule", "best-tuned (ceiling)"))
d4$panel <- as.character(d4$method)
ref4 <- do.call(rbind, lapply(c("CorrectedLasso", "GMUS", "CoCoLasso"), function(b) transform(subset(summary_df, method == "2D-SS" & K == 1), panel = b)))
f4 <- ggplot(d4, aes(delta_obs, F1_mean, color = rule, linetype = rule, shape = rule)) +
  geom_line() + geom_point(size = 2) +
  geom_line(data = ref4, aes(delta_obs, F1_mean), color = "#0072B2", linewidth = 1.3, inherit.aes = FALSE) +
  facet_wrap(~ panel) +
  scale_color_manual(values = c("package-tuned, min rule" = "#E69F00", "package-tuned, 1-SE rule" = "#CC79A7",
                                "best-tuned (ceiling)" = "#000000")) +
  scale_linetype_manual(values = c("package-tuned, min rule" = "dashed", "package-tuned, 1-SE rule" = "dotdash",
                                   "best-tuned (ceiling)" = "dotted")) +
  scale_shape_manual(values = c("package-tuned, min rule" = 17, "package-tuned, 1-SE rule" = 15,
                                "best-tuned (ceiling)" = 16)) +
  labs(x = TeX("$\\delta_{obs}$"), y = "F1 (thick blue line: 2D-SS)", color = NULL, linetype = NULL, shape = NULL) + th
ggsave("plot_study4_tuning_rules.pdf", f4, width = 11, height = 5)

cat("Saved: section3_study4_{raw_results,summary,paired,errors}.csv and plot_study4_{f1,confounder,sensitivity,tuning_rules}.pdf\n")
