# Study 6: elastic-net-path baseline ("ENet-SS") against SS and 2D-SS.
#
# Question (reviewer objection): by Proposition 1 a jittered Lasso fit is, in
# expectation, a Lasso with an added ridge term delta^2. Is 2D-SS then just
# Stability Selection run along a deterministic ridge path?
#
# Three procedures, same data, same half-sample size, same grid, same B per
# grid point, same lambda, same largest-gap rule:
#   SS       : delta = 0 only, B_SS = m * B half-samples (as in Study 2)
#   ENet-SS  : deterministic elastic net over the grid delta_j; NO random noise.
#              The ridge weight is delta_j^2, matched to the expected jittered fit
#              (see enet_pars below).
#   2D-SS    : Gaussian design jittering over the same grid (as in Study 2).
# SS vs ENet-SS   = effect of the ridge term.
# ENet-SS vs 2D-SS = effect of the random noise beyond its expectation.
#
# Design: same as Study 2 (Example 2 embedded in p = 1000).
# Parts: A  = largest-gap selection metrics over (n, K), clean training data
#        A2 = complementary-pairs estimate of the confounder inclusion
#             probability p_k(delta)
#        B  = prediction MSE on fresh data at future noise levels delta_test
#             (training clean), after an OLS refit on the selected set; Oracle added
#        C  = selection metrics when the TRAINING data already carry measurement
#             error, delta_obs in {0.5, 1}
# Each part has a switch (RUN_A, RUN_A2, RUN_B, RUN_C). A part that is switched off
# is read back from its CSV if the plots need it, so finished parts need not be rerun.
#
# Output: section3_enetss_{selection_raw,selection_summary,cp_summary,
#         mse_raw,mse_summary,obs_raw,obs_summary}.csv and plot_enetss_*.pdf
# Set QUICK <- TRUE for a fast smoke test (a few minutes).

suppressWarnings(options(warn = -1))

required_pkgs <- c("glmnet", "parallel", "ggplot2")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs, dependencies = TRUE)
invisible(lapply(required_pkgs, library, character.only = TRUE))

RNGkind("L'Ecuyer-CMRG")
set.seed(2026)
n_cores <- max(1, min(7, parallel::detectCores() - 1))

QUICK  <- FALSE
RUN_A  <- TRUE
RUN_A2 <- TRUE
RUN_B  <- TRUE
RUN_C  <- TRUE

# ---------------------------------------------------------------------------
# Settings (identical to Study 2)
# ---------------------------------------------------------------------------
p <- 1000
s <- 5
beta_S  <- c(-4, -2, 0.5, 0.6, 0.7)
a_conf  <- c(1/8, 1/4, 1/2, 1/2, 1/2)
idio_sd <- sqrt(1 - sum(a_conf^2))

delta_obs_train <- 0
delta_grid <- seq(0, 1, length.out = 5)
m <- length(delta_grid)
B <- 50

n_grid <- c(100, 200, 400)
K_grid <- c(1, 5)
n_rep_select <- 200
N_pairs      <- 200

n_mse           <- 100                 # Part B
delta_test_grid <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)
n_rep_mse       <- 200
n_test          <- 500

n_grid_obs     <- c(200, 400)          # Part C
delta_obs_grid <- c(0.5, 1)
n_rep_obs      <- 100

if (QUICK) { n_grid <- 100; K_grid <- 1; n_rep_select <- 10; N_pairs <- 10; B <- 10
n_rep_mse <- 5; n_test <- 100; n_grid_obs <- 200; delta_obs_grid <- 0.5; n_rep_obs <- 5 }

methods <- c("SS", "ENet-SS", "2D-SS")
method_colors <- c("SS" = "#D55E00", "ENet-SS" = "#009E73", "2D-SS" = "#0072B2", "Oracle" = "#000000")
method_shapes <- c("SS" = 18, "ENet-SS" = 17, "2D-SS" = 8, "Oracle" = 1)
method_lty    <- c("SS" = "solid", "ENet-SS" = "dashed", "2D-SS" = "solid", "Oracle" = "longdash")

make_cfg <- function(n, K) {
  n1 <- floor(n / 2)
  list(n = n, K = K, n1 = n1, p_plain = p - s - K, conf_idx = (s + 1):(s + K),
       lam_half = sqrt(log(p) / n1), B_SS = m * B)
}

# ---------------------------------------------------------------------------
# Data-generating process and helpers (copied from Study 2)
# ---------------------------------------------------------------------------
generate_data <- function(cfg, n_obs, delta_obs) {
  X_S <- matrix(rnorm(n_obs * s), n_obs, s)
  base <- as.numeric(X_S %*% a_conf)
  X_C <- sapply(seq_len(cfg$K), function(k) base + idio_sd * rnorm(n_obs))
  X_plain <- matrix(rnorm(n_obs * cfg$p_plain), n_obs, cfg$p_plain)
  X0 <- cbind(X_S, X_C, X_plain)
  y <- as.numeric(X_S %*% beta_S + rnorm(n_obs))
  X <- if (delta_obs > 0) X0 + matrix(rnorm(n_obs * p, sd = delta_obs), n_obs, p) else X0
  list(X = X, y = y)
}

fit_beta <- function(X, y, lam, alpha)
  as.numeric(coef(glmnet(X, y, alpha = alpha, lambda = lam, standardize = TRUE)))[-1]

largest_gap_select <- function(fbar) {
  pp <- length(fbar)
  ord <- order(fbar, decreasing = TRUE)
  sf <- fbar[ord]
  khat <- which.max(sf[-pp] - sf[-1])
  which(fbar >= 0.5 * (sf[khat] + sf[khat + 1]))
}

# ---------------------------------------------------------------------------
# ENet counterpart of a jittered Lasso fit at noise level delta.
# The jittered design is standardised after jitter: Z = (X + W) / sqrt(1 + delta^2),
# so with b' = b / sqrt(1 + delta^2) (coefficients on the un-jittered standardised X),
#   E[ (1/2n)||y - Z b||^2 ] + lam ||b||_1
#     = (1/2n)||y - X b'||^2 + (delta^2/2) ||b'||^2 + lam sqrt(1 + delta^2) ||b'||_1 .
# glmnet (gaussian) standardises y internally, which leaves the L1 part on the
# original scale but divides the ridge part by sd(y) (checked numerically), i.e. it
# minimises (1/2n)RSS + L [ (1-a)/(2 sd_y) ||b||^2 + a ||b||_1 ], hence
#   L (1 - a) = delta^2 * sd_y,   L a = lam sqrt(1 + delta^2),
# with sd_y the (1/n) standard deviation of the y used in the fit.
# At delta = 0 this is the plain Lasso (alpha = 1, L = lam), i.e. SS.
# ---------------------------------------------------------------------------
enet_pars <- function(lam, delta, y) {
  sd_y <- sqrt(mean((y - mean(y))^2))
  l1 <- lam * sqrt(1 + delta^2)
  L  <- l1 + delta^2 * sd_y
  list(lambda = L, alpha = l1 / L)
}

# Selection frequencies averaged over the grid. type = "jitter" (random noise,
# Lasso) or "enet" (no noise, deterministic elastic net with ridge delta^2).
aggregate_frequencies <- function(X, y, grid, n1_use, B_use, lam, type = c("jitter", "enet")) {
  type <- match.arg(type)
  n_obs <- nrow(X); p_full <- ncol(X)
  freq_mat <- matrix(0, length(grid), p_full)
  for (i in seq_along(grid)) {
    delta <- grid[i]
    cnt <- numeric(p_full)
    for (b in seq_len(B_use)) {
      idx <- if (n1_use < n_obs) sample(n_obs, n1_use) else seq_len(n_obs)
      Xj <- X[idx, , drop = FALSE]
      if (type == "jitter") {
        if (delta > 0) Xj <- Xj + matrix(rnorm(n1_use * p_full, sd = delta), n1_use, p_full)
        bh <- fit_beta(scale(Xj), y[idx], lam, 1)
      } else {
        pr <- enet_pars(lam, delta, y[idx])
        bh <- fit_beta(scale(Xj), y[idx], pr$lambda, pr$alpha)
      }
      cnt <- cnt + (bh != 0)
    }
    freq_mat[i, ] <- cnt / B_use
  }
  colMeans(freq_mat)
}

summarise_col <- function(x) {
  x <- x[!is.na(x) & is.finite(x)]
  c(mean = mean(x), se = sd(x) / sqrt(length(x)))
}

# ===========================================================================
# Part A: selection metrics at training time
# ===========================================================================
sel_metrics <- function(Shat, cfg) {
  c(conf_rate = mean(cfg$conf_idx %in% Shat),
    any_conf  = as.numeric(any(cfg$conf_idx %in% Shat)),
    tp        = sum(Shat %in% 1:s),
    fp        = sum(!(Shat %in% 1:s)),
    exact     = as.numeric(all(1:s %in% Shat) && all(Shat %in% 1:s)))
}

run_one_selection <- function(rep_id, cfg, dobs = delta_obs_train) {
  dat <- generate_data(cfg, cfg$n, dobs)
  X <- dat$X; y <- dat$y
  sets <- list(
    `SS`      = largest_gap_select(aggregate_frequencies(X, y, 0,          cfg$n1, cfg$B_SS, cfg$lam_half, "jitter")),
    `ENet-SS` = largest_gap_select(aggregate_frequencies(X, y, delta_grid, cfg$n1, B,        cfg$lam_half, "enet")),
    `2D-SS`   = largest_gap_select(aggregate_frequencies(X, y, delta_grid, cfg$n1, B,        cfg$lam_half, "jitter")))
  do.call(rbind, lapply(methods, function(mm)
    data.frame(rep = rep_id, n = cfg$n, K = cfg$K, delta_obs = dobs, method = mm,
               t(sel_metrics(sets[[mm]], cfg)))))
}

cells <- expand.grid(n = n_grid, K = K_grid)
vars <- c("conf_rate", "any_conf", "tp", "fp", "exact")
summarise_sel <- function(df, by) {
  out <- do.call(rbind, lapply(split(df, df[by], drop = TRUE), function(d) {
    r <- sapply(vars, function(v) summarise_col(d[[v]]))
    o <- d[1, by, drop = FALSE]
    for (v in vars) { o[[v]] <- r["mean", v]; o[[paste0(v, "_se")]] <- r["se", v] }
    o
  }))
  rownames(out) <- NULL
  out$method <- factor(out$method, levels = methods)
  out
}

if (RUN_A) {
  select_df <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    cfg <- make_cfg(cells$n[i], cells$K[i])
    cat(sprintf("Part A: n = %d, K = %d\n", cfg$n, cfg$K))
    do.call(rbind, mclapply(seq_len(n_rep_select), run_one_selection, cfg = cfg, mc.cores = n_cores))
  }))
  rownames(select_df) <- NULL
  write.csv(select_df, "section3_enetss_selection_raw.csv", row.names = FALSE)
  sel_summary <- summarise_sel(select_df, c("n", "K", "method"))
  sel_summary <- sel_summary[order(sel_summary$K, sel_summary$n, sel_summary$method), ]
  write.csv(sel_summary, "section3_enetss_selection_summary.csv", row.names = FALSE)
  print(sel_summary, digits = 3, row.names = FALSE)
  cat("Part A done.\n\n")
} else {
  sel_summary <- read.csv("section3_enetss_selection_summary.csv")
  sel_summary$method <- factor(sel_summary$method, levels = methods)
}

# ===========================================================================
# Part A2: complementary-pairs estimate of p_confounder(delta)
# (fresh independent dataset per pair, as in Study 2)
# ===========================================================================
cp_one <- function(cfg, grid, lam, type) {
  Reduce(c, lapply(seq_len(N_pairs), function(b) {
    dat <- generate_data(cfg, cfg$n, delta_obs_train)
    perm <- sample(cfg$n)
    halves <- list(perm[seq_len(cfg$n1)], perm[(cfg$n1 + 1):(2 * cfg$n1)])
    mean(sapply(halves, function(idxA) {
      Xb <- dat$X[idxA, , drop = FALSE]; yb <- dat$y[idxA]
      mean(sapply(grid, function(delta) {
        if (type == "jitter") {
          Wb <- if (delta > 0) matrix(rnorm(cfg$n1 * p, sd = delta), cfg$n1, p) else 0
          bh <- fit_beta(scale(Xb + Wb), yb, lam, 1)
        } else {
          pr <- enet_pars(lam, delta, yb)
          bh <- fit_beta(scale(Xb), yb, pr$lambda, pr$alpha)
        }
        mean(bh[cfg$conf_idx] != 0)
      }))
    }))
  }))
}

if (RUN_A2) {
  cp_summary <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    cfg <- make_cfg(cells$n[i], cells$K[i])
    cat(sprintf("Part A2: n = %d, K = %d\n", cfg$n, cfg$K))
    res <- mclapply(1:3, function(j) switch(j,
                                            cp_one(cfg, 0,          cfg$lam_half, "jitter"),
                                            cp_one(cfg, delta_grid, cfg$lam_half, "enet"),
                                            cp_one(cfg, delta_grid, cfg$lam_half, "jitter")),
                    mc.cores = min(3, n_cores))
    do.call(rbind, lapply(1:3, function(j)
      data.frame(n = cfg$n, K = cfg$K, method = methods[j],
                 rate = mean(res[[j]]), se = sd(res[[j]]) / sqrt(N_pairs))))
  }))
  write.csv(cp_summary, "section3_enetss_cp_summary.csv", row.names = FALSE)
  print(cp_summary, digits = 3, row.names = FALSE)
  cat("Part A2 done.\n\n")
} else if (file.exists("section3_enetss_cp_summary.csv")) {
  cp_summary <- read.csv("section3_enetss_cp_summary.csv")
}

# ===========================================================================
# Part B: prediction MSE on fresh data at future noise levels delta_test.
# Training data clean (delta_obs = 0). The selected set is refitted by OLS, as in
# Study 2 Part B. Oracle = OLS on the true support.
# ===========================================================================
ols_refit <- function(X_train, y_train, Shat) {
  beta_full <- rep(0, ncol(X_train))
  if (length(Shat) == 0 || length(Shat) >= nrow(X_train) - 1)
    return(list(beta = beta_full, intercept = mean(y_train)))
  cf <- coef(lm(y_train ~ X_train[, Shat, drop = FALSE]))
  cf[is.na(cf)] <- 0
  beta_full[Shat] <- cf[-1]
  list(beta = beta_full, intercept = cf[1])
}
predict_mse <- function(beta_full, intercept, X_test, y_test)
  mean((y_test - intercept - as.numeric(X_test %*% beta_full))^2)

run_one_mse <- function(rep_id, cfg) {
  train <- generate_data(cfg, cfg$n, 0)
  X <- train$X; y <- train$y
  models <- list(
    `SS`      = ols_refit(X, y, largest_gap_select(aggregate_frequencies(X, y, 0,          cfg$n1, cfg$B_SS, cfg$lam_half, "jitter"))),
    `ENet-SS` = ols_refit(X, y, largest_gap_select(aggregate_frequencies(X, y, delta_grid, cfg$n1, B,        cfg$lam_half, "enet"))),
    `2D-SS`   = ols_refit(X, y, largest_gap_select(aggregate_frequencies(X, y, delta_grid, cfg$n1, B,        cfg$lam_half, "jitter"))),
    `Oracle`  = ols_refit(X, y, 1:s))
  do.call(rbind, lapply(delta_test_grid, function(dt) {
    test <- generate_data(cfg, n_test, dt)
    do.call(rbind, lapply(names(models), function(nm)
      data.frame(rep = rep_id, K = cfg$K, delta_test = dt, method = nm,
                 mse = predict_mse(models[[nm]]$beta, models[[nm]]$intercept, test$X, test$y))))
  }))
}

if (RUN_B) {
  mse_df <- do.call(rbind, lapply(K_grid, function(K) {
    cfg <- make_cfg(n_mse, K)
    cat(sprintf("Part B: n = %d, K = %d\n", cfg$n, cfg$K))
    do.call(rbind, mclapply(seq_len(n_rep_mse), run_one_mse, cfg = cfg, mc.cores = n_cores))
  }))
  rownames(mse_df) <- NULL
  write.csv(mse_df, "section3_enetss_mse_raw.csv", row.names = FALSE)
  mse_summary <- do.call(rbind, lapply(split(mse_df, list(mse_df$K, mse_df$delta_test, mse_df$method), drop = TRUE),
                                       function(d) { r <- summarise_col(d$mse)
                                       data.frame(K = d$K[1], delta_test = d$delta_test[1], method = d$method[1],
                                                  mse_mean = r["mean"], mse_se = r["se"]) }))
  rownames(mse_summary) <- NULL
  write.csv(mse_summary, "section3_enetss_mse_summary.csv", row.names = FALSE)
  print(mse_summary, digits = 4, row.names = FALSE)
  cat("Part B done.\n\n")
}

# ===========================================================================
# Part C: measurement error already present in the training data.
# X_obs = X + N(0, delta_obs^2) on every feature (y is generated from the clean X_S).
# Same selectors, same grid, same lambda; selection metrics as in Part A.
# ===========================================================================
cells_obs <- expand.grid(n = n_grid_obs, K = K_grid, delta_obs = delta_obs_grid)
if (RUN_C) {
  obs_df <- do.call(rbind, lapply(seq_len(nrow(cells_obs)), function(i) {
    cfg <- make_cfg(cells_obs$n[i], cells_obs$K[i])
    cat(sprintf("Part C: n = %d, K = %d, delta_obs = %.1f\n", cfg$n, cfg$K, cells_obs$delta_obs[i]))
    do.call(rbind, mclapply(seq_len(n_rep_obs), run_one_selection, cfg = cfg,
                            dobs = cells_obs$delta_obs[i], mc.cores = n_cores))
  }))
  rownames(obs_df) <- NULL
  write.csv(obs_df, "section3_enetss_obs_raw.csv", row.names = FALSE)
  obs_summary <- summarise_sel(obs_df, c("n", "K", "delta_obs", "method"))
  obs_summary <- obs_summary[order(obs_summary$delta_obs, obs_summary$K, obs_summary$n, obs_summary$method), ]
  write.csv(obs_summary, "section3_enetss_obs_summary.csv", row.names = FALSE)
  print(obs_summary, digits = 3, row.names = FALSE)
  cat("Part C done.\n\n")
}

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
th <- theme_bw(base_size = 13) + theme(legend.position = "bottom", strip.text = element_text(size = 10))
lab_n <- function(x) factor(paste0("n = ", x), levels = paste0("n = ", n_grid))
lab_K <- function(x) factor(ifelse(x == 1, "1 confounder", paste0(x, " confounders")),
                            levels = ifelse(K_grid == 1, "1 confounder", paste0(K_grid, " confounders")))

d1 <- sel_summary; d1$n_f <- lab_n(d1$n); d1$K_f <- lab_K(d1$K)
p1 <- ggplot(d1, aes(method, conf_rate, fill = method)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = pmax(conf_rate - conf_rate_se, 0), ymax = pmin(conf_rate + conf_rate_se, 1)), width = 0.2) +
  geom_text(aes(label = sprintf("%.2f", conf_rate)), vjust = -0.6, size = 3) +
  facet_grid(K_f ~ n_f) + scale_fill_manual(values = method_colors) +
  scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "Mean fraction of confounders included") + th + theme(legend.position = "none")
ggsave("plot_enetss_confounder_inclusion.pdf", p1, width = 9, height = 5.5)

d3 <- sel_summary
d3l <- rbind(
  data.frame(d3[, c("n", "K", "method")], metric = "True positives (of 5)", value = d3$tp, se = d3$tp_se),
  data.frame(d3[, c("n", "K", "method")], metric = "False positives (count)", value = d3$fp, se = d3$fp_se),
  data.frame(d3[, c("n", "K", "method")], metric = "Exact support recovery", value = d3$exact, se = d3$exact_se))
d3l$metric <- factor(d3l$metric, levels = c("True positives (of 5)", "False positives (count)", "Exact support recovery"))
d3l$K_f <- lab_K(d3l$K)
p3 <- ggplot(d3l, aes(n, value, color = method, shape = method, linetype = method)) +
  geom_line() + geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(value - se, 0), ymax = value + se), width = 0, alpha = 0.5) +
  facet_grid(metric ~ K_f, scales = "free_y") + scale_x_log10(breaks = n_grid) +
  scale_color_manual(values = method_colors) + scale_shape_manual(values = method_shapes) +
  scale_linetype_manual(values = method_lty) +
  labs(x = "n (log scale)", y = NULL, color = NULL, shape = NULL, linetype = NULL) + th
ggsave("plot_enetss_selection_metrics.pdf", p3, width = 8, height = 8)

if (RUN_A2) {
  d2 <- cp_summary; d2$method <- factor(d2$method, levels = methods)
  d2$n_f <- lab_n(d2$n); d2$K_f <- lab_K(d2$K)
  p2 <- ggplot(d2, aes(method, rate, fill = method)) +
    geom_col(width = 0.6) +
    geom_errorbar(aes(ymin = pmax(rate - se, 0), ymax = rate + se), width = 0.15) +
    geom_text(aes(label = sprintf("%.3f", rate)), vjust = -0.6, size = 3) +
    facet_grid(K_f ~ n_f) + scale_fill_manual(values = method_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = "Complementary-pairs estimate of p_k(delta)") + th + theme(legend.position = "none")
  ggsave("plot_enetss_confounder_inclusion_cp.pdf", p2, width = 9, height = 5.5)
}

if (RUN_B || file.exists("section3_enetss_mse_summary.csv")) {
  if (!RUN_B) mse_summary <- read.csv("section3_enetss_mse_summary.csv")
  d4 <- mse_summary; d4$method <- factor(d4$method, levels = c(methods, "Oracle")); d4$K_f <- lab_K(d4$K)
  p4 <- ggplot(d4, aes(delta_test, mse_mean, color = method, shape = method, linetype = method, fill = method)) +
    geom_ribbon(aes(ymin = pmax(mse_mean - mse_se, 1e-6), ymax = mse_mean + mse_se), alpha = 0.12, color = NA) +
    geom_line() + geom_point(size = 2) +
    facet_wrap(~ K_f) + scale_y_log10() +
    scale_color_manual(values = method_colors) + scale_fill_manual(values = method_colors) +
    scale_shape_manual(values = method_shapes) + scale_linetype_manual(values = method_lty) +
    labs(x = "delta_test (unknown at training time)", y = "Prediction MSE on fresh data (log scale)",
         color = NULL, fill = NULL, shape = NULL, linetype = NULL) + th
  ggsave("plot_enetss_surprise_curve.pdf", p4, width = 11, height = 5)
}

if (RUN_C || file.exists("section3_enetss_obs_summary.csv")) {
  if (!RUN_C) { obs_summary <- read.csv("section3_enetss_obs_summary.csv")
  obs_summary$method <- factor(obs_summary$method, levels = methods) }
  d5 <- obs_summary
  d5$n_f <- factor(paste0("n = ", d5$n)); d5$K_f <- lab_K(d5$K)
  d5$dobs_f <- factor(paste0("delta_obs = ", d5$delta_obs))
  p5 <- ggplot(d5, aes(method, conf_rate, fill = method)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = pmax(conf_rate - conf_rate_se, 0), ymax = pmin(conf_rate + conf_rate_se, 1)), width = 0.2) +
    geom_text(aes(label = sprintf("%.2f", conf_rate)), vjust = -0.6, size = 3) +
    facet_grid(K_f ~ dobs_f + n_f) + scale_fill_manual(values = method_colors) +
    scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25)) +
    labs(x = NULL, y = "Mean fraction of confounders included") + th + theme(legend.position = "none")
  ggsave("plot_enetss_obs_confounder_inclusion.pdf", p5, width = 11, height = 5.5)
  p6 <- ggplot(d5, aes(method, exact, fill = method)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = pmax(exact - exact_se, 0), ymax = pmin(exact + exact_se, 1)), width = 0.2) +
    geom_text(aes(label = sprintf("%.2f", exact)), vjust = -0.6, size = 3) +
    facet_grid(K_f ~ dobs_f + n_f) + scale_fill_manual(values = method_colors) +
    scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25)) +
    labs(x = NULL, y = "Exact support recovery") + th + theme(legend.position = "none")
  ggsave("plot_enetss_obs_exact_recovery.pdf", p6, width = 11, height = 5.5)
}

cat("Saved section3_enetss_*.csv and plot_enetss_*.pdf\n")