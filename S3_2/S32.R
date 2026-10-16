# Study 2: robustness of the aggregated procedure. Runs everything and draws the figures.
#
# Design: Example 2 (Sigma_SS = I_5, beta_S = (-4,-2,.5,.6,.7), confounder
# covariance a = (1/8,1/4,1/2,1/2,1/2)), embedded in p = 1000 features, trained at
# delta_obs = 0. Varies n in {100,200,400} and K in {1,5}; reports confounder
# inclusion, true and false positives, exact recovery, and MSE against an oracle.
#
#   Part A  : confounder inclusion and selection metrics at training time
#             (largest-gap rule for SS / Jitter-only / 2D-SS; non-zero coefficients
#             for Lasso / Ridge / ENet), all six methods, over (n, K).
#   Part A2 : complementary-pairs estimate of the confounder inclusion probability
#             p_k(delta) for SS and 2D-SS (the quantity of Theorem 8), over (n, K).
#   Part B  : prediction MSE on fresh data at future noise levels delta_test
#             ("surprise curve"), n = 100, K in {1,5}, with an oracle baseline.
#   Part B2 : the same MSE for SS split by whether SS included a confounder.
#
# Output: section3_study2_*.csv and plot_study2_*.pdf

suppressWarnings(options(warn = -1))

required_pkgs <- c("glmnet", "parallel", "ggplot2", "latex2exp")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs, dependencies = TRUE)
invisible(lapply(required_pkgs, library, character.only = TRUE))

RNGkind("L'Ecuyer-CMRG")
set.seed(2026)
n_cores <- max(1, min(7, parallel::detectCores() - 1))

# ---------------------------------------------------------------------------
# Settings
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

n_grid <- c(100, 200, 400)     # Parts A, A2
K_grid <- c(1, 5)              # number of confounders
n_rep_select <- 200
N_pairs      <- 200

n_mse           <- 100         # Part B
delta_test_grid <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)
n_rep_mse       <- 200
n_test          <- 500

methods <- c("Lasso", "Ridge", "ENet", "SS", "Jitter-only", "2D-SS")

# Okabe-Ito palette
method_colors <- c("Lasso" = "#56B4E9", "Ridge" = "#999999", "ENet" = "#009E73",
                   "SS" = "#D55E00", "Jitter-only" = "#CC79A7", "2D-SS" = "#0072B2",
                   "Oracle" = "#000000")
method_shapes <- c("Lasso" = 16, "Ridge" = 15, "ENet" = 17, "SS" = 18,
                   "Jitter-only" = 4, "2D-SS" = 8, "Oracle" = 1)
method_lty <- c("Lasso" = "solid", "Ridge" = "dotted", "ENet" = "dashed", "SS" = "solid",
                "Jitter-only" = "dotdash", "2D-SS" = "solid", "Oracle" = "longdash")

make_cfg <- function(n, K) {
  n1 <- floor(n / 2)
  list(n = n, K = K, n1 = n1, p_plain = p - s - K, conf_idx = (s + 1):(s + K),
       lam_half = sqrt(log(p) / n1), lam_full = sqrt(log(p) / n), B_SS = m * B)
}

# ---------------------------------------------------------------------------
# Data-generating process and selectors
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

fit_coef <- function(X, y, lam, alpha) {
  fit <- glmnet(X, y, alpha = alpha, lambda = lam, standardize = TRUE)
  list(beta = as.numeric(coef(fit))[-1], intercept = as.numeric(coef(fit))[1],
       center = attr(X, "scaled:center"), scale = attr(X, "scaled:scale"))
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
      cnt <- cnt + (fit_coef(scale(Xj), y[idx], lam, alpha)$beta != 0)
    }
    freq_mat[i, ] <- cnt / B_use
  }
  colMeans(freq_mat)    # uniform weights over the grid
}

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

summarise_col <- function(x) {
  x <- x[!is.na(x) & is.finite(x)]
  c(mean = mean(x), se = sd(x) / sqrt(length(x)))
}

# ===========================================================================
# Part A: confounder inclusion and selection metrics at training time
# ===========================================================================
sel_metrics <- function(Shat, cfg) {
  c(conf_rate = mean(cfg$conf_idx %in% Shat),
    any_conf  = as.numeric(any(cfg$conf_idx %in% Shat)),
    tp        = sum(Shat %in% 1:s),
    fp        = sum(!(Shat %in% 1:s)),
    exact     = as.numeric(all(1:s %in% Shat) && all(Shat %in% 1:s)))
}

run_one_selection <- function(rep_id, cfg) {
  dat <- generate_data(cfg, cfg$n, delta_obs_train)
  X <- dat$X; y <- dat$y
  Xs <- scale(X)
  S_lasso <- which(fit_coef(Xs, y, cfg$lam_full, 1)$beta   != 0)
  S_ridge <- which(fit_coef(Xs, y, cfg$lam_full, 0)$beta   != 0)
  S_enet  <- which(fit_coef(Xs, y, cfg$lam_full, 0.5)$beta != 0)
  S_SS  <- largest_gap_select(aggregate_frequencies(X, y, 0, cfg$n1, cfg$B_SS, cfg$lam_half, 1))
  S_JIT <- largest_gap_select(aggregate_frequencies(X, y, delta_grid, cfg$n, B, cfg$lam_full, 1))
  S_2D  <- largest_gap_select(aggregate_frequencies(X, y, delta_grid, cfg$n1, B, cfg$lam_half, 1))
  sets <- list(Lasso = S_lasso, Ridge = S_ridge, ENet = S_enet,
               SS = S_SS, `Jitter-only` = S_JIT, `2D-SS` = S_2D)
  do.call(rbind, lapply(methods, function(mm)
    data.frame(rep = rep_id, n = cfg$n, K = cfg$K, method = mm, t(sel_metrics(sets[[mm]], cfg)))))
}

cells <- expand.grid(n = n_grid, K = K_grid)
select_df <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  cfg <- make_cfg(cells$n[i], cells$K[i])
  cat(sprintf("Part A: n = %d, K = %d\n", cfg$n, cfg$K))
  do.call(rbind, mclapply(seq_len(n_rep_select), run_one_selection, cfg = cfg, mc.cores = n_cores))
}))
rownames(select_df) <- NULL
write.csv(select_df, "section3_study2_selection_raw.csv", row.names = FALSE)

sel_summary <- do.call(rbind, lapply(split(select_df, list(select_df$n, select_df$K, select_df$method), drop = TRUE),
                                     function(d) {
                                       r <- sapply(c("conf_rate", "any_conf", "tp", "fp", "exact"), function(v) summarise_col(d[[v]]))
                                       data.frame(n = d$n[1], K = d$K[1], method = d$method[1],
                                                  conf_rate = r["mean", "conf_rate"], conf_rate_se = r["se", "conf_rate"],
                                                  any_conf = r["mean", "any_conf"],   any_conf_se = r["se", "any_conf"],
                                                  tp = r["mean", "tp"],               tp_se = r["se", "tp"],
                                                  fp = r["mean", "fp"],               fp_se = r["se", "fp"],
                                                  exact = r["mean", "exact"],         exact_se = r["se", "exact"])
                                     }))
rownames(sel_summary) <- NULL
write.csv(sel_summary, "section3_study2_selection_summary.csv", row.names = FALSE)
print(sel_summary, digits = 3, row.names = FALSE)
cat("Part A done.\n\n")

# ===========================================================================
# Part A2: complementary-pairs estimate of p_confounder(delta)
# (each pair uses a fresh, independent dataset: p_k(delta) is marginal over data,
# half-sample and noise)
# ===========================================================================
cp_one <- function(cfg, grid, lam) {
  Reduce(c, lapply(seq_len(N_pairs), function(b) {
    dat <- generate_data(cfg, cfg$n, delta_obs_train)
    perm <- sample(cfg$n)
    halves <- list(perm[seq_len(cfg$n1)], perm[(cfg$n1 + 1):(2 * cfg$n1)])
    mean(sapply(halves, function(idxA) {
      Xb <- dat$X[idxA, , drop = FALSE]; yb <- dat$y[idxA]
      mean(sapply(grid, function(delta) {
        Wb <- if (delta > 0) matrix(rnorm(cfg$n1 * p, sd = delta), cfg$n1, p) else 0
        mean(fit_coef(scale(Xb + Wb), yb, lam, 1)$beta[cfg$conf_idx] != 0)
      }))
    }))
  }))
}

cp_summary <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  cfg <- make_cfg(cells$n[i], cells$K[i])
  cat(sprintf("Part A2: n = %d, K = %d\n", cfg$n, cfg$K))
  res <- mclapply(c("SS", "2D-SS"), function(mm)
    if (mm == "SS") cp_one(cfg, 0, cfg$lam_half) else cp_one(cfg, delta_grid, cfg$lam_half),
    mc.cores = min(2, n_cores))
  do.call(rbind, lapply(1:2, function(j)
    data.frame(n = cfg$n, K = cfg$K, method = c("SS", "2D-SS")[j],
               rate = mean(res[[j]]), se = sd(res[[j]]) / sqrt(N_pairs))))
}))
write.csv(cp_summary, "section3_study2_cp_summary.csv", row.names = FALSE)
print(cp_summary, digits = 3, row.names = FALSE)
cat("Part A2 done.\n\n")

# ===========================================================================
# Part B / B2: prediction MSE at future noise levels
# ===========================================================================
run_one_mse <- function(rep_id, cfg) {
  train <- generate_data(cfg, cfg$n, delta_obs_train)
  X_train <- train$X; y_train <- train$y
  Xs <- scale(X_train)
  rb <- function(f) f$beta / f$scale
  ri <- function(f) f$intercept - sum(f$beta * f$center / f$scale)
  models <- list()
  for (nm in c("Lasso", "Ridge", "ENet")) {
    a <- c(Lasso = 1, Ridge = 0, ENet = 0.5)[[nm]]
    f <- fit_coef(Xs, y_train, cfg$lam_full, a)
    models[[nm]] <- list(beta = rb(f), intercept = ri(f))
  }
  S_SS  <- largest_gap_select(aggregate_frequencies(X_train, y_train, 0, cfg$n1, cfg$B_SS, cfg$lam_half, 1))
  S_JIT <- largest_gap_select(aggregate_frequencies(X_train, y_train, delta_grid, cfg$n, B, cfg$lam_full, 1))
  S_2D  <- largest_gap_select(aggregate_frequencies(X_train, y_train, delta_grid, cfg$n1, B, cfg$lam_half, 1))
  models[["SS"]]          <- ols_refit(X_train, y_train, S_SS)
  models[["Jitter-only"]] <- ols_refit(X_train, y_train, S_JIT)
  models[["2D-SS"]]       <- ols_refit(X_train, y_train, S_2D)
  models[["Oracle"]]      <- ols_refit(X_train, y_train, 1:s)
  ss_conf <- any(cfg$conf_idx %in% S_SS)
  do.call(rbind, lapply(delta_test_grid, function(dt) {
    test <- generate_data(cfg, n_test, dt)
    do.call(rbind, lapply(names(models), function(nm)
      data.frame(rep = rep_id, K = cfg$K, delta_test = dt, method = nm,
                 mse = predict_mse(models[[nm]]$beta, models[[nm]]$intercept, test$X, test$y),
                 ss_conf_included = ss_conf)))
  }))
}

mse_df <- do.call(rbind, lapply(K_grid, function(K) {
  cfg <- make_cfg(n_mse, K)
  cat(sprintf("Part B: n = %d, K = %d\n", cfg$n, cfg$K))
  do.call(rbind, mclapply(seq_len(n_rep_mse), run_one_mse, cfg = cfg, mc.cores = n_cores))
}))
rownames(mse_df) <- NULL
write.csv(mse_df, "section3_study2_mse_raw.csv", row.names = FALSE)

mse_summary <- do.call(rbind, lapply(split(mse_df, list(mse_df$K, mse_df$delta_test, mse_df$method), drop = TRUE),
                                     function(d) { r <- summarise_col(d$mse)
                                     data.frame(K = d$K[1], delta_test = d$delta_test[1], method = d$method[1],
                                                mse_mean = r["mean"], mse_se = r["se"]) }))
rownames(mse_summary) <- NULL
write.csv(mse_summary, "section3_study2_mse_summary.csv", row.names = FALSE)
print(mse_summary, digits = 4, row.names = FALSE)

# Part B2: SS split by whether it included a confounder; 2D-SS as reference
ss_rows <- mse_df[mse_df$method == "SS", ]
ss_rows$group <- ifelse(ss_rows$ss_conf_included, "SS (confounder included)", "SS (confounder excluded)")
td_rows <- mse_df[mse_df$method == "2D-SS", ]; td_rows$group <- "2D-SS"
cond_df <- rbind(ss_rows[, c("K", "delta_test", "group", "mse")], td_rows[, c("K", "delta_test", "group", "mse")])
cond_summary <- do.call(rbind, lapply(split(cond_df, list(cond_df$K, cond_df$delta_test, cond_df$group), drop = TRUE),
                                      function(d) { r <- summarise_col(d$mse)
                                      data.frame(K = d$K[1], delta_test = d$delta_test[1], group = d$group[1],
                                                 mse_mean = r["mean"], mse_se = r["se"], n = nrow(d)) }))
rownames(cond_summary) <- NULL
write.csv(cond_summary, "section3_study2_mse_conditioned_summary.csv", row.names = FALSE)
print(cond_summary, digits = 4, row.names = FALSE)
cat("Part B done.\n\n")

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
th <- theme_bw(base_size = 13) +
  theme(legend.position = "bottom", strip.text = element_text(size = 10))
lab_n <- function(x) factor(paste0("n = ", x), levels = paste0("n = ", n_grid))
lab_K <- function(x) factor(ifelse(x == 1, "1 confounder", paste0(x, " confounders")),
                            levels = ifelse(K_grid == 1, "1 confounder", paste0(K_grid, " confounders")))

# Fig 1: confounder inclusion (largest-gap / non-zero), all six methods
d1 <- sel_summary; d1$method <- factor(d1$method, levels = methods)
d1$n_f <- lab_n(d1$n); d1$K_f <- lab_K(d1$K)
p1 <- ggplot(d1, aes(method, conf_rate, fill = method)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = pmax(conf_rate - conf_rate_se, 0), ymax = pmin(conf_rate + conf_rate_se, 1)), width = 0.2) +
  geom_text(aes(label = sprintf("%.2f", conf_rate)), vjust = -0.6, size = 2.8) +
  facet_grid(K_f ~ n_f) +
  scale_fill_manual(values = method_colors) + scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "Mean fraction of confounders included") + th +
  theme(legend.position = "none", axis.text.x = element_text(angle = 40, hjust = 1))
ggsave("plot_study2_confounder_inclusion.pdf", p1, width = 11, height = 6)

# Fig 2: complementary-pairs estimate for SS vs 2D-SS
d2 <- cp_summary; d2$method <- factor(d2$method, levels = c("SS", "2D-SS"))
d2$n_f <- lab_n(d2$n); d2$K_f <- lab_K(d2$K)
p2 <- ggplot(d2, aes(method, rate, fill = method)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = pmax(rate - se, 0), ymax = rate + se), width = 0.15) +
  geom_text(aes(label = sprintf("%.3f", rate)), vjust = -0.6, size = 3) +
  facet_grid(K_f ~ n_f) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = TeX("Complementary-pairs estimate of $p_k(\\delta)$")) + th +
  theme(legend.position = "none")
ggsave("plot_study2_confounder_inclusion_cp.pdf", p2, width = 9, height = 5.5)

# Fig 3: selection metrics (Ridge omitted: it selects every feature)
d3 <- sel_summary[sel_summary$method != "Ridge", ]
d3l <- rbind(
  data.frame(d3[, c("n", "K", "method")], metric = "True positives (of 5)", value = d3$tp, se = d3$tp_se),
  data.frame(d3[, c("n", "K", "method")], metric = "False positives (count)", value = d3$fp, se = d3$fp_se),
  data.frame(d3[, c("n", "K", "method")], metric = "Exact support recovery", value = d3$exact, se = d3$exact_se))
d3l$method <- factor(d3l$method, levels = setdiff(methods, "Ridge"))
d3l$metric <- factor(d3l$metric, levels = c("True positives (of 5)", "False positives (count)", "Exact support recovery"))
d3l$K_f <- lab_K(d3l$K)
p3 <- ggplot(d3l, aes(n, value, color = method, shape = method, linetype = method)) +
  geom_line() + geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(value - se, 0), ymax = value + se), width = 0, alpha = 0.5) +
  facet_grid(metric ~ K_f, scales = "free_y") + scale_x_log10(breaks = n_grid) +
  scale_color_manual(values = method_colors) + scale_shape_manual(values = method_shapes) +
  scale_linetype_manual(values = method_lty) +
  labs(x = "n (log scale)", y = NULL, color = NULL, shape = NULL, linetype = NULL) + th
ggsave("plot_study2_selection_metrics.pdf", p3, width = 9, height = 8)

# Fig 4: surprise curve
d4 <- mse_summary; d4$method <- factor(d4$method, levels = c(methods, "Oracle")); d4$K_f <- lab_K(d4$K)
p4 <- ggplot(d4, aes(delta_test, mse_mean, color = method, shape = method, linetype = method, fill = method)) +
  geom_ribbon(aes(ymin = pmax(mse_mean - mse_se, 1e-6), ymax = mse_mean + mse_se), alpha = 0.12, color = NA) +
  geom_line() + geom_point(size = 2) +
  facet_wrap(~ K_f) + scale_y_log10() +
  scale_color_manual(values = method_colors) + scale_fill_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) + scale_linetype_manual(values = method_lty) +
  labs(x = TeX("$\\delta_{test}$ (unknown at training time)"), y = "Prediction MSE on fresh data (log scale)",
       color = NULL, fill = NULL, shape = NULL, linetype = NULL) + th
ggsave("plot_study2_surprise_curve.pdf", p4, width = 11, height = 5)

# Fig 5: MSE of SS conditioned on confounder inclusion
grp_levels <- c("SS (confounder excluded)", "SS (confounder included)", "2D-SS")
grp_cols <- c("SS (confounder excluded)" = "#E69F00", "SS (confounder included)" = "#D55E00", "2D-SS" = "#0072B2")
grp_lty  <- c("SS (confounder excluded)" = "dashed", "SS (confounder included)" = "solid", "2D-SS" = "solid")
grp_shp  <- c("SS (confounder excluded)" = 17, "SS (confounder included)" = 16, "2D-SS" = 8)
d5 <- cond_summary; d5$group <- factor(d5$group, levels = grp_levels); d5$K_f <- lab_K(d5$K)
p5 <- ggplot(d5, aes(delta_test, mse_mean, color = group, shape = group, linetype = group, fill = group)) +
  geom_ribbon(aes(ymin = pmax(mse_mean - mse_se, 1e-6), ymax = mse_mean + mse_se), alpha = 0.12, color = NA) +
  geom_line() + geom_point(size = 2) +
  facet_wrap(~ K_f) + scale_y_log10() +
  scale_color_manual(values = grp_cols) + scale_fill_manual(values = grp_cols) +
  scale_shape_manual(values = grp_shp) + scale_linetype_manual(values = grp_lty) +
  labs(x = TeX("$\\delta_{test}$"), y = "Prediction MSE (log scale)",
       color = NULL, fill = NULL, shape = NULL, linetype = NULL) + th
ggsave("plot_study2_surprise_curve_conditioned.pdf", p5, width = 11, height = 5)

cat("Saved CSVs (selection raw/summary, cp, mse raw/summary/conditioned) and plot_study2_*.pdf\n")
