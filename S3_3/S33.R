# Study 3: does the PFER bound of Theorem 8 hold, and where does the
# confounder rate differ? Runs everything and draws the figures.
#
# Design: Example 2 (Sigma_SS = I_5, beta_S = (-4,-2,.5,.6,.7), confounder
# covariance a = (1/8,1/4,1/2,1/2,1/2)), K confounders, embedded in p = 1000.
# n in {100,200}, K in {1,5}. The bound is checked over q in {5,10,20,40} and
# tau in {0.6,0.75,0.9}, all from the same complementary-pair runs (N = 150 pairs).
# Features that never enter the Lasso path are ranked at random, and glmnet runs
# to dfmax so that q features enter the path.
#
#   Figure 1: empirical PFER vs delta_obs (q = 5, tau = 0.6) against the bound.
#   Figure 2: empirical PFER vs bound for EVERY (n,K,delta_obs,q,tau,method) cell;
#             all points on or below the diagonal means the bound holds.
#   Figure 3: confounder's complementary-pairs rate vs delta_obs (q = 5).
#   Figure 4: confounder's rate vs q at delta_obs = 0 (reduction persists in q).
#
# Output: section3_study3_summary.csv and plot_study3_*.pdf

suppressWarnings(options(warn = -1))

required_pkgs <- c("glmnet", "parallel", "ggplot2", "latex2exp")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs, dependencies = TRUE)
invisible(lapply(required_pkgs, library, character.only = TRUE))

RNGkind("L'Ecuyer-CMRG")
set.seed(2026)
n_cores <- max(1, min(7, parallel::detectCores() - 1))

# let the path run until dfmax instead of stopping early on deviance criteria
glmnet.control(fdev = 0, devmax = 1)

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
p <- 1000
s <- 5
beta_S  <- c(-4, -2, 0.5, 0.6, 0.7)
a_conf  <- c(1/8, 1/4, 1/2, 1/2, 1/2)
idio_sd <- sqrt(1 - sum(a_conf^2))
true_support <- 1:s

n_grid         <- c(100, 200)
K_grid         <- c(1, 5)
delta_obs_grid <- c(0, 0.5, 1, 1.5, 2)
delta_grid     <- seq(0, 1, length.out = 5)    # 2D-SS noise grid
q_grid         <- c(5, 10, 20, 40)
tau_grid       <- c(0.6, 0.75, 0.9)
N_pairs        <- 150

q_main   <- 5      # q used in Figures 1 and 3
tau_main <- 0.6

# Okabe-Ito
method_colors <- c("SS" = "#D55E00", "2D-SS" = "#0072B2")
method_shapes <- c("SS" = 18, "2D-SS" = 16)
method_lty    <- c("SS" = "solid", "2D-SS" = "solid")

# ---------------------------------------------------------------------------
# Data-generating process
# ---------------------------------------------------------------------------
generate_data <- function(n_obs, K, delta_obs) {
  p_plain <- p - s - K
  X_S <- matrix(rnorm(n_obs * s), n_obs, s)
  base <- as.numeric(X_S %*% a_conf)
  X_C <- sapply(seq_len(K), function(k) base + idio_sd * rnorm(n_obs))
  X_plain <- matrix(rnorm(n_obs * p_plain), n_obs, p_plain)
  X0 <- cbind(X_S, X_C, X_plain)
  y <- as.numeric(X_S %*% beta_S + rnorm(n_obs))
  X <- if (delta_obs > 0) X0 + matrix(rnorm(n_obs * p, sd = delta_obs), n_obs, p) else X0
  list(X = X, y = y)
}

# ---------------------------------------------------------------------------
# Size-matched selector. Returns the RANK (1 = first to enter the Lasso path) of
# every feature; "selected with size q" is rank <= q, so all q are obtained from
# one fit. Ties and never-entered features are ranked at random.
# ---------------------------------------------------------------------------
path_ranks <- function(X, y, q_max) {
  fit <- glmnet(X, y, alpha = 1, dfmax = q_max + 10)
  nz <- as.matrix(fit$beta) != 0
  entered <- rowSums(nz) > 0
  entry <- ifelse(entered, max.col(nz * 1, ties.method = "first"), Inf)
  ord <- order(entry, runif(length(entry)))
  r <- integer(length(entry)); r[ord] <- seq_along(ord)
  r
}

# One complementary pair at fixed (n, K, delta_obs): fresh dataset, two disjoint
# halves, delta-grid average, pair average. Returns a p x length(q_grid) matrix.
one_pair_fbar <- function(n, K, delta_obs, grid) {
  n1 <- floor(n / 2)
  dat <- generate_data(n, K, delta_obs)
  X <- dat$X; y <- dat$y
  perm <- sample(n)
  halves <- list(perm[seq_len(n1)], perm[(n1 + 1):(2 * n1)])
  q_max <- max(q_grid)
  acc <- matrix(0, p, length(q_grid))
  for (idx in halves) {
    Xb <- X[idx, , drop = FALSE]; yb <- y[idx]
    for (delta in grid) {
      Xj <- if (delta > 0) Xb + matrix(rnorm(n1 * p, sd = delta), n1, p) else Xb
      r <- path_ranks(Xj, yb, q_max)
      acc <- acc + sapply(q_grid, function(q) as.numeric(r <= q))
    }
  }
  acc / (length(halves) * length(grid))
}

# ---------------------------------------------------------------------------
# One cell (n, K, delta_obs, method): N_pairs pairs, then summaries for every
# (q, tau).
# ---------------------------------------------------------------------------
run_cell <- function(n, K, delta_obs, method) {
  grid <- if (method == "SS") 0 else delta_grid
  conf_idx <- (s + 1):(s + K)
  fb <- lapply(seq_len(N_pairs), function(b) one_pair_fbar(n, K, delta_obs, grid))
  out <- list()
  for (iq in seq_along(q_grid)) {
    q <- q_grid[iq]
    fmat <- t(sapply(fb, function(M) M[, iq]))      # N_pairs x p
    conf_pair <- rowMeans(fmat[, conf_idx, drop = FALSE])
    for (tau in tau_grid) {
      V <- rowSums(fmat[, -true_support, drop = FALSE] >= tau)
      out[[length(out) + 1]] <- data.frame(
        n = n, K = K, delta_obs = delta_obs, method = method, q = q, tau = tau,
        empirical_PFER = mean(V), PFER_se = sd(V) / sqrt(N_pairs),
        theoretical_bound = q^2 / ((p - s) * (2 * tau - 1)),
        confounder_rate = mean(conf_pair), confounder_se = sd(conf_pair) / sqrt(N_pairs))
    }
  }
  do.call(rbind, out)
}

tasks <- expand.grid(n = n_grid, K = K_grid, delta_obs = delta_obs_grid,
                     method = c("SS", "2D-SS"), stringsAsFactors = FALSE)
cat("Running", nrow(tasks), "cells x", N_pairs, "pairs on", n_cores, "cores...\n")
res <- mclapply(seq_len(nrow(tasks)), function(i)
  run_cell(tasks$n[i], tasks$K[i], tasks$delta_obs[i], tasks$method[i]),
  mc.cores = n_cores)
results_df <- do.call(rbind, res); rownames(results_df) <- NULL
write.csv(results_df, "section3_study3_summary.csv", row.names = FALSE)

# Does the bound hold? (point estimate and +2 SE)
results_df$holds     <- results_df$empirical_PFER <= results_df$theoretical_bound
results_df$holds_2se <- results_df$empirical_PFER - 2 * results_df$PFER_se <= results_df$theoretical_bound
cat(sprintf("Bound holds (point estimate) in %d / %d cells; max ratio empirical/bound = %.3f\n",
            sum(results_df$holds), nrow(results_df),
            max(results_df$empirical_PFER / results_df$theoretical_bound)))
print(aggregate(I(empirical_PFER / theoretical_bound) ~ q + tau, results_df, max))

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
results_df$method <- factor(results_df$method, levels = c("SS", "2D-SS"))
lab_n <- function(x) factor(paste0("n = ", x), levels = paste0("n = ", n_grid))
lab_K <- function(x) factor(ifelse(x == 1, "1 confounder", paste0(x, " confounders")),
                            levels = ifelse(K_grid == 1, "1 confounder", paste0(K_grid, " confounders")))
results_df$n_f <- lab_n(results_df$n); results_df$K_f <- lab_K(results_df$K)
th <- theme_bw(base_size = 13) + theme(legend.position = "bottom", strip.text = element_text(size = 10))

# Fig 1: PFER vs delta_obs, q = q_main, tau = tau_main
d1 <- subset(results_df, q == q_main & tau == tau_main)
bound1 <- unique(d1$theoretical_bound)
f1 <- ggplot(d1, aes(delta_obs, empirical_PFER, color = method, shape = method, fill = method)) +
  geom_ribbon(aes(ymin = pmax(empirical_PFER - PFER_se, 0), ymax = empirical_PFER + PFER_se),
              alpha = 0.15, color = NA) +
  geom_line() + geom_point(size = 2) +
  geom_hline(yintercept = bound1, linetype = "dashed") +
  facet_grid(K_f ~ n_f) +
  scale_color_manual(values = method_colors) + scale_fill_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  labs(x = TeX("$\\delta_{obs}$"), y = "Empirical PFER (dashed: Theorem 8 bound)",
       color = NULL, fill = NULL, shape = NULL) + th
ggsave("plot_study3_pfer_bound.pdf", f1, width = 9, height = 6)

# Fig 2: empirical vs bound, every cell
d2 <- results_df
lim <- max(d2$theoretical_bound)
f2 <- ggplot(d2, aes(theoretical_bound, empirical_PFER, color = method, shape = factor(tau))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = c("0.6" = 16, "0.75" = 17, "0.9" = 15)) +
  scale_x_sqrt(limits = c(0, lim)) + scale_y_sqrt(limits = c(0, lim)) +
  coord_fixed() +
  labs(x = "Theorem 8 bound", y = "Empirical PFER (both axes square-root scale)",
       color = NULL, shape = expression(tau)) + th
ggsave("plot_study3_bound_all_cells.pdf", f2, width = 6.5, height = 6.5)

# Fig 3: confounder rate vs delta_obs, q = q_main
d3 <- subset(results_df, q == q_main & tau == tau_main)
f3 <- ggplot(d3, aes(delta_obs, confounder_rate, color = method, shape = method, fill = method)) +
  geom_ribbon(aes(ymin = pmax(confounder_rate - confounder_se, 0),
                  ymax = confounder_rate + confounder_se), alpha = 0.15, color = NA) +
  geom_line() + geom_point(size = 2) +
  facet_grid(K_f ~ n_f) +
  scale_color_manual(values = method_colors) + scale_fill_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  labs(x = TeX("$\\delta_{obs}$"),
       y = TeX("Confounder's complementary-pairs rate, $\\bar{f}^{cp}_k$"),
       color = NULL, fill = NULL, shape = NULL) + th
ggsave("plot_study3_confounder_rate.pdf", f3, width = 9, height = 6)

# Fig 4: confounder rate vs q at delta_obs = 0
d4 <- subset(results_df, delta_obs == 0 & tau == tau_main)
f4 <- ggplot(d4, aes(q, confounder_rate, color = method, shape = method, fill = method)) +
  geom_ribbon(aes(ymin = pmax(confounder_rate - confounder_se, 0),
                  ymax = confounder_rate + confounder_se), alpha = 0.15, color = NA) +
  geom_line() + geom_point(size = 2) +
  facet_grid(K_f ~ n_f) + scale_x_log10(breaks = q_grid) +
  scale_color_manual(values = method_colors) + scale_fill_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  labs(x = "q (selector size, log scale)",
       y = TeX("Confounder's complementary-pairs rate, $\\bar{f}^{cp}_k$"),
       color = NULL, fill = NULL, shape = NULL) + th
ggsave("plot_study3_confounder_vs_q.pdf", f4, width = 9, height = 6)

cat("Saved: section3_study3_summary.csv and plot_study3_{pfer_bound,bound_all_cells,confounder_rate,confounder_vs_q}.pdf\n")
