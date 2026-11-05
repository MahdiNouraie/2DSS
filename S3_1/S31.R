# Study 1: single plain fit vs single jittered fit (Corollary "benefit").
# One script: runs the simulation, writes the CSVs and draws the paper figures.
#
# Designs (each with one confounder)
#   A  Example 2: Sigma_SS = I_5, a = (1/8,1/4,1/2,1/2,1/2); Lasso IC violated
#      and repaired by isotropic jitter.
#   B  Example 1: correlated Sigma_SS (rho = 0.5), a = (0.8,1.2); jitter hurts.
#   C  Control: Sigma_SS = I_5, a = (.05,.1,.1,.1,.1); Lasso IC holds.
# Outcome: P(exact recovery of S with correct signs, confounder excluded) of one
# Lasso fit at lambda = c * sqrt(log(p)/n1).
# Part 1: grid over (design, n1, nu), c = 1.   Part 2: tuning constant c at n1 = 500.
# Colours: Okabe-Ito / viridis (colour-blind safe), with shapes and line types.
#
# Output: section3_study0_{results,summary,cnu,lambda_sensitivity}.csv and
#         plot_corollary1_{curves,vs_n1,lambda,cnu}.pdf

suppressWarnings(options(warn = -1))

required_pkgs <- c("glmnet", "parallel", "ggplot2", "latex2exp", "viridisLite")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs, dependencies = TRUE)
invisible(lapply(required_pkgs, library, character.only = TRUE))

RNGkind("L'Ecuyer-CMRG")
set.seed(2026)
n_cores <- max(1, min(7, parallel::detectCores() - 1))

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
n1_grid      <- c(25, 50, 100, 200, 500, 1000, 2000, 5000)
nu_grid      <- seq(0, 1.5, by = 0.1)
n_trials     <- 1000
n1_lam       <- 500
c_lam_grid   <- c(0.5, 0.75, 1, 1.5, 2)
n_trials_lam <- 500

# ---------------------------------------------------------------------------
# Designs
# ---------------------------------------------------------------------------
make_design <- function(name, Sigma_SS, beta_S, a_conf) {
  s <- length(beta_S)
  w <- as.numeric(solve(Sigma_SS, a_conf))
  idio_sd <- sqrt(max(1 - as.numeric(a_conf %*% w), 0.1))
  list(name = name, s = s, p = s + 1, Sigma_SS = Sigma_SS, L = chol(Sigma_SS),
       beta_S = beta_S, z = sign(beta_S), a_conf = a_conf, w = w, idio_sd = idio_sd)
}
b_ex2 <- c(-4, -2, 0.5, 0.6, 0.7)
designs <- list(
  make_design("Example 2: IC violated, repairable", diag(5), b_ex2, c(1/8, 1/4, 1/2, 1/2, 1/2)),
  make_design("Example 1: IC violated, jitter hurts",
              matrix(c(1, 0.5, 0.5, 1), 2, 2), c(2.5, 1.5), c(0.8, 1.2)),
  make_design("Control: IC holds", diag(5), b_ex2, c(0.05, 0.1, 0.1, 0.1, 0.1))
)
design_names <- vapply(designs, function(d) d$name, character(1))

# Population IC statistic c_nu = |a' T^{-1} (nu^2 beta + lambda z)| / lambda,
# T = Sigma_SS + nu^2 I  (= |a'(z + nu^2 b)|/(1+nu^2) when Sigma_SS = I)
c_nu_fun <- function(d, nu, lam) {
  T <- d$Sigma_SS + nu^2 * diag(d$s)
  abs(sum(d$a_conf * solve(T, nu^2 * d$beta_S + lam * d$z))) / lam
}
cnu_df <- do.call(rbind, lapply(designs, function(d)
  do.call(rbind, lapply(n1_grid, function(n1) {
    lam <- sqrt(log(d$p) / n1)
    data.frame(design = d$name, n1 = n1, nu = nu_grid, lambda = lam,
               c_nu = sapply(nu_grid, function(nu) c_nu_fun(d, nu, lam)))
  }))))
write.csv(cnu_df, "section3_study0_cnu.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------
one_trial <- function(d, nu, n1, lam) {
  s <- d$s
  X_S <- matrix(rnorm(n1 * s), n1, s) %*% d$L
  X_C <- as.numeric(X_S %*% d$w) + d$idio_sd * rnorm(n1)
  X <- cbind(X_S, X_C)
  y <- as.numeric(X_S %*% d$beta_S + rnorm(n1))
  if (nu > 0) X <- X + matrix(rnorm(n1 * d$p, sd = nu), n1, d$p)
  bh <- as.numeric(coef(glmnet(X, y, alpha = 1, lambda = lam, standardize = TRUE)))[-1]
  sign_ok <- all(sign(bh[1:s]) == d$z)
  no_fp   <- bh[s + 1] == 0
  c(as.numeric(sign_ok && no_fp), as.numeric(sign_ok), as.numeric(no_fp))
}
binom_se <- function(p_hat, n) sqrt(pmax(p_hat * (1 - p_hat), 0) / n)

run_cell <- function(d, n1, c_lam, trials) {
  lam <- c_lam * sqrt(log(d$p) / n1)
  do.call(rbind, lapply(nu_grid, function(nu) {
    m <- rowMeans(vapply(seq_len(trials), function(i) one_trial(d, nu, n1, lam), numeric(3)))
    data.frame(design = d$name, n1 = n1, c_lam = c_lam, nu = nu, lambda = lam,
               exact_recovery = m[1], signs_ok = m[2], no_false_positive = m[3],
               se = binom_se(m[1], trials), trials = trials)
  }))
}

# Part 1
tasks <- expand.grid(di = seq_along(designs), n1 = n1_grid)
results_df <- do.call(rbind, mclapply(seq_len(nrow(tasks)), function(i)
  run_cell(designs[[tasks$di[i]]], tasks$n1[i], 1, n_trials), mc.cores = n_cores))
rownames(results_df) <- NULL
write.csv(results_df, "section3_study0_results.csv", row.names = FALSE)

# Part 2
tasks2 <- expand.grid(di = seq_along(designs), c_lam = c_lam_grid)
lam_df <- do.call(rbind, mclapply(seq_len(nrow(tasks2)), function(i)
  run_cell(designs[[tasks2$di[i]]], n1_lam, tasks2$c_lam[i], n_trials_lam), mc.cores = n_cores))
rownames(lam_df) <- NULL
write.csv(lam_df, "section3_study0_lambda_sensitivity.csv", row.names = FALSE)

# Summary: plain fit vs best jittered fit (best nu chosen post hoc)
summary_df <- do.call(rbind, lapply(design_names, function(dn) do.call(rbind, lapply(n1_grid, function(n1) {
  sub <- results_df[results_df$design == dn & results_df$n1 == n1, ]
  r0 <- sub[sub$nu == 0, ]; pos <- sub[sub$nu > 0, ]; bi <- which.max(pos$exact_recovery)
  data.frame(design = dn, n1 = n1,
             recovery_at_nu0 = r0$exact_recovery, se_at_nu0 = r0$se,
             best_nu = pos$nu[bi], recovery_at_best_nu = pos$exact_recovery[bi],
             se_at_best_nu = pos$se[bi],
             gap = pos$exact_recovery[bi] - r0$exact_recovery,
             gap_se = sqrt(pos$se[bi]^2 + r0$se^2))
}))))
write.csv(summary_df, "section3_study0_summary.csv", row.names = FALSE)
print(summary_df, digits = 3, row.names = FALSE)

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
ok <- c(blue = "#0072B2", vermillion = "#D55E00")
pal_fit <- c("Plain fit (nu = 0)" = unname(ok["vermillion"]),
             "Jittered fit (best nu)" = unname(ok["blue"]))
n1_cols  <- setNames(viridisLite::viridis(length(n1_grid), option = "D", begin = 0, end = 0.9), n1_grid)
lam_cols <- setNames(viridisLite::viridis(length(c_lam_grid), option = "C", begin = 0.05, end = 0.85), c_lam_grid)
th <- theme_bw(base_size = 13) +
  theme(legend.position = "bottom", strip.text = element_text(size = 10))
dfac <- function(x) factor(x, levels = design_names)

# Figure A: recovery vs nu, by n1
results_df$n1_f <- factor(results_df$n1, levels = n1_grid)
results_df$design_f <- dfac(results_df$design)
p1 <- ggplot(results_df, aes(nu, exact_recovery, color = n1_f, fill = n1_f, shape = n1_f)) +
  geom_ribbon(aes(ymin = pmax(exact_recovery - se, 0), ymax = pmin(exact_recovery + se, 1)),
              alpha = 0.12, color = NA) +
  geom_line() + geom_point(size = 1.6) +
  facet_wrap(~ design_f, ncol = 3) +
  scale_color_manual(values = n1_cols) + scale_fill_manual(values = n1_cols) +
  scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 4, 18)) + ylim(0, 1) +
  labs(x = TeX("$\\nu$"), y = "P(exact sign-consistent recovery)",
       color = TeX("$n_1$"), fill = TeX("$n_1$"), shape = TeX("$n_1$")) + th
ggsave("plot_corollary1_curves.pdf", p1, width = 13, height = 4.8)

# Figure B: plain vs best jittered fit, by n1
ls <- data.frame(design_f = dfac(rep(summary_df$design, 2)), n1 = rep(summary_df$n1, 2),
                 fit = factor(rep(names(pal_fit), each = nrow(summary_df)), levels = names(pal_fit)),
                 recovery = c(summary_df$recovery_at_nu0, summary_df$recovery_at_best_nu),
                 se = c(summary_df$se_at_nu0, summary_df$se_at_best_nu))
fit_tex <- TeX(c("Plain fit ($\\nu = 0$)", "Jittered fit (best $\\nu$)"))
p2 <- ggplot(ls, aes(n1, recovery, color = fit, fill = fit, shape = fit, linetype = fit)) +
  geom_ribbon(aes(ymin = pmax(recovery - se, 0), ymax = pmin(recovery + se, 1)),
              alpha = 0.15, color = NA) +
  geom_line() + geom_point(size = 2) +
  facet_wrap(~ design_f, ncol = 3) + scale_x_log10() +
  scale_color_manual(values = pal_fit, labels = fit_tex) +
  scale_fill_manual(values = pal_fit, labels = fit_tex) +
  scale_shape_manual(values = c(16, 17), labels = fit_tex) +
  scale_linetype_manual(values = c("solid", "longdash"), labels = fit_tex) + ylim(0, 1) +
  labs(x = TeX("$n_1$ (log scale)"), y = "P(exact sign-consistent recovery)",
       color = NULL, fill = NULL, shape = NULL, linetype = NULL) + th
ggsave("plot_corollary1_vs_n1.pdf", p2, width = 13, height = 4.8)

# Figure C: sensitivity to the tuning constant (n1 = n1_lam)
lam_df$c_lam_f <- factor(lam_df$c_lam, levels = c_lam_grid)
lam_df$design_f <- dfac(lam_df$design)
p3 <- ggplot(lam_df, aes(nu, exact_recovery, color = c_lam_f, shape = c_lam_f)) +
  geom_line() + geom_point(size = 1.6) +
  facet_wrap(~ design_f, ncol = 3) + scale_color_manual(values = lam_cols) + ylim(0, 1) +
  labs(x = TeX("$\\nu$"), y = "P(exact sign-consistent recovery)",
       color = TeX("$c = \\lambda / \\sqrt{\\log(p)/n_1}$"),
       shape = TeX("$c = \\lambda / \\sqrt{\\log(p)/n_1}$")) + th
ggsave("plot_corollary1_lambda.pdf", p3, width = 13, height = 4.8)

# Figure D (optional): population c_nu, capped at 3 (Example 1 grows past 100)
cnu_df$n1_f <- factor(cnu_df$n1, levels = n1_grid)
cnu_df$design_f <- dfac(cnu_df$design)
p4 <- ggplot(cnu_df, aes(nu, c_nu, color = n1_f)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35") + geom_line() +
  facet_wrap(~ design_f, ncol = 3) + scale_color_manual(values = n1_cols) +
  coord_cartesian(ylim = c(0, 3)) +
  labs(x = TeX("$\\nu$"), y = TeX("$c_\\nu$"), color = TeX("$n_1$")) + th
ggsave("plot_corollary1_cnu.pdf", p4, width = 13, height = 4.8)

cat("Saved 4 CSVs and plot_corollary1_{curves,vs_n1,lambda,cnu}.pdf\n")