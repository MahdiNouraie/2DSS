# Study 5: when does the largest-gap rule fail, and what should replace it?
#
# Theorem "gap" needs the boundary gap f_(s) - f_(s+1) to exceed every
# within-group gap. With strong AND weak signals the largest gap can fall
# inside S, and the rule underestimates s. Here we
#   (a) construct the failure (weak strength w = 0.5 is the original design),
#   (b) sweep the weak-signal strength to show where it appears and disappears,
#   (c) compare the largest-gap rule with a known-s rule (top-s) and with
#       fixed thresholds tau in {0.15, 0.3, 0.6},
#   (d) do this for both the plain aggregation (SS, delta = 0 only) and the
#       2D-SS aggregation (grid 0, .25, .5, .75, 1). No confounder, no extra
#       measurement error: this isolates the thresholding step.
#
# Design: n = 100, p = 500, s = 5, Sigma = I,
#         beta = (6, 5, w*1, w*0.9, w*0.8), w in {0.3, 0.5, 0.8, 1.2, 2}.
#         (w = 0.5 gives (0.5, 0.45, 0.4).)
#
# Colours: Okabe-Ito (colour-blind safe). SS #D55E00, 2D-SS #0072B2 as in the
# other studies; rules use the remaining Okabe-Ito colours plus linetype/shape.

suppressWarnings(options(warn = -1))

required_pkgs <- c("glmnet", "parallel", "ggplot2", "latex2exp")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs, dependencies = TRUE)
invisible(lapply(required_pkgs, library, character.only = TRUE))

RNGkind("L'Ecuyer-CMRG")
set.seed(2026)

n_cores <- 7

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
n <- 100
p <- 500
s <- 5
true_support <- 1:s
strong_idx <- 1:2
weak_idx <- 3:5

weak_grid <- c(0.3, 0.5, 0.8, 1.2, 2)   # w: weak-signal strength
weak_shape <- c(1, 0.9, 0.8)            # beta_weak = w * weak_shape
beta_strong <- c(6, 5)

n1 <- floor(n / 2)
lambda_half <- sqrt(log(p) / n1)

B <- 100
delta_grid <- seq(0, 1, length.out = 5)  # 2D-SS grid (SS uses delta = 0 only)
n_rep <- 100

tau_set <- c(0.15, 0.3, 0.6)             # pre-specified fixed thresholds
profile_len <- 12                        # ranks kept for the profile plot

# ---------------------------------------------------------------------------
# Data and aggregation
# ---------------------------------------------------------------------------
generate_data <- function(w) {
  beta <- c(beta_strong, w * weak_shape)
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X[, true_support] %*% beta + rnorm(n))
  list(X = X, y = y)
}

fit_coef <- function(X, y, lam) {
  fit <- glmnet(X, y, alpha = 1, lambda = lam, standardize = TRUE)
  as.numeric(coef(fit))[-1]
}

# Same jitter convention as Study 4: N(0, delta^2) noise added to the half-sample
# design, then standardise and fit at the fixed lambda.
aggregate_frequencies <- function(X, y, grid, n1_use, B_use, lam) {
  n_obs <- nrow(X); p_full <- ncol(X)
  freq_mat <- matrix(0, length(grid), p_full)
  for (i in seq_along(grid)) {
    delta <- grid[i]
    cnt <- numeric(p_full)
    for (b in seq_len(B_use)) {
      idx <- sample(n_obs, n1_use)
      Xj <- X[idx, , drop = FALSE]
      if (delta > 0) Xj <- Xj + matrix(rnorm(n1_use * p_full, sd = delta), n1_use, p_full)
      cnt <- cnt + (fit_coef(scale(Xj), y[idx], lam) != 0)
    }
    freq_mat[i, ] <- cnt / B_use
  }
  colMeans(freq_mat)
}

largest_gap_select <- function(fbar) {
  pp <- length(fbar)
  ord <- order(fbar, decreasing = TRUE)
  sorted_f <- fbar[ord]
  gaps <- sorted_f[-pp] - sorted_f[-1]
  khat <- which.max(gaps)
  tau_hat <- 0.5 * (sorted_f[khat] + sorted_f[khat + 1])
  list(Shat = which(fbar >= tau_hat), shat = khat)
}

top_s_select <- function(fbar, s_known) {
  # known model size: take the s_known most frequently selected features
  # (random tie-break so ties cannot favour the true indices)
  ord <- order(fbar, runif(length(fbar)), decreasing = TRUE)
  ord[seq_len(s_known)]
}

metrics_of <- function(Shat, S) {
  tp <- length(intersect(Shat, S))
  fp <- length(setdiff(Shat, S))
  fn <- length(setdiff(S, Shat))
  c(TPR = tp / length(S),
    FDR = if (length(Shat) > 0) fp / length(Shat) else 0,
    F1 = if (tp > 0) 2 * tp / (2 * tp + fp + fn) else 0,
    size = length(Shat),
    exact = as.numeric(setequal(Shat, S)),
    strong = sum(strong_idx %in% Shat),
    weak = sum(weak_idx %in% Shat))
}

rule_list <- function(fbar) {
  gp <- largest_gap_select(fbar)
  out <- list("Largest-gap" = gp$Shat, "Known s (top-s)" = top_s_select(fbar, s))
  for (tt in tau_set) out[[paste0("Fixed tau=", tt)]] <- which(fbar >= tt)
  out
}

# ---------------------------------------------------------------------------
# One replicate (one weak level, one rep): both aggregations, all rules
# ---------------------------------------------------------------------------
run_one <- function(job) {
  w <- job$w; rep_id <- job$rep
  dat <- generate_data(w)
  aggs <- list(SS = aggregate_frequencies(dat$X, dat$y, 0, n1, B, lambda_half),
               "2D-SS" = aggregate_frequencies(dat$X, dat$y, delta_grid, n1, B, lambda_half))
  rows <- list(); prof <- list()
  for (an in names(aggs)) {
    fbar <- aggs[[an]]
    rl <- rule_list(fbar)
    for (rn in names(rl)) {
      m <- metrics_of(rl[[rn]], true_support)
      rows[[length(rows) + 1]] <- data.frame(
        weak = w, rep = rep_id, aggregation = an, rule = rn,
        TPR = m["TPR"], FDR = m["FDR"], F1 = m["F1"], size = m["size"],
        exact = m["exact"], strong_found = m["strong"], weak_found = m["weak"],
        row.names = NULL)
    }
    sf <- sort(fbar, decreasing = TRUE)
    prof[[length(prof) + 1]] <- data.frame(
      weak = w, rep = rep_id, aggregation = an,
      rank = seq_len(profile_len), freq = sf[seq_len(profile_len)])
  }
  list(res = do.call(rbind, rows), prof = do.call(rbind, prof))
}

jobs <- expand.grid(w = weak_grid, rep = seq_len(n_rep))
jobs <- split(jobs, seq_len(nrow(jobs)))
jobs <- lapply(jobs, as.list)

cat("Running", length(jobs), "jobs (", length(weak_grid), "weak levels x", n_rep, "reps ) on",
    n_cores, "cores\n")
out <- mclapply(jobs, run_one, mc.cores = n_cores, mc.set.seed = TRUE)
ok <- !vapply(out, function(z) inherits(z, "try-error") || is.null(z$res), logical(1))
if (any(!ok)) cat("WARNING:", sum(!ok), "jobs failed and are dropped\n")
out <- out[ok]
res <- do.call(rbind, lapply(out, `[[`, "res"))
prof <- do.call(rbind, lapply(out, `[[`, "prof"))
rownames(res) <- NULL

# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------
se <- function(x) sd(x) / sqrt(length(x))
key <- interaction(res$weak, res$aggregation, res$rule, drop = TRUE)
summary_df <- do.call(rbind, lapply(split(res, key), function(d) {
  data.frame(weak = d$weak[1], aggregation = d$aggregation[1], rule = d$rule[1], n = nrow(d),
             TPR_mean = mean(d$TPR), TPR_se = se(d$TPR),
             FDR_mean = mean(d$FDR), FDR_se = se(d$FDR),
             F1_mean = mean(d$F1), F1_se = se(d$F1),
             size_mean = mean(d$size), size_se = se(d$size),
             size_correct = mean(d$size == s),
             exact_mean = mean(d$exact),
             strong_mean = mean(d$strong_found), weak_mean = mean(d$weak_found))
}))
rownames(summary_df) <- NULL
summary_df <- summary_df[order(summary_df$aggregation, summary_df$weak, summary_df$rule), ]

# Paired differences in F1: largest-gap minus each other rule, by replicate
wide <- reshape(res[, c("weak", "rep", "aggregation", "rule", "F1")],
                idvar = c("weak", "rep", "aggregation"), timevar = "rule", direction = "wide")
paired <- do.call(rbind, lapply(setdiff(unique(res$rule), "Largest-gap"), function(rn) {
  d <- wide[["F1.Largest-gap"]] - wide[[paste0("F1.", rn)]]
  agg <- aggregate(d, by = list(weak = wide$weak, aggregation = wide$aggregation),
                   FUN = function(x) c(mean = mean(x), se = se(x)))
  data.frame(weak = agg$weak, aggregation = agg$aggregation,
             versus = rn, diff_mean = agg$x[, "mean"], diff_se = agg$x[, "se"])
}))

print(summary_df[summary_df$aggregation == "SS", c("weak", "rule", "TPR_mean", "FDR_mean",
                                                   "F1_mean", "size_mean", "size_correct")],
      digits = 3, row.names = FALSE)

write.csv(res, "section3_study5_raw_results.csv", row.names = FALSE)
write.csv(summary_df, "section3_study5_summary.csv", row.names = FALSE)
write.csv(paired, "section3_study5_paired.csv", row.names = FALSE)

prof_sum <- do.call(rbind, lapply(split(prof, interaction(prof$weak, prof$aggregation, prof$rank, drop = TRUE)),
                                  function(d) data.frame(weak = d$weak[1], aggregation = d$aggregation[1], rank = d$rank[1],
                                                         mean = mean(d$freq), lo = quantile(d$freq, 0.1), hi = quantile(d$freq, 0.9))))
rownames(prof_sum) <- NULL
write.csv(prof_sum, "section3_study5_profile.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# Figures (Okabe-Ito)
# ---------------------------------------------------------------------------
agg_cols <- c("SS" = "#D55E00", "2D-SS" = "#0072B2")
rule_levels <- c("Largest-gap", "Known s (top-s)", paste0("Fixed tau=", tau_set))
rule_cols <- c("Largest-gap" = "#000000", "Known s (top-s)" = "#009E73",
               "Fixed tau=0.15" = "#999999", "Fixed tau=0.3" = "#E69F00",
               "Fixed tau=0.6" = "#CC79A7")
rule_lty <- c("Largest-gap" = "solid", "Known s (top-s)" = "solid",
              "Fixed tau=0.15" = "dotted", "Fixed tau=0.3" = "dashed", "Fixed tau=0.6" = "dotdash")
rule_shp <- c("Largest-gap" = 16, "Known s (top-s)" = 15,
              "Fixed tau=0.15" = 4, "Fixed tau=0.3" = 17, "Fixed tau=0.6" = 18)
rule_lab <- c("Largest-gap" = "Largest-gap", "Known s (top-s)" = "Known s (top-s)",
              "Fixed tau=0.15" = "Fixed 0.15", "Fixed tau=0.3" = "Fixed 0.3", "Fixed tau=0.6" = "Fixed 0.6")

theme_s5 <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"))

# 1. Sorted-frequency profile (mean and 10-90% band over replicates)
prof_sum$wlab <- factor(paste0("w = ", prof_sum$weak), levels = paste0("w = ", weak_grid))
prof_sum$aggregation <- factor(prof_sum$aggregation, levels = c("SS", "2D-SS"))
p1 <- ggplot(prof_sum, aes(x = rank, y = mean, colour = aggregation, fill = aggregation)) +
  geom_vline(xintercept = s + 0.5, linetype = "dotted", colour = "grey30") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
  facet_grid(aggregation ~ wlab) +
  scale_colour_manual(values = agg_cols) + scale_fill_manual(values = agg_cols) +
  scale_x_continuous(breaks = c(1, 5, 10)) +
  labs(x = "Rank (sorted by selection frequency)",
       y = TeX("Mean sorted frequency $\\bar{f}_{(k)}$ (band: 10--90%)")) +
  theme_s5 + theme(legend.position = "none")
ggsave("plot_study5_profile.pdf", p1, width = 11, height = 5)

# 2. TPR / FDR / F1 against weak strength, by rule (rows: metric, cols: aggregation)
long <- do.call(rbind, lapply(c("TPR", "FDR", "F1"), function(m)
  data.frame(summary_df[, c("weak", "aggregation", "rule")], metric = m,
             mean = summary_df[[paste0(m, "_mean")]], se = summary_df[[paste0(m, "_se")]])))
long$metric <- factor(long$metric, levels = c("TPR", "FDR", "F1"))
long$rule <- factor(long$rule, levels = rule_levels)
long$aggregation <- factor(long$aggregation, levels = c("SS", "2D-SS"))
p2 <- ggplot(long, aes(x = weak, y = mean, colour = rule, linetype = rule, shape = rule)) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(mean - se, 0), ymax = pmin(mean + se, 1)), width = 0.03, linewidth = 0.3) +
  facet_grid(metric ~ aggregation) +
  scale_colour_manual(values = rule_cols, labels = rule_lab) +
  scale_linetype_manual(values = rule_lty, labels = rule_lab) +
  scale_shape_manual(values = rule_shp, labels = rule_lab) +
  scale_x_continuous(breaks = weak_grid) +
  labs(x = "Weak-signal strength w", y = NULL, colour = NULL, linetype = NULL, shape = NULL) +
  theme_s5 + guides(colour = guide_legend(nrow = 1))
ggsave("plot_study5_rules.pdf", p2, width = 9, height = 8)

# 3. Selected size against weak strength (true s = 5)
sz <- summary_df
sz$rule <- factor(sz$rule, levels = rule_levels)
sz$aggregation <- factor(sz$aggregation, levels = c("SS", "2D-SS"))
p3 <- ggplot(sz, aes(x = weak, y = size_mean, colour = rule, linetype = rule, shape = rule)) +
  geom_hline(yintercept = s, colour = "grey30", linetype = "dotted") +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = size_mean - size_se, ymax = size_mean + size_se), width = 0.03, linewidth = 0.3) +
  facet_wrap(~ aggregation) +
  scale_colour_manual(values = rule_cols, labels = rule_lab) +
  scale_linetype_manual(values = rule_lty, labels = rule_lab) +
  scale_shape_manual(values = rule_shp, labels = rule_lab) +
  scale_x_continuous(breaks = weak_grid) + scale_y_log10() +
  labs(x = "Weak-signal strength w", y = "Mean number selected (log scale; dotted: s = 5)",
       colour = NULL, linetype = NULL, shape = NULL) +
  theme_s5 + guides(colour = guide_legend(nrow = 1))
ggsave("plot_study5_size.pdf", p3, width = 9, height = 4.5)

cat("Saved: section3_study5_{raw_results,summary,paired,profile}.csv,",
    "plot_study5_{profile,rules,size}.pdf\n")