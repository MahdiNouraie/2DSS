# Study 4 add-on (one script: runs the experiment, then writes all tables and figures).
# Is 2D-SS's lower F1 at low noise caused by the largest-gap rule rather than by the aggregation?
#
# Same design, settings and generator as S34.R (n = 100, p = 200, K in {1, 5, 0},
# delta_obs in {0, .5, 1, 1.5, 2}, B = 50, 2D-SS grid seq(0, 1, 5), 40 replicates).
# Only SS and 2D-SS are run (no baselines), each with four selection rules:
#   Largest-gap      : the rule used in Study 4 (reproduces the existing rows)
#   Known s (top-s)  : the s most frequent features (uses the true s; compare with the
#                      baselines' CEILINGS, not their package-tuned versions)
#   Fixed tau = 0.3  : pre-specified, no knowledge of s
#   Fixed tau = 0.6  : pre-specified, the value compatible with the error bound (tau > 1/2)
# Colours: Okabe-Ito. If section3_study4_summary.csv (from S34.R) is in the working
# directory it is used for GMUS (elbow) and the best baseline ceiling; otherwise the same
# values, stored below from the Study 4 run, are used. Either way the figures include them.
#
# Outputs: section3_study4b_raw_results.csv, section3_study4b_summary.csv,
#   section3_study4b_combined.csv, section3_study4b_paired.csv,
#   plot_study4b_f1.pdf, plot_study4b_confounder.pdf

suppressWarnings(options(warn = -1))
required_pkgs <- c("glmnet", "parallel", "ggplot2")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs, dependencies = TRUE)
invisible(lapply(required_pkgs, library, character.only = TRUE))

RNGkind("L'Ecuyer-CMRG")
set.seed(2026)
n_cores <- max(1, min(7, parallel::detectCores() - 1))

n <- 100; p <- 200; s <- 5
beta_S  <- c(-4, -2, 0.5, 0.6, 0.7)
a_conf  <- c(1/8, 1/4, 1/2, 1/2, 1/2)
idio_sd <- sqrt(1 - sum(a_conf^2))
true_support <- 1:s

K_grid         <- c(1, 5, 0)
delta_obs_grid <- c(0, 0.5, 1, 1.5, 2)
delta_grid     <- seq(0, 1, length.out = 5)
m <- length(delta_grid)
B <- 50
B_SS <- m * B
n1 <- floor(n / 2)
lambda_half <- sqrt(log(p) / n1)
n_rep <- 40
tau_set <- c(0.3, 0.6)

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

fit_coef <- function(X, y, lam) {
  fit <- glmnet(X, y, alpha = 1, lambda = lam, standardize = TRUE)
  as.numeric(coef(fit))[-1]
}

aggregate_frequencies <- function(X, y, grid, n1_use, B_use, lam) {
  n_obs <- nrow(X); p_full <- ncol(X)
  freq_mat <- matrix(0, length(grid), p_full)
  for (i in seq_along(grid)) {
    delta <- grid[i]
    cnt <- numeric(p_full)
    for (b in seq_len(B_use)) {
      idx <- if (n1_use < n_obs) sample(n_obs, n1_use) else seq_len(n_obs)
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
  sf <- fbar[ord]
  khat <- which.max(sf[-pp] - sf[-1])
  which(fbar >= 0.5 * (sf[khat] + sf[khat + 1]))
}
top_s_select <- function(fbar, s_known) {
  ord <- order(fbar, runif(length(fbar)), decreasing = TRUE)   # random tie-break
  ord[seq_len(s_known)]
}

selection_metrics <- function(Shat, K) {
  tp <- length(intersect(Shat, true_support))
  fp <- length(setdiff(Shat, true_support))
  fn <- length(setdiff(true_support, Shat))
  conf_idx <- if (K > 0) (s + 1):(s + K) else integer(0)
  c(TPR = tp / s,
    FDR = if (length(Shat) > 0) fp / length(Shat) else 0,
    F1 = if (tp > 0) 2 * tp / (2 * tp + fp + fn) else 0,
    exact = as.numeric(tp == s && fp == 0),
    size = length(Shat),
    conf_any = if (K > 0) as.numeric(any(conf_idx %in% Shat)) else NA_real_)
}

run_one <- function(K, delta_obs, rep_id) {
  dat <- generate_data(K, delta_obs)
  X <- dat$X; y <- dat$y
  aggs <- list(SS = aggregate_frequencies(X, y, 0, n1, B_SS, lambda_half),
               "2D-SS" = aggregate_frequencies(X, y, delta_grid, n1, B, lambda_half))
  rows <- list()
  for (an in names(aggs)) {
    f <- aggs[[an]]
    sel <- list("Largest-gap" = largest_gap_select(f),
                "Known s (top-s)" = top_s_select(f, s))
    for (tt in tau_set) sel[[paste0("Fixed tau=", tt)]] <- which(f >= tt)
    for (rn in names(sel))
      rows[[length(rows) + 1]] <- data.frame(t(selection_metrics(sel[[rn]], K)), K = K,
                                             delta_obs = delta_obs, rep = rep_id, aggregation = an, rule = rn, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

tasks <- expand.grid(K = K_grid, delta_obs = delta_obs_grid, rep_id = seq_len(n_rep))
cat("Running", nrow(tasks), "replicates on", n_cores, "cores...\n")
out <- mclapply(seq_len(nrow(tasks)), function(i)
  run_one(tasks$K[i], tasks$delta_obs[i], tasks$rep_id[i]), mc.cores = n_cores, mc.set.seed = TRUE)
ok <- !sapply(out, function(r) inherits(r, "try-error") || is.null(r))
if (any(!ok)) cat("Warning:", sum(!ok), "replicates failed and were dropped.\n")
raw <- do.call(rbind, out[ok]); rownames(raw) <- NULL
write.csv(raw, "section3_study4b_raw_results.csv", row.names = FALSE)

se <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA else sd(x) / sqrt(length(x)) }
mn <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA else mean(x) }
mets <- c("TPR", "FDR", "F1", "exact", "size", "conf_any")
grp <- split(raw, list(raw$K, raw$delta_obs, raw$aggregation, raw$rule), drop = TRUE)
summ <- do.call(rbind, lapply(grp, function(d) {
  o <- data.frame(K = d$K[1], delta_obs = d$delta_obs[1], aggregation = d$aggregation[1],
                  rule = d$rule[1], n = nrow(d))
  for (mm in mets) { o[[paste0(mm, "_mean")]] <- mn(d[[mm]]); o[[paste0(mm, "_se")]] <- se(d[[mm]]) }
  o
}))
rownames(summ) <- NULL
write.csv(summ, "section3_study4b_summary.csv", row.names = FALSE)
print(summ[summ$aggregation == "2D-SS" & summ$K == 1,
           c("delta_obs", "rule", "TPR_mean", "FDR_mean", "F1_mean", "size_mean")], digits = 3, row.names = FALSE)

# ---------------------------------------------------------------------------
# Series: add-on (SS / 2D-SS under two rules), GMUS (elbow), best baseline ceiling
# ---------------------------------------------------------------------------
sm <- summ
keep <- sm[sm$rule %in% c("Largest-gap", "Known s (top-s)"), ]
keep$series <- paste0(keep$aggregation, ": ", ifelse(keep$rule == "Largest-gap", "largest-gap", "known s"))
cols_add <- c("K", "delta_obs", "series", "TPR_mean", "FDR_mean", "F1_mean", "F1_se",
              "size_mean", "conf_any_mean", "conf_any_se")
add <- keep[, cols_add]

have_s4 <- file.exists("section3_study4_summary.csv")
if (have_s4) {
  s4 <- read.csv("section3_study4_summary.csv")
  pick <- function(lbl, nm) {
    g <- s4[s4$label == lbl, ]
    data.frame(K = g$K, delta_obs = g$delta_obs, series = nm, TPR_mean = g$TPR_mean, FDR_mean = g$FDR_mean,
               F1_mean = g$F1_mean, F1_se = g$F1_se, size_mean = NA_real_,
               conf_any_mean = g$conf_any_mean, conf_any_se = g$conf_any_se)
  }
  gm <- pick("GMUS (elbow)", "GMUS (elbow)")
  ce <- s4[grepl("\\(ceiling\\)", s4$label), ]
  best <- do.call(rbind, lapply(split(ce, list(ce$K, ce$delta_obs), drop = TRUE), function(d) {
    d <- d[which.max(d$F1_mean), ]
    data.frame(K = d$K, delta_obs = d$delta_obs, series = "Best baseline ceiling", TPR_mean = d$TPR_mean,
               FDR_mean = d$FDR_mean, F1_mean = d$F1_mean, F1_se = d$F1_se, size_mean = NA_real_,
               conf_any_mean = d$conf_any_mean, conf_any_se = d$conf_any_se)
  }))
} else {
  cat("section3_study4_summary.csv not found: using the GMUS / ceiling values stored from the Study 4 run (40 replicates).\n")
  ref <- data.frame(
    K = c(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 5.0, 5.0, 5.0, 5.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 5.0, 5.0, 5.0, 5.0, 5.0),
    delta_obs = c(0.0, 0.5, 1.0, 1.5, 2.0, 0.0, 0.5, 1.0, 1.5, 2.0, 0.0, 0.5, 1.0, 1.5, 2.0, 0.0, 0.5, 1.0, 1.5, 2.0, 0.0, 0.5, 1.0, 1.5, 2.0, 0.0, 0.5, 1.0, 1.5, 2.0),
    series = c("Best baseline ceiling", "Best baseline ceiling", "Best baseline ceiling", "Best baseline ceiling", "Best baseline ceiling", "Best baseline ceiling", "Best baseline ceiling", "Best baseline ceiling", "Best baseline ceiling", "Best baseline ceiling", "Best baseline ceiling", "Best baseline ceiling", "Best baseline ceiling", "Best baseline ceiling", "Best baseline ceiling", "GMUS (elbow)", "GMUS (elbow)", "GMUS (elbow)", "GMUS (elbow)", "GMUS (elbow)", "GMUS (elbow)", "GMUS (elbow)", "GMUS (elbow)", "GMUS (elbow)", "GMUS (elbow)", "GMUS (elbow)", "GMUS (elbow)", "GMUS (elbow)", "GMUS (elbow)", "GMUS (elbow)"),
    TPR_mean = c(0.98, 0.665, 0.43, 0.345, 0.305, 0.84, 0.61, 0.46, 0.41, 0.325, 0.65, 0.59, 0.465, 0.365, 0.295, 0.765, 0.485, 0.375, 0.3, 0.235, 0.645, 0.46, 0.405, 0.34, 0.265, 0.565, 0.485, 0.385, 0.315, 0.255),
    FDR_mean = c(0.054018, 0.120694, 0.071369, 0.202024, 0.281388, 0.134256, 0.109236, 0.117381, 0.204286, 0.235759, 0.107143, 0.063185, 0.1075, 0.245, 0.255208, 0.039524, 0.123333, 0.256518, 0.430817, 0.46501, 0.121131, 0.171548, 0.333968, 0.361483, 0.453898, 0.182083, 0.169286, 0.230575, 0.42004, 0.438362),
    F1_mean = c(0.958498, 0.729491, 0.566726, 0.442204, 0.388169, 0.831486, 0.695019, 0.577397, 0.51542, 0.39913, 0.733929, 0.699756, 0.583636, 0.456064, 0.388235, 0.829939, 0.610496, 0.45282, 0.338929, 0.278256, 0.724161, 0.566494, 0.46988, 0.385928, 0.291602, 0.647309, 0.58991, 0.476207, 0.358299, 0.292885),
    F1_se = c(0.009334, 0.018745, 0.015221, 0.021549, 0.017225, 0.014028, 0.018987, 0.015941, 0.014052, 0.017085, 0.015472, 0.017319, 0.018916, 0.015531, 0.022264, 0.024263, 0.016373, 0.019512, 0.024764, 0.021145, 0.018368, 0.019776, 0.018713, 0.018415, 0.02135, 0.01852, 0.015385, 0.0215, 0.014973, 0.025376),
    conf_any_mean = c(NA, NA, NA, NA, NA, 0.4, 0.025, 0.0, 0.0, 0.0, 0.325, 0.025, 0.05, 0.025, 0.025, NA, NA, NA, NA, NA, 0.45, 0.0, 0.0, 0.0, 0.025, 0.55, 0.075, 0.075, 0.05, 0.05),
    conf_any_se = c(NA, NA, NA, NA, NA, 0.078446, 0.025, 0.0, 0.0, 0.0, 0.075, 0.025, 0.034899, 0.025, 0.025, NA, NA, NA, NA, NA, 0.079663, 0.0, 0.0, 0.0, 0.025, 0.079663, 0.042176, 0.042176, 0.034899, 0.034899),
    size_mean = NA_real_, stringsAsFactors = FALSE)
  ref <- ref[, names(add)]
  gm <- ref[ref$series == "GMUS (elbow)", ]; best <- ref[ref$series == "Best baseline ceiling", ]
}
comb <- rbind(add, gm, best)
rownames(comb) <- NULL
comb <- comb[order(comb$K, comb$delta_obs, comb$series), ]
write.csv(comb, "section3_study4b_combined.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# Paired differences within the add-on run (same data in each replicate)
# ---------------------------------------------------------------------------
se <- function(x) sd(x) / sqrt(length(x))
# pull the columns needed for each comparison
get <- function(a, r, v) {
  d <- raw[raw$aggregation == a & raw$rule == r, c("K", "delta_obs", "rep", v)]
  names(d)[4] <- paste(a, r, v, sep = "|"); d
}
mrg <- Reduce(function(x, y) merge(x, y, by = c("K", "delta_obs", "rep")),
              list(get("SS", "Largest-gap", "F1"), get("SS", "Known s (top-s)", "F1"),
                   get("2D-SS", "Largest-gap", "F1"), get("2D-SS", "Known s (top-s)", "F1"),
                   get("SS", "Known s (top-s)", "conf_any"), get("2D-SS", "Known s (top-s)", "conf_any")))
comps <- list(
  "SS minus 2D-SS, largest-gap (F1)"    = mrg[["SS|Largest-gap|F1"]] - mrg[["2D-SS|Largest-gap|F1"]],
  "SS minus 2D-SS, known s (F1)"        = mrg[["SS|Known s (top-s)|F1"]] - mrg[["2D-SS|Known s (top-s)|F1"]],
  "2D-SS known s minus largest-gap (F1)" = mrg[["2D-SS|Known s (top-s)|F1"]] - mrg[["2D-SS|Largest-gap|F1"]],
  "SS minus 2D-SS, known s (confounder rate)" =
    mrg[["SS|Known s (top-s)|conf_any"]] - mrg[["2D-SS|Known s (top-s)|conf_any"]])
paired <- do.call(rbind, lapply(names(comps), function(nm) {
  d <- cbind(mrg[, c("K", "delta_obs")], diff = comps[[nm]])
  if (grepl("confounder", nm)) d <- d[d$K > 0, ]
  a <- aggregate(diff ~ K + delta_obs, d, function(x) c(mean = mean(x), se = se(x), n = length(x)))
  data.frame(comparison = nm, K = a$K, delta_obs = a$delta_obs,
             diff_mean = a$diff[, "mean"], diff_se = a$diff[, "se"], n = a$diff[, "n"])
}))
write.csv(paired, "section3_study4b_paired.csv", row.names = FALSE)
print(paired[paired$delta_obs %in% c(0, 1), ], digits = 3, row.names = FALSE)

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
lv <- c("2D-SS: largest-gap", "2D-SS: known s", "SS: largest-gap", "SS: known s",
        "GMUS (elbow)", "Best baseline ceiling")
cl <- c("2D-SS: largest-gap" = "#0072B2", "2D-SS: known s" = "#0072B2",
        "SS: largest-gap" = "#D55E00", "SS: known s" = "#D55E00",
        "GMUS (elbow)" = "#E69F00", "Best baseline ceiling" = "#999999")
lt <- c("2D-SS: largest-gap" = "solid", "2D-SS: known s" = "dashed",
        "SS: largest-gap" = "solid", "SS: known s" = "dashed",
        "GMUS (elbow)" = "solid", "Best baseline ceiling" = "dotted")
sh <- c("2D-SS: largest-gap" = 8, "2D-SS: known s" = 15, "SS: largest-gap" = 18, "SS: known s" = 17,
        "GMUS (elbow)" = 4, "Best baseline ceiling" = 1)
comb$series <- factor(comb$series, levels = lv)
comb$design <- factor(ifelse(comb$K == 0, "No confounder (control)",
                             ifelse(comb$K == 1, "1 confounder", "5 confounders")),
                      levels = c("1 confounder", "5 confounders", "No confounder (control)"))
base <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"))

p1 <- ggplot(comb, aes(delta_obs, F1_mean, colour = series, linetype = series, shape = series)) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(F1_mean - F1_se, 0), ymax = pmin(F1_mean + F1_se, 1)),
                width = 0.04, linewidth = 0.3) +
  facet_wrap(~ design) +
  scale_colour_manual(values = cl) + scale_linetype_manual(values = lt) + scale_shape_manual(values = sh) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2)) +
  labs(x = expression(delta[obs]), y = expression(F[1]), colour = NULL, linetype = NULL, shape = NULL) +
  base + guides(colour = guide_legend(nrow = 2))
ggsave("plot_study4b_f1.pdf", p1, width = 11, height = 5.2)

cf <- comb[comb$K > 0 & is.finite(comb$conf_any_mean), ]
p2 <- ggplot(cf, aes(delta_obs, conf_any_mean, colour = series, linetype = series, shape = series)) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(conf_any_mean - conf_any_se, 0), ymax = pmin(conf_any_mean + conf_any_se, 1)),
                width = 0.04, linewidth = 0.3) +
  facet_wrap(~ design) +
  scale_colour_manual(values = cl) + scale_linetype_manual(values = lt) + scale_shape_manual(values = sh) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2)) +
  labs(x = expression(delta[obs]), y = "P(any confounder selected)", colour = NULL, linetype = NULL, shape = NULL) +
  base + guides(colour = guide_legend(nrow = 2))
ggsave("plot_study4b_confounder.pdf", p2, width = 8, height = 5)

cat("Saved: section3_study4b_combined.csv, section3_study4b_paired.csv,",
    "plot_study4b_f1.pdf, plot_study4b_confounder.pdf\n")