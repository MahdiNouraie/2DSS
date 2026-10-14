# Sensitivity of 2D-SS to its free inputs (grid and weights), confounder design.
# Design: as Study 2 (Sigma_SS = I_5, beta_S = (-4,-2,.5,.6,.7), confounder
# covariance a = (1/8,1/4,1/2,1/2,1/2), p = 1000, delta_obs = 0), n = 200, K in {1,5}.
#
# All grid/weight configurations are formed from ONE set of selection frequencies per
# replicate, computed on a fine grid delta = 0, 0.125, ..., 3 (B = 50 half-sample fits
# per level). The configurations therefore see the same data and the same fits.
# Reference: Stability Selection with B = 250 half-sample fits at delta = 0.
#
# Families of configurations (m levels equally spaced from 0 to delta_max unless stated)
#   upper : delta_max in {0.5, 1, 2, 3}, m = 5, uniform weights
#   levels: m in {3, 5, 9}, delta_max = 1, uniform weights
#   weights: w_delta ~ exp(-kappa * delta / delta_max), kappa in {-2, 0, 2, 4},
#            m = 5, delta_max in {1, 2}  (kappa < 0 favours large delta)
#   noSS  : delta in {0.25, 0.5, 0.75, 1} (no delta = 0 slice), uniform
# Rules: largest-gap, and known s (top-5 frequencies).
# Outputs: sens_raw.csv, sens_summary.csv, plot_sens_f1.pdf,
#          plot_sens_confounder.pdf, plot_sens_margin.pdf

suppressWarnings(options(warn = -1))
for (pk in c("glmnet", "parallel", "ggplot2"))
  if (!requireNamespace(pk, quietly = TRUE)) install.packages(pk)
library(glmnet); library(parallel); library(ggplot2)

quick_test <- FALSE             # TRUE: 4 replicates, for a speed check
RNGkind("L'Ecuyer-CMRG"); set.seed(2027)
n_cores <- if (.Platform$OS.type == "windows") 1 else max(1, min(7, detectCores() - 1))

p <- 1000; s <- 5
beta_S <- c(-4, -2, 0.5, 0.6, 0.7)
a_conf <- c(1/8, 1/4, 1/2, 1/2, 1/2)
idio_sd <- sqrt(1 - sum(a_conf^2))
n <- 200; n1 <- n / 2
lam <- sqrt(log(p) / n1)
B <- 50; B_SS <- 250
K_grid <- c(1, 5)
n_rep <- if (quick_test) 4 else 100
fine <- seq(0, 3, by = 0.125)

## ---- configurations ----
mk <- function(family, label, grid, kappa = 0) {
  w <- exp(-kappa * grid / max(grid, 1e-12)); w <- w / sum(w)
  out <- data.frame(family = family, label = label, stringsAsFactors = FALSE)
  out$grid <- I(list(grid)); out$w <- I(list(w))
  out
}
cfg <- do.call(rbind, c(
  lapply(c(0.5, 1, 2, 3), function(d) mk("Upper end (m = 5)", as.character(d), seq(0, d, length.out = 5))),
  lapply(c(3, 5, 9),      function(m) mk("Number of levels (max 1)", as.character(m), seq(0, 1, length.out = m))),
  lapply(c(-2, 0, 2, 4),  function(k) mk("Weights (max 1)", as.character(k), seq(0, 1, length.out = 5), k)),
  lapply(c(-2, 0, 2, 4),  function(k) mk("Weights (max 2)", as.character(k), seq(0, 2, length.out = 5), k)),
  list(mk("No delta = 0 slice", "0.25-1", c(0.25, 0.5, 0.75, 1)))
))
cfg$id <- seq_len(nrow(cfg))
cfg$idx <- lapply(cfg$grid, function(g) match(round(g, 6), round(fine, 6)))
stopifnot(!anyNA(unlist(cfg$idx)))

## ---- data and selectors (as Study 2) ----
generate_data <- function(K) {
  X_S <- matrix(rnorm(n * s), n, s)
  base <- as.numeric(X_S %*% a_conf)
  X_C <- sapply(seq_len(K), function(k) base + idio_sd * rnorm(n))
  X <- cbind(X_S, X_C, matrix(rnorm(n * (p - s - K)), n, p - s - K))
  list(X = X, y = as.numeric(X_S %*% beta_S + rnorm(n)))
}
level_freq <- function(X, y, delta, B_use) {
  cnt <- numeric(p)
  for (b in seq_len(B_use)) {
    idx <- sample(n, n1)
    Xj <- X[idx, , drop = FALSE]
    if (delta > 0) Xj <- Xj + matrix(rnorm(n1 * p, sd = delta), n1, p)
    fit <- glmnet(scale(Xj), y[idx], alpha = 1, lambda = lam, standardize = TRUE)
    cnt <- cnt + (as.numeric(coef(fit))[-1] != 0)
  }
  cnt / B_use
}
gap_select <- function(f) {
  o <- order(f, decreasing = TRUE); sf <- f[o]
  k <- which.max(sf[-p] - sf[-1]); which(f >= 0.5 * (sf[k] + sf[k + 1]))
}
metrics <- function(S, K) {
  tp <- sum(S %in% 1:s); conf <- (s + 1):(s + K)
  data.frame(size = length(S), TPR = tp / s, FDR = if (length(S)) (length(S) - tp) / length(S) else 0,
             F1 = 2 * tp / (length(S) + s), exact = as.numeric(setequal(S, 1:s)),
             conf_any = as.numeric(any(S %in% conf)))
}

one_rep <- function(r, K) {
  d <- generate_data(K)
  FM <- sapply(fine, function(dl) level_freq(d$X, d$y, dl, B))     # p x length(fine)
  f_ss <- level_freq(d$X, d$y, 0, B_SS)
  conf <- (s + 1):(s + K)
  rows <- list()
  add <- function(id, f) {
    top <- order(-f, runif(p))[1:s]
    margin <- min(f[1:s]) - max(f[conf])
    for (ru in c("gap", "known_s")) {
      S <- if (ru == "gap") gap_select(f) else top
      rows[[length(rows) + 1]] <<- cbind(data.frame(rep = r, K = K, id = id, rule = ru, margin = margin),
                                         metrics(S, K))
    }
  }
  add(0, f_ss)
  for (i in seq_len(nrow(cfg))) {
    f <- as.numeric(FM[, cfg$idx[[i]], drop = FALSE] %*% cfg$w[[i]])
    add(cfg$id[i], f)
  }
  do.call(rbind, rows)
}

t0 <- Sys.time()
res <- lapply(K_grid, function(K) {
  out <- mclapply(seq_len(n_rep), one_rep, K = K, mc.cores = n_cores)
  bad <- vapply(out, inherits, logical(1), what = "try-error")
  if (any(bad)) stop("worker failed: ", as.character(out[[which(bad)[1]]]))
  do.call(rbind, out)
})
raw <- do.call(rbind, res)
cat("elapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "\n")
raw$family <- c("SS (B = 250)", cfg$family)[match(raw$id, c(0, cfg$id))]
raw$label  <- c("SS", cfg$label)[match(raw$id, c(0, cfg$id))]
write.csv(raw, "sens_raw.csv", row.names = FALSE)

se <- function(v) sd(v) / sqrt(length(v))
agg <- aggregate(cbind(size, TPR, FDR, F1, exact, conf_any, margin) ~ K + id + family + label + rule,
                 data = raw, FUN = mean)
agg_se <- aggregate(cbind(F1, conf_any, margin) ~ K + id + rule, data = raw, FUN = se)
names(agg_se)[4:6] <- c("F1_se", "conf_any_se", "margin_se")
summ <- merge(agg, agg_se, by = c("K", "id", "rule"))
write.csv(summ[order(summ$K, summ$rule, summ$id), ], "sens_summary.csv", row.names = FALSE)
print(summ[order(summ$K, summ$rule, summ$id), c("K", "family", "label", "rule", "F1", "conf_any", "margin")], digits = 3)

## ---- plots (Okabe-Ito) ----
cols <- c(gap = "#0072B2", known_s = "#D55E00")
rl <- c(gap = "Largest-gap", known_s = "Known s")
summ$rule_lab <- rl[summ$rule]
summ$Klab <- ifelse(summ$K == 1, "1 confounder", "5 confounders")
ref <- summ[summ$id == 0, ]; dat <- summ[summ$id != 0, ]
fam_levels <- unique(dat$family[order(dat$id)])
dat$family <- factor(dat$family, levels = fam_levels)
ord <- dat[!duplicated(dat$id), ]; ord <- ord[order(ord$id), ]
dat$xlab <- factor(paste(dat$family, dat$label, sep = "|"),
                   levels = paste(ord$family, ord$label, sep = "|"))
mkplot <- function(y, se_col, ylab, file) {
  ref_y <- ref[, c("Klab", "rule_lab", y)]      # no facet column: drawn in every panel
  g <- ggplot(dat, aes(xlab, .data[[y]], colour = rule_lab, shape = rule_lab, group = rule_lab)) +
    geom_hline(data = ref_y, aes(yintercept = .data[[y]], colour = rule_lab), linetype = "dashed") +
    geom_line() + geom_point(size = 2.2) +
    {if (!is.null(se_col)) geom_errorbar(aes(ymin = .data[[y]] - .data[[se_col]], ymax = .data[[y]] + .data[[se_col]]), width = .15)} +
    facet_grid(Klab ~ family, scales = "free_x", space = "free_x",
               labeller = labeller(family = label_wrap_gen(10))) +
    scale_x_discrete(labels = function(x) sub(".*\\|", "", x)) +
    scale_colour_manual(values = setNames(cols, rl)) +
    labs(x = "Grid setting (delta_max, number of levels, kappa, or fixed levels)", y = ylab, colour = NULL, shape = NULL) +
    theme_bw(base_size = 9) + theme(legend.position = "bottom", strip.text = element_text(size = 7))
  ggsave(file, g, width = 12, height = 4.8)
}
tryCatch({
  mkplot("F1", "F1_se", "F1 (dashed: Stability Selection, B = 250)", "plot_sens_f1.pdf")
  mkplot("conf_any", "conf_any_se", "Confounder inclusion rate", "plot_sens_confounder.pdf")
  mkplot("margin", "margin_se", "min relevant freq - max confounder freq", "plot_sens_margin.pdf")
}, error = function(e) message("Plotting failed (CSVs are saved): ", conditionMessage(e)))
cat("Done.\n")