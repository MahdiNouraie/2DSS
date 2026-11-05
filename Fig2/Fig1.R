suppressWarnings(options(warn = -1))

# ---------------- packages ----------------
required_pkgs <- c("glmnet", "ggplot2", "scales")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs, dependencies = TRUE)
invisible(lapply(required_pkgs, library, character.only = TRUE))

set.seed(26)

# ---------------- settings ----------------
n  <- 100
p  <- 1000
n1 <- floor(n / 2)
s  <- 5
K  <- 10                        # number of confounders
rho_irr <- 0.1
lambda  <- sqrt(log(p) / n1)    # fixed, not cross-validated
B <- 100
delta_grid <- c(0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3)

n_weak   <- 4      # number of weak relevant features (at least 2)
m_weak   <- 4      # weak   beta = m_weak   * lambda
m_strong <- 15     # strong beta = m_strong * lambda
x_conf   <- 0.4    # confounder covariance with each weak feature
B_target <- -0.3   # target B = a'b (negative => repair at large nu)

# ---------------- design ----------------
n_strong <- s - n_weak
eps_conf <- (n_weak * x_conf * m_weak - B_target) / (n_strong * m_strong)
a        <- c(rep(x_conf, n_weak), rep(-eps_conf, n_strong))
beta_rel <- c(rep(m_weak, n_weak), rep(m_strong, n_strong)) * lambda
stopifnot(sum(a^2) < 1)
cat("A = a'z =", sum(a), " B = a'b =", sum(a * beta_rel / lambda),
    " norm(a) =", sqrt(sum(a^2)), "\n")
stopifnot(sum(a) > 1)           # Lasso IC violated at delta = 0

# ---------------- data ----------------
XS <- matrix(rnorm(n * s), n, s)
Xc <- sapply(1:K, function(k) as.numeric(XS %*% a) + sqrt(1 - sum(a^2)) * rnorm(n))
f0 <- rnorm(n)
Xo <- sapply(1:(p - s - K), function(k) sqrt(rho_irr) * f0 + sqrt(1 - rho_irr) * rnorm(n))
X  <- scale(cbind(XS, Xc, Xo))
beta <- c(beta_rel, rep(0, p - s))
y <- as.numeric(X %*% beta + rnorm(n))

# ---------------- sub-sampled, jittered Lasso fits ----------------
sel_freq <- matrix(0, nrow = length(delta_grid), ncol = p)
for (i in seq_along(delta_grid)) {
  sel_count <- numeric(p)
  for (b in 1:B) {
    idx_b <- sample(seq_len(n), size = n1, replace = FALSE)
    Xj <- X[idx_b, , drop = FALSE] +
      matrix(rnorm(n1 * p, sd = delta_grid[i]), n1, p)
    fit <- glmnet(Xj, y[idx_b], alpha = 1, lambda = lambda)
    sel_count <- sel_count + (as.numeric(coef(fit))[-1] != 0)
  }
  sel_freq[i, ] <- sel_count / B
}

# ---------------- groups ----------------
grp <- rep("Other irrelevant", p)
grp[1:n_weak] <- "Weak relevant"
grp[(n_weak + 1):s] <- "Strong relevant"
grp[(s + 1):(s + K)] <- "Confounder"
grp_levels <- c("Strong relevant", "Weak relevant", "Confounder", "Other irrelevant")
grp <- factor(grp, levels = grp_levels)

df <- data.frame(
  delta = rep(delta_grid, times = p),
  var   = rep(1:p, each = length(delta_grid)),
  freq  = as.vector(sel_freq),
  group = rep(grp, each = length(delta_grid))
)

# ---------------- save results as CSV ----------------
write.csv(df, "fig1_ic_violated_frequencies.csv", row.names = FALSE)

summ <- do.call(rbind, lapply(grp_levels, function(g) {
  idx <- which(grp == g)
  F <- sel_freq[, idx, drop = FALSE]
  data.frame(group = g, delta = delta_grid, n_features = length(idx),
             mean = rowMeans(F), min = apply(F, 1, min), max = apply(F, 1, max))
}))
write.csv(summ, "fig1_ic_violated_summary.csv", row.names = FALSE)

b <- beta_rel / lambda
z <- sign(beta_rel)
cnu <- sapply(delta_grid, function(d) abs(sum(a * (z + d^2 * b))) / (1 + d^2))
write.csv(data.frame(delta = delta_grid, c_nu = cnu), "fig1_ic_violated_cnu.csv", row.names = FALSE)

write.csv(data.frame(n = n, p = p, n1 = n1, K = K, rho_irr = rho_irr,
                     lambda = lambda, B = B, n_weak = n_weak, m_weak = m_weak,
                     m_strong = m_strong, x_conf = x_conf,
                     A = sum(a * z), B_stat = sum(a * b), norm_a = sqrt(sum(a^2))),
          "fig1_ic_violated_design.csv", row.names = FALSE)

# ---------------- plot ----------------
draw_order <- c("Other irrelevant", "Confounder", "Weak relevant", "Strong relevant")
df <- df[order(match(as.character(df$group), draw_order)), ]
df$var <- factor(df$var, levels = unique(df$var))

p1 <- ggplot(df, aes(x = delta, y = freq, group = var, color = group,
                     linewidth = group, alpha = group)) +
  geom_line() +
  scale_color_manual(values = c("Strong relevant" = "red", "Weak relevant" = "orange",
                                "Confounder" = "blue", "Other irrelevant" = "black")) +
  scale_linewidth_manual(values = c("Strong relevant" = 1.0, "Weak relevant" = 1.0,
                                    "Confounder" = 0.7, "Other irrelevant" = 0.3)) +
  scale_alpha_manual(values = c("Strong relevant" = 1, "Weak relevant" = 1,
                                "Confounder" = 1, "Other irrelevant" = 0.35)) +
  ylim(0, 1) +
  scale_x_continuous(trans = scales::compose_trans("sqrt", "reverse"),
                     breaks = c(0, 0.1, 0.25, 0.5, 1, 2, 3)) +
  theme_bw(base_size = 16) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "right") +
  labs(x = expression(delta), y = "Selection frequency")

print(p1)
ggsave("fig1_ic_violated.pdf", plot = p1, width = 8, height = 5, device = "pdf")
cat("Saved: fig1_ic_violated.pdf and four CSV files (frequencies, summary, cnu, design)\n")