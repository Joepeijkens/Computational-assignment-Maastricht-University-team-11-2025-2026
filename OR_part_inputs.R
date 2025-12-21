# OR_part_inputs.R
# Computes OR input-model quantities FROM ScanRecords.csv using lecture-consistent bootstrap:
# - Type 1 duration: Normal assumption -> parametric bootstrap
# - Type 2 duration: unknown -> nonparametric bootstrap (EDF)
# - Type 1 arrivals: Poisson per day -> lambda + bootstrap CI + implied Exp interarrival mean
# - Type 2 interarrival times: unknown -> empirical distribution + nonparametric bootstrap stats
#
# Needs: ScanRecords.csv in working directory

# -----------------------------
# 0) Settings
# -----------------------------
CSV_PATH <- "ScanRecords.csv"

WORK_START_H <- 8
WORK_END_H   <- 17
DAY_LEN_MIN  <- (WORK_END_H - WORK_START_H) * 60  # 540 minutes

B <- 5000
set.seed(123)

# -----------------------------
# 1) Load + preprocess (base R)
# -----------------------------
df <- read.csv(CSV_PATH, stringsAsFactors = FALSE)
df$Date <- as.Date(df$Date)
df$Day  <- df$Date

# Convert to minutes since 08:00
df$CallMin     <- (df$Time - WORK_START_H) * 60
df$DurationMin <- df$Duration * 60

df1 <- df[df$PatientType == "Type 1", ]
df2 <- df[df$PatientType == "Type 2", ]

dur1 <- df1$DurationMin
dur2 <- df2$DurationMin

# Daily arrival counts by type
days <- sort(unique(df$Day))
n_days <- length(days)
c1 <- integer(n_days)
c2 <- integer(n_days)
for (i in seq_along(days)) {
  d <- days[i]
  c1[i] <- sum(df$Day == d & df$PatientType == "Type 1")
  c2[i] <- sum(df$Day == d & df$PatientType == "Type 2")
}

# -----------------------------
# 2) Bootstrap helpers (lecture algorithm)
# -----------------------------
# Nonparametric bootstrap: resample from EDF (sample with replacement)
boot_np <- function(x, statfun, B = 5000) {
  n <- length(x)
  out <- numeric(B)
  for (b in 1:B) {
    xb <- sample(x, n, replace = TRUE)
    out[b] <- statfun(xb)
  }
  out
}

# Parametric bootstrap (Normal, truncated > 0): resample from N(mu_hat, sd_hat)
boot_param_truncnorm <- function(n, mu, sd, statfun, B = 5000) {
  out <- numeric(B)
  for (b in 1:B) {
    x <- rnorm(n, mu, sd)
    while (any(x <= 0)) {
      idx <- which(x <= 0)
      x[idx] <- rnorm(length(idx), mu, sd)
    }
    out[b] <- statfun(x)
  }
  out
}

ci95 <- function(samples) {
  quantile(samples, probs = c(0.025, 0.975), names = FALSE)
}

# -----------------------------
# 3) TYPE 1 duration (Normal) -> mean, sd + parametric bootstrap
# -----------------------------
mu1_hat <- mean(dur1)
sd1_hat <- sd(dur1)

boot_mu1 <- boot_param_truncnorm(length(dur1), mu1_hat, sd1_hat, mean, B)
boot_sd1 <- boot_param_truncnorm(length(dur1), mu1_hat, sd1_hat, sd,   B)

mu1_ci <- ci95(boot_mu1)
sd1_ci <- ci95(boot_sd1)

# -----------------------------
# 4) TYPE 2 duration (unknown) -> nonparametric bootstrap
# -----------------------------
mu2_hat <- mean(dur2)
sd2_hat <- sd(dur2)

boot_mu2 <- boot_np(dur2, mean, B)
boot_sd2 <- boot_np(dur2, sd,   B)

mu2_ci <- ci95(boot_mu2)
sd2_ci <- ci95(boot_sd2)

# Optional (for DES only): choose a parametric family by AIC among Gamma/Lognormal/Weibull
# This does NOT replace NP bootstrap inference; it’s only for random variate generation.
loglik_gamma <- function(x, shape, scale) sum(dgamma(x, shape = shape, scale = scale, log = TRUE))
loglik_lnorm <- function(x, meanlog, sdlog) sum(dlnorm(x, meanlog = meanlog, sdlog = sdlog, log = TRUE))
loglik_weib  <- function(x, shape, scale) sum(dweibull(x, shape = shape, scale = scale, log = TRUE))
AIC_val <- function(loglik, k) 2*k - 2*loglik

# MOM starting values (simple and stable)
gamma_shape0 <- (mu2_hat^2) / (sd2_hat^2)
gamma_scale0 <- (sd2_hat^2) / mu2_hat

# Estimate parameters (simple MLE via optim)
negll_gamma <- function(par) -loglik_gamma(dur2, shape = exp(par[1]), scale = exp(par[2]))
negll_lnorm <- function(par) -loglik_lnorm(dur2, meanlog = par[1], sdlog = exp(par[2]))
negll_weib  <- function(par) -loglik_weib(dur2, shape = exp(par[1]), scale = exp(par[2]))

fit_g <- optim(c(log(gamma_shape0), log(gamma_scale0)), negll_gamma)
fit_l <- optim(c(mean(log(dur2)), log(sd(log(dur2)))), negll_lnorm)
fit_w <- optim(c(log(2), log(mu2_hat)), negll_weib)

g_shape <- exp(fit_g$par[1]); g_scale <- exp(fit_g$par[2])
l_meanlog <- fit_l$par[1];    l_sdlog <- exp(fit_l$par[2])
w_shape <- exp(fit_w$par[1]); w_scale <- exp(fit_w$par[2])

aic_g <- AIC_val(-fit_g$value, 2)
aic_l <- AIC_val(-fit_l$value, 2)
aic_w <- AIC_val(-fit_w$value, 2)

type2_duration_parametric_choice <- c(Gamma = aic_g, Lognormal = aic_l, Weibull = aic_w)
best_type2_duration_family <- names(which.min(type2_duration_parametric_choice))

# -----------------------------
# 5) TYPE 1 arrivals: Poisson per day -> lambda + bootstrap CI + implied Exp interarrival mean
# -----------------------------
lambda1_hat <- mean(c1)
boot_lam1 <- boot_np(c1, mean, B)   # bootstrap over days (EDF)
lambda1_ci <- ci95(boot_lam1)

rate1_per_min <- lambda1_hat / DAY_LEN_MIN  # constant rate during 08:00-17:00
mean_interarrival1_min <- 1 / rate1_per_min # Exponential mean (minutes)

# -----------------------------
# 6) TYPE 2 interarrival times: empirical distribution + NP bootstrap stats
# -----------------------------
type2_inter <- c()
for (d in days) {
  times <- df2$CallMin[df2$Day == d]
  if (length(times) >= 2) {
    times <- sort(times)
    type2_inter <- c(type2_inter, diff(times))
  }
}

i2_mean <- mean(type2_inter)
i2_sd   <- sd(type2_inter)
i2_p50  <- as.numeric(quantile(type2_inter, 0.50))
i2_p90  <- as.numeric(quantile(type2_inter, 0.90))
i2_p95  <- as.numeric(quantile(type2_inter, 0.95))

boot_i2_mean <- boot_np(type2_inter, mean, B)
boot_i2_sd   <- boot_np(type2_inter, sd,   B)

i2_mean_ci <- ci95(boot_i2_mean)
i2_sd_ci   <- ci95(boot_i2_sd)

# -----------------------------
# 7) Print the four OR items (plus CIs)
# -----------------------------
cat("\n==================== OR INPUTS (from ScanRecords.csv) ====================\n")

cat("\n[1] Patient 1 duration (Normal)\n")
cat(sprintf("Mean (min): %.3f   95%% param-bootstrap CI: [%.3f, %.3f]\n", mu1_hat, mu1_ci[1], mu1_ci[2]))
cat(sprintf("SD   (min): %.3f   95%% param-bootstrap CI: [%.3f, %.3f]\n", sd1_hat, sd1_ci[1], sd1_ci[2]))

cat("\n[2] Patient 2 duration (unknown -> EDF for inference)\n")
cat(sprintf("Mean (min): %.3f   95%% NP-bootstrap CI: [%.3f, %.3f]\n", mu2_hat, mu2_ci[1], mu2_ci[2]))
cat(sprintf("SD   (min): %.3f   95%% NP-bootstrap CI: [%.3f, %.3f]\n", sd2_hat, sd2_ci[1], sd2_ci[2]))
cat("Distribution type (for inference): Empirical/EDF (nonparametric)\n")
cat(sprintf("Optional DES parametric family by AIC (Gamma/Lognormal/Weibull): best = %s\n", best_type2_duration_family))

cat("\n[3] Patient 1 arrivals (Poisson per day)\n")
cat(sprintf("lambda (calls/day): %.3f   95%% bootstrap CI: [%.3f, %.3f]\n", lambda1_hat, lambda1_ci[1], lambda1_ci[2]))
cat(sprintf("Implied Exp mean interarrival during 08:00-17:00: %.3f minutes\n", mean_interarrival1_min))

cat("\n[4] Patient 2 interarrival times\n")
cat("Distribution type: Empirical/EDF (nonparametric)\n")
cat(sprintf("Mean (min): %.3f   95%% NP-bootstrap CI: [%.3f, %.3f]\n", i2_mean, i2_mean_ci[1], i2_mean_ci[2]))
cat(sprintf("SD   (min): %.3f   95%% NP-bootstrap CI: [%.3f, %.3f]\n", i2_sd,   i2_sd_ci[1],   i2_sd_ci[2]))
cat(sprintf("Quantiles (min): p50=%.3f, p90=%.3f, p95=%.3f\n", i2_p50, i2_p90, i2_p95))

cat("\n=========================================================================\n\n")

# -----------------------------
# 8) Export a compact CSV for your report
# -----------------------------
out <- data.frame(
  item = c(
    "Type1_duration_mean_min","Type1_duration_sd_min",
    "Type2_duration_mean_min","Type2_duration_sd_min",
    "Type2_duration_inference_dist","Type2_duration_best_parametric_family_AIC",
    "Type1_lambda_calls_per_day","Type1_mean_interarrival_min",
    "Type2_interarrival_dist","Type2_interarrival_mean_min","Type2_interarrival_sd_min",
    "Type2_interarrival_p50_min","Type2_interarrival_p90_min","Type2_interarrival_p95_min"
  ),
  value = c(
    mu1_hat, sd1_hat,
    mu2_hat, sd2_hat,
    "Empirical/EDF (nonparametric)", best_type2_duration_family,
    lambda1_hat, mean_interarrival1_min,
    "Empirical/EDF (nonparametric)", i2_mean, i2_sd,
    i2_p50, i2_p90, i2_p95
  )
)
write.csv(out, "OR_inputs_summary.csv", row.names = FALSE)
