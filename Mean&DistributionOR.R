# OR_inputs_bootstrap.R
# Purpose: Provide OR input-model answers WITH bootstraps (as the assignment expects).
# Input:   ScanRecords.csv in the working directory
# Output:  prints results + saves OR_bootstrap_summary.csv

# -----------------------------
# 0) Settings
# -----------------------------
CSV_PATH <- "ScanRecords.csv"

WORK_START_H <- 8
WORK_END_H   <- 17
DAY_LEN_MIN  <- (WORK_END_H - WORK_START_H) * 60  # 540

B <- 5000
set.seed(123)

# -----------------------------
# 1) Load data (base R)
# -----------------------------
df <- read.csv(CSV_PATH, stringsAsFactors = FALSE)

# parse date (assumes yyyy-mm-dd or similar)
df$Date <- as.Date(df$Date)
df$Day  <- df$Date

# convert to minutes
df$CallMin     <- (df$Time - WORK_START_H) * 60
df$DurationMin <- df$Duration * 60

# split by type
df1 <- df[df$PatientType == "Type 1", ]
df2 <- df[df$PatientType == "Type 2", ]

dur1 <- df1$DurationMin
dur2 <- df2$DurationMin

# daily counts by type
# (tabulate counts per day per type)
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
# 2) Bootstrap helpers (simple)
# -----------------------------
boot_np <- function(x, statfun, B = 5000) {
  n <- length(x)
  out <- numeric(B)
  for (b in 1:B) {
    xb <- sample(x, n, replace = TRUE)
    out[b] <- statfun(xb)
  }
  return(out)
}

boot_param_truncnorm <- function(n, mu, sd, statfun, B = 5000) {
  out <- numeric(B)
  for (b in 1:B) {
    x <- rnorm(n, mu, sd)
    # truncate at >0 by resampling negatives
    while (any(x <= 0)) {
      idx <- which(x <= 0)
      x[idx] <- rnorm(length(idx), mu, sd)
    }
    out[b] <- statfun(x)
  }
  return(out)
}

ci95 <- function(samples) {
  quantile(samples, probs = c(0.025, 0.975), names = FALSE)
}

# -----------------------------
# 3) TYPE 1 duration: Normal + parametric bootstrap
# -----------------------------
mu1_hat <- mean(dur1)
sd1_hat <- sd(dur1)

# Parametric bootstrap for mean and for SD (and optionally 90% quantile)
boot_mu1 <- boot_param_truncnorm(length(dur1), mu1_hat, sd1_hat, mean, B)
boot_sd1 <- boot_param_truncnorm(length(dur1), mu1_hat, sd1_hat, sd,   B)
boot_q1  <- boot_param_truncnorm(length(dur1), mu1_hat, sd1_hat,
                                 function(x) as.numeric(quantile(x, 0.90)),
                                 B)

mu1_ci <- ci95(boot_mu1)
sd1_ci <- ci95(boot_sd1)
q1_ci  <- ci95(boot_q1)

# -----------------------------
# 4) TYPE 2 duration: Nonparametric bootstrap
# -----------------------------
mu2_hat <- mean(dur2)
sd2_hat <- sd(dur2)
q2_hat  <- as.numeric(quantile(dur2, 0.90))

boot_mu2 <- boot_np(dur2, mean, B)
boot_sd2 <- boot_np(dur2, sd,   B)
boot_q2  <- boot_np(dur2, function(x) as.numeric(quantile(x, 0.90)), B)

mu2_ci <- ci95(boot_mu2)
sd2_ci <- ci95(boot_sd2)
q2_ci  <- ci95(boot_q2)

# -----------------------------
# 5) TYPE 1 arrivals: Poisson per day
# -----------------------------
lambda1_hat <- mean(c1)
boot_lam1 <- boot_np(c1, mean, B)     # bootstrap over days
lambda1_ci <- ci95(boot_lam1)

# implied exponential interarrival mean during working hours
rate1_per_min <- lambda1_hat / DAY_LEN_MIN
mean_interarrival1_min <- 1 / rate1_per_min

# -----------------------------
# 6) TYPE 2 interarrival times: empirical + nonparametric bootstrap
# -----------------------------
# compute within-day interarrival times from observed call times
type2_inter <- c()
for (d in days) {
  times <- df2$CallMin[df2$Day == d]
  if (length(times) >= 2) {
    times <- sort(times)
    type2_inter <- c(type2_inter, diff(times))
  }
}

# basic stats
i2_mean <- mean(type2_inter)
i2_sd   <- sd(type2_inter)
i2_q90  <- as.numeric(quantile(type2_inter, 0.90))

# nonparametric bootstrap on interarrivals
boot_i2_mean <- boot_np(type2_inter, mean, B)
boot_i2_sd   <- boot_np(type2_inter, sd,   B)
boot_i2_q90  <- boot_np(type2_inter, function(x) as.numeric(quantile(x, 0.90)), B)

i2_mean_ci <- ci95(boot_i2_mean)
i2_sd_ci   <- ci95(boot_i2_sd)
i2_q90_ci  <- ci95(boot_i2_q90)

# Also check whether Type 2 daily counts look Poisson-ish (mean vs variance)
c2_mean <- mean(c2)
c2_var  <- var(c2)

# -----------------------------
# 7) OPTIONAL: Gamma fit for Type 2 durations (DES input only)
# -----------------------------
# This is NOT needed to answer “unknown distribution” in Part 1;
# it is only useful if you want a parametric generator for simulation.
# (Method of moments as a simple, readable fit.)
gamma_shape_mom <- (mu2_hat^2) / (sd2_hat^2)
gamma_scale_mom <- (sd2_hat^2) / mu2_hat

# -----------------------------
# 8) Print answers (OR-style)
# -----------------------------
cat("\n==================== OR INPUTS (WITH BOOTSTRAPS) ====================\n")

cat("\n[Type 1 duration]  Normal assumption (truncated >0 in simulation)\n")
cat(sprintf("Mean (min): %.3f   95%% param-bootstrap CI: [%.3f, %.3f]\n", mu1_hat, mu1_ci[1], mu1_ci[2]))
cat(sprintf("SD   (min): %.3f   95%% param-bootstrap CI: [%.3f, %.3f]\n", sd1_hat, sd1_ci[1], sd1_ci[2]))
cat(sprintf("90%% quantile (min): %.3f   95%% param-bootstrap CI: [%.3f, %.3f]\n",
            as.numeric(quantile(dur1, 0.90)), q1_ci[1], q1_ci[2]))

cat("\n[Type 2 duration]  Distribution unknown -> Nonparametric bootstrap\n")
cat(sprintf("Mean (min): %.3f   95%% NP-bootstrap CI: [%.3f, %.3f]\n", mu2_hat, mu2_ci[1], mu2_ci[2]))
cat(sprintf("SD   (min): %.3f   95%% NP-bootstrap CI: [%.3f, %.3f]\n", sd2_hat, sd2_ci[1], sd2_ci[2]))
cat(sprintf("90%% quantile (min): %.3f   95%% NP-bootstrap CI: [%.3f, %.3f]\n", q2_hat, q2_ci[1], q2_ci[2]))

cat("\n[Type 1 arrivals]  Daily calls Poisson(lambda1)\n")
cat(sprintf("lambda1 (calls/day): %.3f   95%% bootstrap CI: [%.3f, %.3f]\n", lambda1_hat, lambda1_ci[1], lambda1_ci[2]))
cat(sprintf("Implied Exp interarrival mean during 08:00-17:00: %.3f minutes\n", mean_interarrival1_min))

cat("\n[Type 2 interarrival times]  Distribution unknown -> Empirical + NP-bootstrap\n")
cat(sprintf("Interarrival mean (min): %.3f   95%% CI: [%.3f, %.3f]\n", i2_mean, i2_mean_ci[1], i2_mean_ci[2]))
cat(sprintf("Interarrival SD   (min): %.3f   95%% CI: [%.3f, %.3f]\n", i2_sd, i2_sd_ci[1], i2_sd_ci[2]))
cat(sprintf("Interarrival 90%% q (min): %.3f   95%% CI: [%.3f, %.3f]\n", i2_q90, i2_q90_ci[1], i2_q90_ci[2]))
cat(sprintf("Type 2 daily counts: mean = %.3f, var = %.3f (var>mean suggests non-Poisson)\n", c2_mean, c2_var))

cat("\n[Optional for DES only] Gamma method-of-moments fit for Type 2 duration:\n")
cat(sprintf("shape ≈ %.3f, scale ≈ %.3f (mean = shape*scale)\n", gamma_shape_mom, gamma_scale_mom))

cat("\n======================================================================\n\n")

# -----------------------------
# 9) Save a compact CSV for your report
# -----------------------------
out <- data.frame(
  item = c(
    "Type1_duration_mean_min", "Type1_duration_sd_min", "Type1_duration_q90_min",
    "Type1_duration_mean_CI_low", "Type1_duration_mean_CI_high",
    "Type1_duration_sd_CI_low", "Type1_duration_sd_CI_high",
    "Type1_duration_q90_CI_low", "Type1_duration_q90_CI_high",
    "Type2_duration_mean_min", "Type2_duration_sd_min", "Type2_duration_q90_min",
    "Type2_duration_mean_CI_low", "Type2_duration_mean_CI_high",
    "Type2_duration_sd_CI_low", "Type2_duration_sd_CI_high",
    "Type2_duration_q90_CI_low", "Type2_duration_q90_CI_high",
    "Type1_lambda_calls_per_day", "Type1_lambda_CI_low", "Type1_lambda_CI_high",
    "Type1_interarrival_mean_min",
    "Type2_interarrival_mean_min", "Type2_interarrival_sd_min", "Type2_interarrival_q90_min",
    "Type2_interarrival_mean_CI_low", "Type2_interarrival_mean_CI_high",
    "Type2_interarrival_sd_CI_low", "Type2_interarrival_sd_CI_high",
    "Type2_interarrival_q90_CI_low", "Type2_interarrival_q90_CI_high",
    "Type2_dailycount_mean", "Type2_dailycount_var",
    "Optional_Type2_duration_gamma_shape_MOM", "Optional_Type2_duration_gamma_scale_MOM"
  ),
  value = c(
    mu1_hat, sd1_hat, as.numeric(quantile(dur1, 0.90)),
    mu1_ci[1], mu1_ci[2],
    sd1_ci[1], sd1_ci[2],
    q1_ci[1], q1_ci[2],
    mu2_hat, sd2_hat, q2_hat,
    mu2_ci[1], mu2_ci[2],
    sd2_ci[1], sd2_ci[2],
    q2_ci[1], q2_ci[2],
    lambda1_hat, lambda1_ci[1], lambda1_ci[2],
    mean_interarrival1_min,
    i2_mean, i2_sd, i2_q90,
    i2_mean_ci[1], i2_mean_ci[2],
    i2_sd_ci[1], i2_sd_ci[2],
    i2_q90_ci[1], i2_q90_ci[2],
    c2_mean, c2_var,
    gamma_shape_mom, gamma_scale_mom
  )
)

write.csv(out, "OR_bootstrap_summary.csv", row.names = FALSE)
