# OR_inputs.R
# Computes the OR “input model” quantities requested from ScanRecords.csv:
# (1) Type 1 duration: Normal mean & sd
# (2) Type 2 duration: mean, sd, and a recommended distribution type (Gamma) + fitted params
# (3) Type 1 arrivals: Poisson lambda (calls/day) and implied interarrival Exp mean
# (4) Type 2 interarrival times: empirical (nonparametric) distribution + key stats
#
# Run: Rscript OR_inputs.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(MASS)        # fitdistr
})

CSV_PATH <- "ScanRecords.csv"

WORK_START_H <- 8
WORK_END_H   <- 17
DAY_LEN_MIN  <- (WORK_END_H - WORK_START_H) * 60  # 540 minutes

# -------------------------
# 1) Load + preprocess
# -------------------------
df <- read_csv(CSV_PATH, show_col_types = FALSE) %>%
  mutate(
    Date = as.Date(Date),
    Day  = Date,
    CallMin     = (Time - WORK_START_H) * 60,  # minutes since 08:00
    DurationMin = Duration * 60
  )

df1 <- df %>% filter(PatientType == "Type 1")
df2 <- df %>% filter(PatientType == "Type 2")

# daily counts
daily_counts <- df %>%
  count(Day, PatientType) %>%
  pivot_wider(names_from = PatientType, values_from = n, values_fill = 0) %>%
  arrange(Day)

c1 <- daily_counts$`Type 1`
c2 <- daily_counts$`Type 2`

# -------------------------
# 2) Type 1 duration: Normal mean & sd
# -------------------------
mu1_min <- mean(df1$DurationMin)
sd1_min <- sd(df1$DurationMin)

# -------------------------
# 3) Type 2 duration: mean, sd, and distribution type + params
#     Recommended: Gamma (positive support, right-skewed)
# -------------------------
mu2_min <- mean(df2$DurationMin)
sd2_min <- sd(df2$DurationMin)

# Fit Gamma by MLE: shape/rate
# (Gamma mean = shape/rate, var = shape/rate^2)
gamma_fit <- suppressWarnings(fitdistr(df2$DurationMin, densfun = "gamma"))
g_shape <- unname(gamma_fit$estimate["shape"])
g_rate  <- unname(gamma_fit$estimate["rate"])
g_scale <- 1 / g_rate

# -------------------------
# 4) Type 1 arrivals: Poisson lambda + implied exponential interarrival
# -------------------------
lambda1_day <- mean(c1)                 # calls per working day
rate1_per_min <- lambda1_day / DAY_LEN_MIN  # calls per minute during working hours
mean_interarrival1_min <- 1 / rate1_per_min # minutes

# (Optional) quick check: Poisson-ish daily counts (mean~var)
mean_c1 <- mean(c1)
var_c1  <- var(c1)

# -------------------------
# 5) Type 2 interarrival times: empirical (nonparametric) + key stats
# -------------------------
# Compute within-day interarrival times (minutes) from observed call times
type2_interarrivals <- df2 %>%
  arrange(Day, CallMin) %>%
  group_by(Day) %>%
  summarise(
    inter = list(diff(sort(CallMin))),
    .groups = "drop"
  ) %>%
  pull(inter) %>%
  unlist()

# If some days have 0/1 call -> diff() gives length 0; unlist drops them automatically
type2_interarrivals <- as.numeric(type2_interarrivals)
type2_interarrivals <- type2_interarrivals[is.finite(type2_interarrivals) & type2_interarrivals >= 0]

# Summary stats for the empirical distribution
summ_stats <- function(x) {
  tibble(
    n = length(x),
    mean = mean(x),
    sd = sd(x),
    p50 = quantile(x, 0.50, names = FALSE),
    p90 = quantile(x, 0.90, names = FALSE),
    p95 = quantile(x, 0.95, names = FALSE)
  )
}

type2_inter_stats <- summ_stats(type2_interarrivals)

# Also: Type 2 daily arrivals mean/var (to justify non-Poisson if overdispersed)
mean_c2 <- mean(c2)
var_c2  <- var(c2)

# -------------------------
# 6) Print OR-ready answers
# -------------------------
cat("\n================ OR INPUT MODEL OUTPUTS ================\n")

cat("\n(1) Patient Type 1 scan duration (Normal assumption)\n")
cat(sprintf("    Mean (min): %.3f\n", mu1_min))
cat(sprintf("    SD   (min): %.3f\n", sd1_min))

cat("\n(2) Patient Type 2 scan duration\n")
cat(sprintf("    Mean (min): %.3f\n", mu2_min))
cat(sprintf("    SD   (min): %.3f\n", sd2_min))
cat("    Recommended distribution type: Gamma (positive, right-skewed)\n")
cat(sprintf("    Fitted Gamma params (MLE): shape = %.4f, scale = %.4f (rate = %.4f)\n",
            g_shape, g_scale, g_rate))

cat("\n(3) Patient Type 1 arrivals\n")
cat("    Daily arrivals modeled as: N1 ~ Poisson(lambda1)\n")
cat(sprintf("    Estimated lambda1 (calls/day): %.3f\n", lambda1_day))
cat(sprintf("    Implied constant arrival rate during working hours: %.6f calls/min\n", rate1_per_min))
cat("    Interarrival time modeled as: T1 ~ Exponential(rate = lambda1/540)\n")
cat(sprintf("    Mean interarrival time (min): %.3f\n", mean_interarrival1_min))
cat(sprintf("    Daily count check (Type 1): mean = %.3f, var = %.3f\n", mean_c1, var_c1))

cat("\n(4) Patient Type 2 time between consecutive arrivals\n")
cat("    Distribution type: Empirical / Nonparametric (resampling from observed interarrivals)\n")
print(type2_inter_stats)
cat(sprintf("    Daily count check (Type 2): mean = %.3f, var = %.3f (var>mean suggests non-Poisson)\n",
            mean_c2, var_c2))

# -------------------------
# 7) Save a clean table for your report
# -------------------------
out_tbl <- tibble(
  item = c("Type1_duration_mean_min", "Type1_duration_sd_min",
           "Type2_duration_mean_min", "Type2_duration_sd_min",
           "Type2_duration_dist", "Type2_gamma_shape", "Type2_gamma_scale",
           "Type1_lambda_calls_per_day", "Type1_arrival_rate_calls_per_min", "Type1_mean_interarrival_min",
           "Type2_interarrival_dist", "Type2_interarrival_mean_min", "Type2_interarrival_sd_min",
           "Type2_interarrival_p50_min", "Type2_interarrival_p90_min", "Type2_interarrival_p95_min"),
  value = c(mu1_min, sd1_min,
            mu2_min, sd2_min,
            "Empirical durations in inference; Gamma fitted for simulation input", g_shape, g_scale,
            lambda1_day, rate1_per_min, mean_interarrival1_min,
            "Empirical (nonparametric) interarrival times from data", type2_inter_stats$mean, type2_inter_stats$sd,
            type2_inter_stats$p50, type2_inter_stats$p90, type2_inter_stats$p95)
)

write_csv(out_tbl, "OR_inputs_summary.csv")

cat("\nSaved: OR_inputs_summary.csv\n")
cat("========================================================\n\n")
