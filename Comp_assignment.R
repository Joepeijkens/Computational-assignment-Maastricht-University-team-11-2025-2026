rm(list=ls())
setwd("~/University Maastricht SBE/EBS4044 Computational Assignment")

library(readr)

set.seed(67)

case_data <- read_csv("ScanRecords.csv")

##### PATIENT TYPE 1 #####

data_p1 <- case_data[case_data$PatientType == "Type 1", ]

summary(data_p1[c(3)])
date_counts <- table(data_p1$Date)
summary(as.numeric(date_counts))
length(as.numeric(date_counts))

#1 will refer to scan times
#2 will refer to number per patients

X <- data_p1$Duration

dur_mean1 <- mean(X)
dur_sd1 <- sd(X)

n <- length(X)
B <- 5000
bootstrap_means <- 1:B
bootstrap_sd <- 1:B

for (i in 1:B){
  bootstrap_sample <- rnorm(n, mean = dur_mean1, sd = dur_sd1)  
  bootstrap_means[i] <- mean(bootstrap_sample)
  bootstrap_sd[i] <- sd(bootstrap_sample)
}

mean1 <- mean(bootstrap_means)
sd1 <- mean(bootstrap_sd)

mean1*60
sd1*60
#this is the amount of minutes it takes to scan on average and the sd
hist(bootstrap_means)
abline(v = mean(bootstrap_means), col = "black", lty = 1, lwd = 2)
abline(v = mean(bootstrap_means)+(1.96*sd(bootstrap_means)),
       col = "red", lty = 2, lwd = 2)
abline(v = mean(bootstrap_means)-(1.96*sd(bootstrap_means)),
       col = "red", lty = 2, lwd = 2)
legend("topleft", legend = c("Mean bootstraps", "95% confidence interval"), 
       col = c("black", "red"), lty = c(2,2), cex = 0.65)
hist(bootstrap_sd)
abline(v = mean(bootstrap_sd), col = "black", lty = 1, lwd = 2)
abline(v = mean(bootstrap_sd)+(1.96*sd(bootstrap_sd)),
       col = "red", lty = 2, lwd = 2)
abline(v = mean(bootstrap_sd)-(1.96*sd(bootstrap_sd)),
       col = "red", lty = 2, lwd = 2)
legend("topleft", legend = c("Mean bootstraps", "95% confidence interval"), 
       col = c("black", "red"), lty = c(2,2), cex = 0.61)
#Note using quantiles for CI does not make a meaningful difference.

data_p1$Date <- as.Date(data_p1$Date)

unique_dates <- unique(data_p1$Date)

first_last_df <- data.frame(
  Date = character(),
  Time = numeric(),
  Duration = numeric(),
  PatientType = character()
)

for (i in seq_along(unique_dates)) {
  current_date <- unique_dates[i]
  date_rows <- data_p1[data_p1$Date == current_date, ]
  first_row <- date_rows[1, ]
  last_row <- date_rows[nrow(date_rows), ]
  first_last_df <- rbind(first_last_df, first_row)
  first_last_df <- rbind(first_last_df, last_row)
}

first_last_df <- first_last_df[2:41, ] #Don't count last and first

calc_results <- numeric(40)

# Loop through i = 1, 3, 5, ..., 39
for (i in seq(1, 39, by = 2)) {
  calc_results[i] <- 17 - first_last_df$Time[i]
}

# Loop through j = 2, 4, 6, ..., 40
for (j in seq(2, 40, by = 2)) {
  calc_results[j] <- first_last_df$Time[j] - 8
}

pair_sums <- numeric()

# Loop through pairs (1-2, 3-4, 5-6, ...)
for (i in seq(1, length(calc_results), by = 2)) {
  # Check if there's a pair
  if (i + 1 <= length(calc_results)) {
    pair_sum <- calc_results[i] + calc_results[i + 1]
    pair_sums <- c(pair_sums, pair_sum)
  } else {
    # If odd number of elements, just use the last single value
    pair_sums <- c(pair_sums, calc_results[i])
  }
}

X <- diff(data_p1$Time)
sum(X < 0) #Check if amount is 20 since that is amount of transitions between workdays in aug.

neg_indices <- which(X < 0)

# Replace negative values with pair_sums in order
if (length(pair_sums) >= length(neg_indices)) {
  X[neg_indices] <- pair_sums[1:length(neg_indices)]
} else {
  # If there aren't enough pair_sums, fill as many as possible
  X[neg_indices[1:length(pair_sums)]] <- pair_sums
  # The remaining negative values stay as they are (or you could set to NA)
}

summary(X)
length(X)

for (i in 1:B){
  bootstrap_sample <- rexp(length(X), rate = 1/mean(X))  
  bootstrap_means[i] <- mean(bootstrap_sample)
  bootstrap_sd[i] <- sd(bootstrap_sample)
}

mean2 <- mean(bootstrap_means)

mean2*60

hist(bootstrap_means)
ci <- quantile(bootstrap_means, c(0.025, 0.975))
abline(v = mean(bootstrap_means), col = "black", lty = 1, lwd = 2)
abline(v = ci[1], col = "red", lty = 2, lwd = 2)
abline(v = ci[2],
       col = "red", lty = 2, lwd = 2)
legend("topleft", legend = c("Mean bootstraps", "95% confidence interval"), 
       col = c("black", "red"), lty = c(2,2), cex = 0.65)

##### PATIENT TYPE 2 #####

data_p2 <- case_data[case_data$PatientType == "Type 2", ]

summary(data_p2[c(3)])
date_counts <- table(data_p2$Date)
summary(as.numeric(date_counts))
length(as.numeric(date_counts))

B <- 5000
n <- nrow(data_p2)
X <- data_p2$Duration

bootstrap_means <- 1:5000
bootstrap_sd <- 1:5000

for (i in 1:B){
  bootstrap_sample <- sample(X, n, replace = TRUE)
  bootstrap_means[i] <- mean(bootstrap_sample)
  bootstrap_sd[i] <- sd(bootstrap_sample)
}

hist(bootstrap_means)
hist(bootstrap_sd)

mean3 <- mean(bootstrap_means)
sd3 <- mean(bootstrap_sd)
mean3*60
sd3*60

plot(ecdf(bootstrap_means))
plot(ecdf(bootstrap_sd))

unique_dates <- unique(data_p2$Date)

first_last_df <- data.frame(
  Date = character(),
  Time = numeric(),
  Duration = numeric(),
  PatientType = character()
)

for (i in seq_along(unique_dates)) {
  current_date <- unique_dates[i]
  date_rows <- data_p2[data_p2$Date == current_date, ]
  first_row <- date_rows[1, ]
  last_row <- date_rows[nrow(date_rows), ]
  first_last_df <- rbind(first_last_df, first_row)
  first_last_df <- rbind(first_last_df, last_row)
}

first_last_df <- first_last_df[2:41, ] #Don't count last and first

calc_results <- numeric(40)

# Loop through i = 1, 3, 5, ..., 39
for (i in seq(1, 39, by = 2)) {
  calc_results[i] <- 17 - first_last_df$Time[i]
}

# Loop through j = 2, 4, 6, ..., 40
for (j in seq(2, 40, by = 2)) {
  calc_results[j] <- first_last_df$Time[j] - 8
}

pair_sums <- numeric()

# Loop through pairs (1-2, 3-4, 5-6, ...)
for (i in seq(1, length(calc_results), by = 2)) {
  # Check if there's a pair
  if (i + 1 <= length(calc_results)) {
    pair_sum <- calc_results[i] + calc_results[i + 1]
    pair_sums <- c(pair_sums, pair_sum)
  } else {
    # If odd number of elements, just use the last single value
    pair_sums <- c(pair_sums, calc_results[i])
  }
}

X <- diff(data_p2$Time)
sum(X < 0) #Check if amount is 20 since that is amount of transitions between workdays in aug.

neg_indices <- which(X < 0)

# Replace negative values with pair_sums in order
if (length(pair_sums) >= length(neg_indices)) {
  X[neg_indices] <- pair_sums[1:length(neg_indices)]
} else {
  # If there aren't enough pair_sums, fill as many as possible
  X[neg_indices[1:length(pair_sums)]] <- pair_sums
  # The remaining negative values stay as they are (or you could set to NA)
}

for (i in 1:B){
  bootstrap_sample <- sample(X, n, replace = TRUE) 
  bootstrap_means[i] <- mean(bootstrap_sample)
  bootstrap_sd[i] <- sd(bootstrap_sample)
}

hist(bootstrap_means)
hist(bootstrap_sd)

mean4 <- mean(bootstrap_means)
sd4 <- mean(bootstrap_sd)

mean4*60
sd4*60

plot(ecdf(bootstrap_means))
plot(ecdf(bootstrap_sd))
