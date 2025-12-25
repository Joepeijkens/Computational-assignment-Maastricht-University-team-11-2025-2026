rm(list=ls())
setwd("~/University Maastricht SBE/EBS4044 Computational Assignment")

library(readr)

set.seed(67)

case_data <- read_csv("ScanRecords.csv")

data_p1 <- case_data[case_data$PatientType == "Type 1", ]

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

hist(bootstrap_means)
hist(bootstrap_sd)

mean1 <- mean(bootstrap_means)
sd1 <- mean(bootstrap_sd)

mean1*60
sd1*60
#this is the amount of minutes it takes to scan on average and the sd

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

for (i in 1:B){
  bootstrap_sample <- rexp(length(X), rate = 1/mean(X))  
  bootstrap_means[i] <- mean(bootstrap_sample)
  bootstrap_sd[i] <- sd(bootstrap_sample)
}

hist(bootstrap_means)
hist(bootstrap_sd)

mean2 <- mean(bootstrap_means)
sd2 <- mean(bootstrap_sd)

mean2*60
sd2*60











data_p2 <- case_data[case_data$PatientType == "Type 2", ]

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
