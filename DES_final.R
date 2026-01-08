install.packages("tidyverse")
install.packages("lubridate")

library(tidyverse)
library(lubridate)

this_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this_dir)
#read data to bootstrap for patient 2
data <- read.csv("./ScanRecords.csv")

data_type2 <- data  %>%
  filter(PatientType =="Type 2")

data_type2_day <- as.numeric(as.Date(data_type2$Date) - min(as.Date(data_type2$Date)) + 1)
data_type2_hour <- data_type2$Time
data_type2_duration <- data_type2$Duration
  
#parameters to test the code
runs <- 500
days <- 21
begin_day <- 8.0
end_day <- 17.0
hours_day <- end_day - begin_day

#for this we need the statistic analysis
arrival_rate_patient1 <- 16.92580386
arrival_rate_patient2 <- 10.01120699
duration_patient1_mean <- 0.428527
duration_patient2_mean <- 0.6784842
duration_patient1_sd <- 0.09723408
duration_patient2_sd <- 0.2001784
#using the RHS of the 95% data interval of normal distribution
length_slot1 <- duration_patient1_mean + 2*duration_patient1_sd
#Using the RHS of the 95% data interval using a bootstrap from the statistical part
length_slot2 <- 0.7061751

#simulation for patients calling
call_simulation <- function(type, lambda = NULL, days = NULL, data_type2_day = NULL, data_type2_hour = NULL){

  if(type == "Type 1") {
    #number of patients of type 1 is poisson distributed
    n_patients <- rpois(1, lambda * days)
  
    #simulation for when patients call
    #creates for each patient a fraction, later rounded down, which represents the day on which they call
    patients_day <- rexp(n_patients, rate=lambda)
    arrival_times <- cumsum(patients_day)
    call_day <- floor(arrival_times) 
    call_hour <- (arrival_times - call_day)*9+8 #converts the fraction of the day to the hour of calling between 8.00 and 17.00
  }
  #since distribution is unkown for type 2 we use bootstrapped data
  else if (type == "Type 2"){
    n_patients <- length(data_type2_day)
    call_day <- sample(data_type2_day, n_patients, replace = TRUE)
    call_hour <- sample(data_type2_hour, n_patients, replace = TRUE)
  }
    
  #create an overview of when patients called
  tibble(
    PatientID = paste0(type, "_", 1:n_patients),
    Patient_type = type,
    Day_called = call_day,
    Hour_called = call_hour
  )
}

#separate scheduling system (old)
schedule_separate <- function(seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  
  #runs simulation twice, for both types
  arrivals_patient1 <- call_simulation("Type 1", arrival_rate_patient1, days)
  arrivals_patient2 <- call_simulation("Type 2", lambda = arrival_rate_patient2, days = days, data_type2_day = data_type2_day, data_type2_hour = data_type2_hour)

#in this part we schedule the patients
schedule <- function(arrivals, length_slot, MRI){
  
  #this vector keeps track of how much of the day is already used by booked appointments
  current_used <- rep(0, days + 100) #we allow for a backlog to be scheduled after the days that we initially run the simulation
  
  #stores the schedule for all patients
  result <- vector("list", nrow(arrivals))
  
  #sets the first possible day that each patient can be scheduled the day after the patient called
  for(i in seq_len(nrow(arrivals))){
    patient <- arrivals[i, ]
    next_day <- patient$Day_called +1
    
    #checks if we can schedule a patient on the day after calling without exceeding 17.00, otherwise it moves to the day after that
    repeat{
      if((current_used[next_day] + length_slot) <= hours_day){
        
        start <- begin_day + current_used[next_day]
        current_used[next_day] <- current_used[next_day] + length_slot
        
        result[[i]]<- tibble(
          PatientID = patient$PatientID,
          Patient_type = patient$Patient_type,
          Day_called = patient$Day_called,
          Day_treatment = next_day,
          Start_time_treatment = start,
          End_time_treatment = start + length_slot,
          Length_Slot = length_slot,
          MRI = MRI
        )
        break
      }
      else{
        next_day <- next_day +1
      }
    }
  }
  bind_rows(result)
}

MRI1_schedule <- schedule(arrivals_patient1, length_slot1, "MRI1")
MRI2_schedule <- schedule(arrivals_patient2, length_slot2, "MRI2")

#chronological table of scheduled patients on both MRI's
complete_schedule <- bind_rows(MRI1_schedule, MRI2_schedule)%>%
  arrange(MRI, Day_treatment, Start_time_treatment)

#create the durations of eacht patient from the normal distribution for type 1 and use bootstrapped data for type 2
data <- complete_schedule %>%
  mutate(
    duration = case_when(
      Patient_type == "Type 1" ~ rnorm(n(), duration_patient1_mean, duration_patient1_sd),
      Patient_type == "Type 2" ~ sample(data_type2_duration, n(), replace = TRUE)
    ),
    #set a lowerbound for the duration since we do not allow the duration to be negative or zero
    duration = pmax(0.1, duration)
  )

#sort by MRI and time to ensure we process the queue in order
data <- data %>% 
  arrange(MRI, Day_treatment, Start_time_treatment)

#initialize columns for actual start and end times
data$Actual_start_time <- data$Start_time_treatment
data$Actual_end_time <- 0 

#loop through every patient to update times based on the previous patient's finish time
for(i in 1:nrow(data)) {
  # check if it is the first patient of the dataset OR a new machine OR a new day (in these cases, the patient starts exactly on their scheduled time)
  if(i == 1 || data$MRI[i] != data$MRI[i-1] || data$Day_treatment[i] != data$Day_treatment[i-1]) {
    data$Actual_start_time[i] <- data$Start_time_treatment[i]
  } else {
    # otherwise, they start at the LATER of: A) their scheduled time or B) when the previous patient actually finished (the ripple effect)
    data$Actual_start_time[i] <- max(data$Start_time_treatment[i], data$Actual_end_time[i-1])
  }
  
  #calculate their actual finish time based on the actual start + random duration
  data$Actual_end_time[i] <- data$Actual_start_time[i] + data$duration[i]
}

#calculate final performance metrics
final_data <- data %>%
  mutate(
    #calculate waiting room delay (difference between actual start and scheduled start)
    Waiting_room_delay = Actual_start_time - Start_time_treatment,
    #check for overtime (did they finish after 17.00?)
    Is_overtime = Actual_end_time > 17.0,
    #calculate how many hours of overtime occurred
    Overtime_amount = pmax(0, Actual_end_time - 17.0),
    #the number of days the patient waited for the appointment
    Access_time_days = Day_treatment - Day_called
  )
list(KPI_patient=final_data)
}

#pooled scheduling system (new)
schedule_pooled <- function(seed=NULL){
  if(!is.null(seed)) set.seed(seed)

  #runs simulation twice, for both types
  arrivals_patient1 <- call_simulation("Type 1", arrival_rate_patient1, days)
  arrivals_patient2 <- call_simulation("Type 2", lambda = arrival_rate_patient2, days = days, data_type2_day = data_type2_day, data_type2_hour = data_type2_hour)
  
  #For the pooled model, we combine EVERYONE into one queue immediately
  #processed first-come-first-served based on when they called
  total_arrivals_pooled <- bind_rows(arrivals_patient1, arrivals_patient2) %>%
    arrange(Day_called, Hour_called)
  
  #implement two trackers, one for each machine
  current_used_mri1 <- rep(0, days + 100)
  current_used_mri2 <- rep(0, days + 100)
  
  result <- vector("list", nrow(total_arrivals_pooled))
  
  for(i in seq_len(nrow(total_arrivals_pooled))){
    patient <- total_arrivals_pooled[i, ]
    
    #determine how long this specific patient takes
    # (Type 1 is always 0.5, Type 2 is always 0.75, regardless of machine)
    needed_slot <- if(patient$Patient_type == "Type 1") length_slot1 else length_slot2
    
    #helper function to find the earliest slot on a specific machine tracker
    find_earliest_slot <- function(tracker, start_search_day, duration) {
      for(check_day in start_search_day:length(tracker)) {
        if((tracker[check_day] + duration) <= hours_day) {
          return(list(day = check_day, start_time = begin_day + tracker[check_day]))
        }
      }
      stop("No feasible slot in tracker")
    }
    
    #best option on MRI 1
    option1 <- find_earliest_slot(current_used_mri1, patient$Day_called + 1, needed_slot)
    #best option on MRI 2
    option2 <- find_earliest_slot(current_used_mri2, patient$Day_called + 1, needed_slot)
    #compare and pick the winner (We prefer the earliest day. If days are equal, we prefer the earliest time.)
    use_MRI1 <- FALSE
    
    if (option1$day < option2$day) {
      use_MRI1 <- TRUE
    } else if (option1$day == option2$day && option1$start_time <= option2$start_time) {
      use_MRI1 <- TRUE
    } else {
      use_MRI1 <- FALSE
    }
    
    #book the appointment
    if(use_MRI1) {
      chosen_MRI <- "MRI1"
      chosen_day <- option1$day
      chosen_start <- option1$start_time
      #update MRI 1 tracker
      current_used_mri1[chosen_day] <- current_used_mri1[chosen_day] + needed_slot
    } else {
      chosen_MRI <- "MRI2"
      chosen_day <- option2$day
      chosen_start <- option2$start_time
      #update MRI 2 tracker
      current_used_mri2[chosen_day] <- current_used_mri2[chosen_day] + needed_slot
    }
    
    result[[i]] <- tibble(
      PatientID = patient$PatientID,
      Patient_type = patient$Patient_type,
      Day_called = patient$Day_called,
      Day_treatment = chosen_day,
      Start_time_treatment = chosen_start,
      End_time_treatment = chosen_start + needed_slot,
      Length_Slot = needed_slot,
      MRI = chosen_MRI
    )
    }
  complete_schedule <- bind_rows(result)%>%
    arrange(MRI, Day_treatment, Start_time_treatment)
  
  data <- complete_schedule %>%
    mutate(
      duration = case_when(
        Patient_type == "Type 1" ~ rnorm(n(), duration_patient1_mean, duration_patient1_sd),
        Patient_type == "Type 2" ~ sample(data_type2_duration, n(), replace = TRUE)
      ),
      #set a lowerbound for the duration since we do not allow the duration to be negative or zero
      duration = pmax(0.1, duration)
    )
  #initialize columns for actual start and end times
  data$Actual_start_time <- data$Start_time_treatment
  data$Actual_end_time <- 0 
  
  #loop through every patient to update times based on the previous patient's finish time
  for(i in 1:nrow(data)) {
    # check if it is the first patient of the dataset OR a new machine OR a new day (in these cases, the patient starts exactly on their scheduled time)
    if(i == 1 || data$MRI[i] != data$MRI[i-1] || data$Day_treatment[i] != data$Day_treatment[i-1]) {
      data$Actual_start_time[i] <- data$Start_time_treatment[i]
    } else {
      # otherwise, they start at the LATER of: A) their scheduled time or B) when the previous patient actually finished (the ripple effect)
      data$Actual_start_time[i] <- max(data$Start_time_treatment[i], data$Actual_end_time[i-1])
    }
    
    #calculate their actual finish time based on the actual start + random duration
    data$Actual_end_time[i] <- data$Actual_start_time[i] + data$duration[i]
  }
  
  #calculate final performance metrics
  final_data <- data %>%
    mutate(
      #calculate waiting room delay (difference between actual start and scheduled start)
      Waiting_room_delay = Actual_start_time - Start_time_treatment,
      #check for overtime (did they finish after 17.00?)
      Is_overtime = Actual_end_time > 17.0,
      #calculate how many hours of overtime occurred
      Overtime_amount = pmax(0, Actual_end_time - 17.0),
      #the number of days the patient waited for the appointment
      Access_time_days = Day_treatment - Day_called
    )
  list(KPI_patient = final_data)
}

#run DES multiple times and compare results
results <- map_dfr(1:runs, function(r){
  
  pooled <- schedule_pooled(seed = 100 + r)$KPI_patient
  separate <- schedule_separate(seed = 100 + r)$KPI_patient
  
  pooled_kpi <- pooled%>%  
  summarise(
      Avg_Waiting_room_delay = mean(Waiting_room_delay),
      Avg_Access_time_days = mean(Access_time_days),
      Fraction_overtime = mean(Is_overtime),
      Avg_Overtime_amount = mean(Overtime_amount)
    ) 
  
  separate_kpi <- separate %>%
    summarise(
      Avg_Waiting_room_delay = mean(Waiting_room_delay),
      Avg_Access_time_days = mean(Access_time_days),
      Fraction_overtime = mean(Is_overtime),
      Avg_Overtime_amount = mean(Overtime_amount)
    ) 
  
  tibble(run = r, KPI = names(pooled_kpi), pooled = as.numeric(pooled_kpi), separate = as.numeric(separate_kpi),difference = as.numeric(pooled_kpi) - as.numeric(separate_kpi))
})

results

summary <- results %>%
  group_by(KPI) %>%
  summarise(mean_pooled= mean(pooled), mean_separate = mean(separate), mean_difference = mean(difference), lb_interval = quantile(difference, 0.025),ub_interval = quantile(difference, 0.975), P_pooled_better = round(mean(difference <0)*100, 1))
summary
