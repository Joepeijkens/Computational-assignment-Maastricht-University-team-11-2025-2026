install.packages("simmer")
library(simmer)

set.seed(42)

env <- simmer("SuperDuperSim")
env
url <- "https://raw.githubusercontent.com/Joepeijkens/Computational-assignment-Maastricht-University-team-11-2025-2026/0cffd7a386b5c0d1e46f01ed35e4fd8cf61d82bb/ScanRecords.csv"

data <- read.csv(url)

patient1 <- trajectory("patients' path") %>%
  seize("call", 1) %>%
  timeout(0) %>%
  release("call",1) %>%
  
  seize("MRI",1)%>%
  timeout(function() rnorm(1, mean = 0.43, sd = 0.097))%>%
  release("MRI",1)

env%>%
  add_resource("call", 1) %>%
  add_resource("MRI", 1) %>%
  add_generator("patient1", patient1, function()rpois(1, lambda = 17))


env %>% run(365)
