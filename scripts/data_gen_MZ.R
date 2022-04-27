### data wrangling for Marginal Zone compartmnet
rm(list = ls()); gc();

library(tidyverse)

## importing data to be fitted 
counts_file <- file.path("datafiles", "counts_MZ.csv")
counts_data <- read.csv(counts_file) %>% 
  arrange(age.at.S1K)

Nfd_file <- file.path("datafiles", "Nfd_MZ.csv")
Nfd_data <- read.csv(Nfd_file) %>% 
  arrange(age.at.S1K) %>%
  filter(Nfd <= 1.2)

ki_file <- file.path("datafiles", "Ki67_MZ.csv")
ki_data <- read.csv(ki_file) %>% 
  arrange(age.at.S1K) %>%
  rename(host_ki = host_ki67_MZ,
         donor_ki = donor_ki67_MZ)

## pooled data
chimera_data <- full_join(counts_data, Nfd_data, 
                          by = c("Lamis.ID", "age.at.S1K", "age.at.bmt", "days.post.bmt")) %>%
  full_join(ki_data, 
            by = c("Lamis.ID", "age.at.S1K", "age.at.bmt", "days.post.bmt")) %>%
  ### Binning the data for easy predictions
  mutate(age_bins = ifelse(age.at.bmt <= 56, "age_bin1",
                           ifelse(age.at.bmt <= 77, "age_bin2", "age_bin3")))

## defining the function to calculate mode of a vector series
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

## mean of age at BMT within each age at BMT bin
chimera_data %>% arrange(age.at.bmt) %>%
  group_by(age_bins) %>%
  summarise("Mode"= getmode(age.at.bmt),
            "Mean"= mean(age.at.bmt))


## Unique time points with indices to map
unique_times_chi <- chimera_data %>% distinct(age.at.S1K, .keep_all = TRUE) 
data_time_chi <- chimera_data$age.at.S1K 
solve_time_chi <- unique_times_chi$age.at.S1K  ## unique time points in the data
## Map of the unique time points on all the timepoints
time_index_chi <- purrr::map_dbl(data_time_chi, function(x) which(x == solve_time_chi))    # keeping track of index of time point in relation to solve_time

## Data to import in Stan
numObs <- length(data_time_chi)
n_shards <- length(solve_time_chi)
solve_time <- solve_time_chi
time_index <- time_index_chi
dpBMT <- chimera_data$age.at.S1K -chimera_data$age.at.bmt
ageAtBMT <- unique_times_chi$age.at.bmt
counts <- chimera_data$total_counts
Nfd <- chimera_data$Nfd
ki_donor <- chimera_data$donor_ki
ki_host <- chimera_data$host_ki
# time sequence for predictions specific to age bins within the data
ts_pred1 <- 10^seq(log10(45), log10(750), length.out = 300)
ts_pred2 <- 10^seq(log10(67), log10(750), length.out = 300)
ts_pred3 <- 10^seq(log10(89), log10(750), length.out = 300)
tb_pred1 <- rep(45, 300)
tb_pred2 <- rep(67, 300)
tb_pred3 <- rep(89, 300)
numPred <- length(ts_pred1)


stan_rdump(c("numObs",  "n_shards", "solve_time", "dpBMT", "ageAtBMT", "time_index",
             "counts",  "Nfd", "ki_donor", "ki_host",
             "ts_pred1", "ts_pred2", "ts_pred3",
             "tb_pred1", "tb_pred2", "tb_pred3", "numPred"),
             file = file.path('datafiles', paste0('MZ_data',".Rdump")))






















