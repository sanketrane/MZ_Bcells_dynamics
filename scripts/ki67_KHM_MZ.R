## Kinetic heterogeneity model 

## clearing the environment and freeing the memory held by previously loaded objects
rm(list = ls()); gc()   

# setting working dirctory
#setwd("/opt/mesh/eigg/sanket/MZ_New_dynamics")


#####!/usr/bin/env Rscript
####args = commandArgs(trailingOnly = TRUE)
####
##### test if there is at least one argument: if not, return an error
####if (length(args)==0) {
####  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
####} else if (length(args)==1) {
####  # default output file
####  args[2] = "out.txt"
####}

####################################################################################
## Installing r-stan pachage on the go:
if(!("rstan" %in% rownames(installed.packages())) ){
  install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies=TRUE)
}

## Installing loo:
if(!("loo" %in% rownames(installed.packages())) ){
  install.packages("loo", repos = "https://cloud.r-project.org/", dependencies=TRUE)
}

## Installing tidyverse:
if(!("tidyverse" %in% rownames(installed.packages())) ){
  install.packages("tidyverse", repos = "https://cloud.r-project.org/", dependencies=TRUE)
}
####################################################################################

# loading libararies
library(rstan)
library(parallel)
library(loo)
library(tidyverse)

# model sepcific details
## data files
data_derived1 <- "T1_MZ.csv"    # name of the file for precursor pop
data_derived2 <- "counts_MZ.csv"  # name of the file for counts data of traget pop
data_derived3 <- "Nfd_MZ.csv"     # name of the file for normalised fd data of target pop
data_derived4 <- "Ki67_MZ.csv"    # name of the file for ki67 data of  pop

## model files (Stan)
modelName <- "ki67_KHM_MZ"        # name of the file for stan model 
#### in some cases name of the stan file is depedent on the name of the precursor pop used
## Setting all the directories for opeartions
projectDir <- getwd()
scriptDir <- file.path(projectDir, "scripts")
modelDir <- file.path(projectDir, "models")
dataDir <- file.path(projectDir, "data")
toolsDir <- file.path(scriptDir, "tools")
outputDir <- file.path(projectDir, "output")
saveDir <- file.path(outputDir, paste(modelName, "_", substr(data_derived1, 1,2), sep=""))  ### varies according to the precursor pop used

# parallelisation for R-Stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# set different seeds for different tasks
### if the script is called on different nodes across the cluster then it will use different seed so as to optimise the sampling 
i <- as.numeric(Sys.getenv("SLURM_PROCID"))
#seed <- args[1]                     ### Use your lucky mumber! haha! 
seed <- 1890 + i                     ### Use your lucky mumber! haha! 

################################################################################################
################################################################################################
# loading data files
#source_sorted <- read_csv(file.path(dataDir, data_derived1))%>% arrange(age.at.S1K)
counts_sorted <- read_csv(file.path(dataDir, data_derived2))%>% arrange(age.at.S1K)
Nfd_sorted <- read_csv(file.path(dataDir, data_derived3))%>% arrange(age.at.S1K)
ki67_sorted <- read_csv(file.path(dataDir, data_derived4))%>% arrange(age.at.S1K)

### solving ODEs only for unique timepoints 
### also Stan ODE solver throws an error for repetitive time-points 
# unique time points in data for odes solver
unique_times_df <- counts_sorted %>% distinct(age.at.S1K, .keep_all = TRUE)

data_time <- counts_sorted$age.at.S1K                                          # data and solver time 
solve_time <- unique_times_df$age.at.S1K                                       # unique time points to solve ode
time_index <- purrr::map_dbl(data_time, function(x) which(x == solve_time))    # keeping track of index of time point in relation to solve_time


# delay time for age correction 
# for each host with different age at BMT -- time zero is different which is becomes another variable in the ODEs -- is it delayed differential equation now?
# tb_time -- for solving for initial conditions at each tb --  we use the earlist host age at BMT to calculate the initial conditions at each tb
solve_ageAtBMT <- unique_times_df$age.at.bmt   
tb_time <- solve_ageAtBMT %>% unique() %>% sort()                              # unique age at BMT points to solve ode
tb_index <- purrr::map_dbl(solve_ageAtBMT, function(x) which(x == tb_time))    #keeping track of index of  delay time in relation to tb_time


# time sequence for predictions specific to age bins within the data
# In this analysis 3 age bins were selected with mean ages for each bin as 48, 72 and 120 respectively.
# time sequence for predictions specific to age bins within the data
# In this analysis 3 age bins were selected with mean ages for each bin as 48, 72 and 120 respectively.
counts_agebmt_sorted <- counts_sorted %>% arrange(age.at.bmt)
Median_agebin1 <- round(median(counts_agebmt_sorted$age.at.bmt[1:10]))
Median_agebin2 <- round(median(counts_agebmt_sorted$age.at.bmt[11:30]))
Median_agebin3 <- round(median(counts_agebmt_sorted$age.at.bmt[31:42]))

ts_pred1 <- seq(Median_agebin1, to = 750, by = 1)
tb_pred1 <- rep(Median_agebin1, length(ts_pred1))
ts_pred2 <- seq(Median_agebin2, to = 750, by = 1)
tb_pred2 <- rep(Median_agebin2, length(ts_pred2))
ts_pred3 <- seq(Median_agebin3, to = 750, by = 1)
tb_pred3 <- rep(Median_agebin3, length(ts_pred3))

# tb_time_pred -- for solving for initial conditions at each tb (here the mean age of each age bin)
tb_time_pred1 <- c(min(counts_sorted$age.at.bmt), Median_agebin1)
tb_time_pred2 <- c(min(counts_sorted$age.at.bmt), Median_agebin2)
tb_time_pred3 <- c(min(counts_sorted$age.at.bmt), Median_agebin3)


# setting data input according to the precurosr pop used
## these values dictate how splines for counts and chimerism change with time.
## obtained by fitting splines separately to counts and chimerism in the respective precurosr pop
source_list = data.frame("nu" = 1.86e-03, "theta0" = 14.36, "chiEst" = 0.76, "qEst" = 0.094, "eps_donor" = 0.99, "eps_host" = 0.95)


## data list that feeds into the Stan model
## Use names that exactly match with the Stan file
data <- list(
  numObs = nrow(counts_sorted),              # Total number of observations within data
  Nd_0 = 0,                                  # donor counts at t0
  solve_time = solve_time,                   # unique observations for ode solver
  num_index = length(solve_time),            # to declare dimensions of the array in Stan
  num_tb <- length(tb_time),                 # to declare dimensions of the array in Stan
  time_index = time_index,
  tb_index = tb_index,
  ageAtBMT = solve_ageAtBMT,
  tb_time = tb_time,
  dpBMT = counts_sorted$age.at.S1K - counts_sorted$age.at.bmt,       #to be used for the spline of chimerism of the precurosr pop
  counts = counts_sorted$total_counts,       ## data to fit
  Nfd = Nfd_sorted$Nfd,                      ## data to fit
  ki_host = ki67_sorted$host_ki67_MZ,        ## data to fit
  ki_donor = ki67_sorted$donor_ki67_MZ,      ## data to fit
  numPred1 = length(ts_pred1),                # to declare dimensions of the array in Stan
  numPred2 = length(ts_pred2),                # to declare dimensions of the array in Stan
  numPred3 = length(ts_pred3),                # to declare dimensions of the array in Stan
  ts_pred1 = ts_pred1,
  ts_pred2 = ts_pred2,
  ts_pred3 = ts_pred3,
  tb_pred1 = tb_pred1,
  tb_pred2 = tb_pred2,
  tb_pred3 = tb_pred3,
  tb_time_pred1 = tb_time_pred1,
  tb_time_pred2 = tb_time_pred2,
  tb_time_pred3 =tb_time_pred3,
  theta0 = exp(source_list$theta0),        #size of precurosr pop ar t0
  nu = source_list$nu,                     #rate of change of precurosr pop with time
  chiEst = source_list$chiEst,             #stabilised value of chimerism for precurosr pop
  qEst = source_list$qEst,                 #rate of change of chimerism in precurosr pop
  eps_donor = source_list$eps_donor,       #proportion of ki67+ cells in donor precurosr pop
  eps_host = source_list$eps_host        #proportion of ki67+ cells in host precurosr pop
)

## create initial guesstimates
init <- function() list(
  psi = exp(rnorm(1, log(0.5), 0.2)),
  Beta = exp(rnorm(1,log(3.5), 0.1)),
  rhoFast_Log = rnorm(1, -4, 0.1),
  rhoSlow_Log = rnorm(1, -4, 0.1),
  lambdaSlow_Log = rnorm(1, -4, 0.1),
  lambdaFast_Log = rnorm(1, -3, 0.1),
  alpha = exp(rnorm(1, log(0.5), 0.2)),
  
  y0_Log = rnorm(1, 14 , 0.1),
  kappaF_0 = exp(rnorm(1, log(0.4), 0.3)),
  kappaS_0 = exp(rnorm(1, log(0.2), 0.3)),
  
  sigma1 = exp(rnorm(1,log(0.5), 1)),
  sigma2 = exp(rnorm(1,log(3), 1)),
  sigma3 = exp(rnorm(1,log(0.5), 1)),
  sigma4 = exp(rnorm(1,log(1), 1)))

## Specify the variables for which you want history and density plots
parametersToPlot <- c("psi", "f_fast", "rhoFast", "rhoSlow", "lambdaSlow", "lambdaFast", "alpha",
                      "Beta", "y0_Log", "kappaF_0", "kappaS_0", "sigma1", "sigma2", "sigma3", "sigma4")

## Additional variables to monitor
otherRVs <- c("y1_mean_pred_age1", "countspred_age1", "y1_mean_pred_age2", "countspred_age2", "y1_mean_pred_age3", "countspred_age3",
              "y2_mean_pred_age1", "fdpred_age1", "y2_mean_pred_age2", "fdpred_age2", "y2_mean_pred_age3", "fdpred_age3", 
              "y3_mean_pred1",  "donor_kiprop_pred1", "y4_mean_pred1",  "host_kiprop_pred1",
              "y3_mean_pred2",  "donor_kiprop_pred2", "y4_mean_pred2",  "host_kiprop_pred2",
              "y3_mean_pred3",  "donor_kiprop_pred3", "y4_mean_pred3",  "host_kiprop_pred3",
              "log_lik1", "log_lik2", "log_lik3", "log_lik4", 
              "lambdaFast_inv", "lambdaSlow_inv", "rhoFast_inv", "deltaSlow", "deltaSlow_inv", "deltaFast", "deltaFast_inv")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
## run Stan

# parameters for running fits
nChains <- 3
nPost <- 1000 ## Number of post-burn-in samples per chain after thinning
nBurn <- 500 ## Number of burn-in samples per chain after thinning
nThin <- 1


nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

fit <- stan(file = file.path(modelDir, paste(modelName, ".stan", sep = "")),
            data = data,
            pars = parameters,
            iter = nIter,
            warmup = nBurnin,
            thin = nThin, 
            init = init,
            control = list(adapt_delta = 0.9),
            chains = nChains)

################################################################################################
# save results in output directory
if (!file.exists(saveDir)){
  dir.create(saveDir)
}

### distinguishing models runs from each other and gathering information about individual tasks

# saving output file as 
output_filename=paste(Sys.getenv("SLURM_JOB_ID"), "_P", Sys.getenv("SLURM_PROCID"), "_",  modelName, ".rds",  sep="")

# saving the stan fit object for individual runs as 
saveRDS(fit, file = file.path(saveDir, output_filename))

## calculating an output from individual runs for the validation of a sucessful run
loo_loglik <- extract_log_lik(fit, parameter_name = "log_lik", merge_chains = TRUE)
loo_ic <- loo(loo_loglik)

## this value is added to the job specific run_info file which is gets written in the save directory after each job copmpletion

# information about the cluster run
cluster_run_info <- paste("Time = ", timestamp(),
                          "|| Job ID = ", Sys.getenv("SLURM_JOB_ID"),
                          "|| Modelname = ", modelName,
                          "|| crude_loo_ic = ", loo_ic$estimates[3],
                          "|| Host = ", Sys.getenv("SLURM_SUBMIT_HOST"),
                          "|| Nodes = ",  Sys.info()["nodename"],
                          "|| Process_ID = P", Sys.getenv("SLURM_PROCID"))

# writing the cluster run info as a data file in the out dir 
write.table(cluster_run_info, file = file.path(saveDir, paste("job_", Sys.getenv("SLURM_JOB_ID"), "_run_info.txt")), append = TRUE)

