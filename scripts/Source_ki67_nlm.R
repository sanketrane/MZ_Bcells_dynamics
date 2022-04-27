## Simple linear ODE model -- Homogeneous model 
## clearing the environment
rm(list = ls())  
gc()    
setwd("/opt/mesh/eigg/sanket/MZ_New_dynamics")


library(rstan)
library(bayesplot)
library(tidyverse)
library(parallel)

modelName <- "Source_ki"
data_derived1 <- "T1_MZ.csv"

## Relative paths assuming the working directory is the script directory
## containing this script
projectDir <- getwd()
scriptDir <- file.path(projectDir, "scripts")
figDir <- file.path(projectDir, "deliv", "figures", modelName)
tabDir <- file.path(projectDir, "deliv", "tables", modelName)
dataDir <- file.path(projectDir, "data", data_derived1)
modelDir <- file.path(projectDir, "models")
outputDir <- file.path(projectDir, "output")
saveDir <- file.path(outputDir, modelName)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#set.seed(1027) ## not required but assures repeatable results

################################################################################################

## import the data set
data_imp <- read_csv(dataDir) %>% arrange(days.post.bmt)

ts_pred <- seq(0, 750)

## create data set
data <- with(data_imp,
             list(
               numObs = nrow(data_imp),
               Time = data_imp$days.post.bmt,
               numPred = length(ts_pred),
               ts_pred = ts_pred,
               ki_host = data_imp$host_ki67_Tra,
               ki_donor = data_imp$donor_ki67_Tra))

  
## create initial estimates
init <- function() list(
  ki_0_host = rnorm(1, 0.8, 0.01),
  ki_slope_HLog = rnorm(1, -5, 0.1),
  ki_0_donor = rnorm(1, 0.9, 0.001),
  ki_slope_DLog = rnorm(1, -5, 0.1),
  
  sigma_ki_host = exp(rnorm(1,log(1.5), 1)),
  sigma_ki_donor = exp(rnorm(1,log(1.5), 1)))

## Specify the variables for which you want history and density plots
parametersToPlot <- c("ki_0_host","ki_slope_Log", "sigma_ki_host", "ki_0_donor","sigma_ki_donor")

## Additional variables to monitor
otherRVs <- c("y_kihost_pred","ki_host_pred", "y_kidonor_pred","ki_donor_pred")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
# run Stan

nChains <- 10
nPost <- 1500 ## Number of post-burn-in samples per chain after thinning
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

# save results in output directory
if (!file.exists(saveDir)){
  dir.create(saveDir)
}

# saving output file as 
output_filename = paste(modelName, "_", substr(data_derived1, 1,2), sep="")


saveRDS(fit, file = file.path(saveDir, paste(output_filename, ".rds", sep = "")))
