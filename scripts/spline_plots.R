## clearing the environment
rm(list = ls())  
gc()    
setwd("/opt/mesh/eigg/sanket/MZ_New_dynamics")


library(rstan)
library(bayesplot)
library(tidyverse)
library(parallel)

## Stan-fit specific details
modelName1 <- "source_counts"
modelName2 <- "Source_chi"
modelName3 <- "Source_ki"
data_derived1 <- "T1_MZ.csv"

## Relative paths assuming the working directory is the script directory
## containing this script
projectDir <- getwd()
scriptDir <- file.path(projectDir, "scripts")
modelDir <- file.path(projectDir, "models")
outputDir <- file.path(projectDir, "output")
saveDir1 <- file.path(outputDir, modelName1)
saveDir2 <- file.path(outputDir, modelName2)
saveDir3 <- file.path(outputDir, modelName3)

## file names for the saved fit objects
output_filename1 = paste(modelName1, "_", substr(data_derived1, 1,2), sep="")
output_filename2 = paste(modelName2, "_", substr(data_derived1, 1,2), sep="")
output_filename3 = paste(modelName3, "_", substr(data_derived1, 1,2), sep="")

## loading saved stanfit object
stanfit1 <- readRDS(file.path(saveDir1, paste(output_filename1, ".rds", sep = "")))
stanfit2 <- readRDS(file.path(saveDir2, paste(output_filename2, ".rds", sep = "")))
stanfit3 <- readRDS(file.path(saveDir3, paste(output_filename3, ".rds", sep = "")))


##################################################################################################
## import the data set
data_imp1 <- read_csv(file.path(projectDir, "data", data_derived1))

# time points in the data
ts_pred1 <- seq(0, 750)
ts_pred2 <- seq(0, 650) 


## Specify the variables for which you want history and density plots
parametersToPlot1 <- c("theta0","nu","sigma_counts")
parametersToPlot2 <- c("chiEst","qEst","sigma_chi")
parametersToPlot3 <- c("ki_0_host","ki_slope_Log", "sigma_ki_host", "ki_0_donor","sigma_ki_donor")

## Additional variables to monitor
otherRVs <- c("y_counts_pred", "countspred", "y_chipred", "chipred",
              "y_kihost_pred","ki_host_pred", "y_kidonor_pred","ki_donor_pred")

parameters <- c(parametersToPlot1, parametersToPlot2, otherRVs)

################################################################################################
## posterior distributions of parameters
### posterior distributions of parameters
ptable1 <- monitor(as.array(stanfit1, pars = parametersToPlot1), warmup = 0, print = FALSE)
ptable1[,c(1,4,8)]
ptable2 <- monitor(as.array(stanfit2, pars = parametersToPlot2), warmup = 0, print = FALSE)
ptable2[,c(1,4,8)]
ptable3 <- monitor(as.array(stanfit3, pars = parametersToPlot3), warmup = 0, print = FALSE)
ptable3[,c(1,4,8)]

################################################################################################
### plotting options
myTheme <- theme(text = element_text(size = 11), axis.text = element_text(size = 11), axis.title =  element_text(size = 11, face = "bold"),
                 plot.title = element_text(size=11,  hjust = 0.5, face = "bold"),
                 legend.background = element_blank(), legend.key = element_blank(),
                 legend.text = element_text(size= 9), legend.title = element_text(9))

# setting ggplot theme for rest fo the plots
theme_set(theme_bw())

####### plotting
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (e^+2 -> e^2)
  l <- gsub("e\\+","e",l)  
  # turn the 'e' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text=l)
}

log10minorbreaks=as.numeric(1:10 %o% 10^(4:8))

################################################################################################
## posterior predictive distributions
Ypred1 <- as.data.frame(stanfit1, pars = "y_counts_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1) %>%  filter(timeseries >=40) %>%  filter(timeseries <=750)

countspred <- as.data.frame(stanfit1, pars = "countspred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1) %>%  filter(timeseries >=40) %>%  filter(timeseries <=750)

Ypred2 <- as.data.frame(stanfit2, pars = "y_chipred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred2) %>%  filter(timeseries <=650)


fdpred <- as.data.frame(stanfit2, pars = "chipred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred2) %>%  filter(timeseries <=650)

Ypred3 <- as.data.frame(stanfit3, pars = "y_kihost_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1) %>% filter(timeseries >=10) %>%  filter(timeseries <= 750)


ki_host_pred <- as.data.frame(stanfit3, pars = "ki_host_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1) %>% filter(timeseries >=10) %>%  filter(timeseries <= 750)

Ypred4 <- as.data.frame(stanfit3, pars = "y_kidonor_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1) %>% filter(timeseries >=10) %>%  filter(timeseries <= 750)

ki_donor_pred <- as.data.frame(stanfit3, pars = "ki_donor_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1) %>% filter(timeseries >=10) %>%  filter(timeseries <= 750)


## plots for counts
p1 <- ggplot() +
  geom_ribbon(data = Ypred1, aes(x = timeseries, ymin = lb, ymax = ub), fill= "#805E80", alpha = 0.4)+
  geom_ribbon(data = countspred, aes(x = timeseries, ymin = lb, ymax = ub), fill= "#805E80", alpha = 0.2)+
  geom_line(data = Ypred1, aes(x = timeseries, y = median), col ="#805E80", size =1) +
  geom_point(data = data_imp1, aes(x = age.at.S1K, y = total_counts)) +
  labs(title=paste('Total numbers of', substr(data_derived1, 1, 2), "B cells"),  y=NULL, x="Host age (days)") + 
  scale_x_continuous(limits = c(40, 750), trans = "log10", breaks = c(75, 150, 300, 600))+
  scale_y_continuous(limits = c(1e5, 3e7), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  guides(color = FALSE)+ myTheme


p2 <- ggplot() +
  geom_hline(aes(yintercept = 1), linetype = 2)+
  geom_ribbon(data = fdpred, aes(x = timeseries, ymin = lb, ymax=ub), fill= "#805E80", alpha = 0.2) +
  geom_ribbon(data = Ypred2, aes(x = timeseries, ymin = lb, ymax = ub), fill= "#805E80", alpha = 0.4)+
  geom_line(data = Ypred2, aes(x = timeseries, y = median), size =1,  col= "#805E80") +
  geom_point(data = data_imp1, aes(x = days.post.bmt, y = fd)) +
  labs(x = "Days post BMT", y = NULL, title = paste("Chimerism within ", substr(data_derived1, 1,2), " B cells",  sep="")) +
  scale_x_continuous(limits = c(0, 675), breaks = c(150, 300, 450, 600))+
  scale_y_continuous(limits =c(-0.1, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1.0))+ 
  guides(color = FALSE)+ myTheme

p3 <- ggplot() +
  geom_ribbon(data = ki_host_pred, aes(x = timeseries, ymin = lb, ymax=ub), fill= "#ff8432", alpha = 0.2) +
  geom_ribbon(data = Ypred3, aes(x = timeseries, ymin = lb, ymax = ub), fill= "#ff8432", alpha = 0.4)+
  geom_line(data = Ypred3, aes(x = timeseries, y = median), size =1, col= "#ff8432") +
  geom_ribbon(data = ki_donor_pred, aes(x = timeseries, ymin = lb, ymax=ub), fill= "#2ca1db", alpha = 0.2) +
  geom_ribbon(data = Ypred4, aes(x = timeseries, ymin = lb, ymax = ub), fill= "#2ca1db", alpha = 0.4)+
  geom_line(data = Ypred4, aes(x = timeseries, y = median), size =1,  col= "#2ca1db") +
  geom_point(data = data_imp1, aes(x = days.post.bmt, y = host_ki67_Tra), col= "#ff8432") +
  geom_point(data = data_imp1, aes(x = days.post.bmt, y = donor_ki67_Tra), col= "#2ca1db") +
  labs(title=paste('Fraction of Ki67High cells within', substr(data_derived1, 1, 2), "B cells"),  y=NULL, x="Host age (days)") + 
  scale_x_continuous(limits = c(10, 750), trans = "log10", breaks = c(20, 60, 200, 600))+ ylim(0,1.2)+
  guides(color = FALSE)+ myTheme


pdf(file = file.path("figures", paste("Spline_plots.pdf", sep = "")),
    width = 8, height = 6, onefile = FALSE, useDingbats = FALSE)

lay1 <- rbind(c(1,1,2,2),
              c(NA,3,3,NA))

gridExtra::grid.arrange(p1, p2, p3, layout_matrix = lay1)


dev.off()

