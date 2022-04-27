### data wrangling for Marginal Zone compartmnet
rm(list = ls()); gc();

library(tidyverse)

## importing data to be fitted 
T2MZP_file <- file.path("datafiles", "T2MZp_data.csv")
T2MZP_data <- read.csv(T2MZP_file) %>% arrange(age.at.S1K) %>%
  mutate(total_counts = T2.MZP_host + T2.MZP_donor,
         fd = (T2.MZP_donor)/total_counts,
         total_ki = (T2.MZP_host_Ki67 + T2.MZP_donor_Ki67)/total_counts,
         host_counts = T2.MZP_host,
         donor_counts = T2.MZP_donor,
         # ki67 fractions without the immature FM compartment
         host_ki67 = T2.MZP_host_Ki67/T2.MZP_host,
         donor_ki67 = T2.MZP_donor_Ki67/T2.MZP_donor) %>%
  select(Lamis.ID, contains('age'), total_counts, fd, total_ki, host_counts, donor_counts, host_ki67, donor_ki67) %>%
  na.omit()


ggplot(T2MZP_data) +
  geom_point(aes(x=age.at.S1K, y=total_counts)) + scale_y_log10(limits=c(1e3, 5e4))+ scale_x_log10(limits=c(50, 1000))

ggplot(T2MZP_data) +
  geom_point(aes(x=age.at.S1K - age.at.BMT, y=fd)) + ylim(0,1) + scale_x_log10(limits=c(10, 1000))

ggplot(T2MZP_data) +
  geom_point(aes(x=age.at.S1K, y=donor_ki67 ), col=2, size=2) +
  geom_point(aes(x=age.at.S1K, y=total_counts * fd * donor_ki67), col=4) + scale_x_log10(limits=c(50, 1000)) 

ggplot(T2MZP_data) +
  geom_point(aes(x=age.at.S1K, y=host_ki67), col=2)+
  geom_point(aes(x=age.at.S1K, y=donor_ki67), col=4) + scale_x_log10(limits=c(50, 1000)) + ylim(0,0.6)



### custom function
theta_spline <- function(Time, nu, theta0){
  t0 = 40
  Theta = exp(theta0) * exp(-nu * (Time-t0))
  return(Theta)
}

chi_spline <- function(Time, chiEst, qEst){
  Chi = ifelse((Time - 10) < 0, 0,
         chiEst * (1-exp(-qEst * (Time - 10))))
  return(Chi)
}

ki_spline <- function(Time, eps_0, eps_f, A){
  Ki_val = exp(- eps_f * (Time + A)) + eps_0
  return(Ki_val)
}

## prediction plots
ts_p <- seq(50, 750)


T1_file <- file.path("datafiles", "T1_as_sourceSpline.csv")
T1_data <- read.csv(T1_file) %>% arrange(age.at.S1K)

T2_file <- file.path("datafiles", "T2_as_sourceSpline.csv")
T2_data <- read.csv(T2_file) %>% arrange(age.at.S1K)

FM_file <- file.path("datafiles", "FM_as_sourceSpline.csv")
FM_data <- read.csv(FM_file) %>% arrange(age.at.S1K)

## counts fit
T1_counts_nlm <- nls(total_counts ~ theta_spline(age.at.S1K, nu, theta0),
                  data =  T1_data, start = list(nu=0.005, theta0=14))
T1_counts_pars <- coef(T1_counts_nlm)

T1_counts_fit <- data.frame(ts_p, "y_sp" = theta_spline(ts_p, T1_counts_pars[1], T1_counts_pars[2]))

ggplot() + geom_point(data=T1_data, aes(age.at.S1K, total_counts), col=4, size =2) +
  geom_line(data = T1_counts_fit, aes(x = ts_p, y = y_sp), col=4, size =1) + 
  scale_y_log10(limits=c(1e5, 1e7))+ 
  scale_x_log10(limits= c(50, 750), breaks = c(75, 150, 300,  600))+
  labs(title = 'Counts of T1 cells',  y=NULL,  x = 'Host age (days)') 

T2_counts_nlm <- nls(total_counts ~ theta_spline(age.at.S1K, nu, theta0),
                     data =  T2_data, start = list(nu=0.005, theta0=14))
T2_counts_pars <- coef(T2_counts_nlm)

T2_counts_fit <- data.frame(ts_p, "y_sp" = theta_spline(ts_p, T2_counts_pars[1], T2_counts_pars[2]))

ggplot() + geom_point(data=T2_data, aes(age.at.S1K, total_counts), col=4, size =2) +
  geom_line(data = T2_counts_fit, aes(x = ts_p, y = y_sp), col=4, size =1) + 
  scale_y_log10(limits=c(1e5, 1e7))+ 
  scale_x_log10(limits= c(50, 750), breaks = c(75, 150, 300,  600))+
  labs(title = 'Counts of T1 cells',  y=NULL,  x = 'Host age (days)') 

FM_counts_nlm <- glm(log(total_counts) ~ age.at.S1K,
                     data =  FM_data)
FM_counts_pars <- coef(FM_counts_nlm)

FM_counts_fit <- data.frame(ts_p, "y_sp" = exp(FM_counts_pars[1]) * exp(0.002* ts_p))
FM_counts_fit_mm <- data.frame(ts_p, "y_sp" = theta_spline(ts_p, -0.0015, FM_counts_pars[2]))

ggplot() + geom_point(data=FM_data, aes(age.at.S1K, total_counts), col=4, size =2) +
  geom_line(data = FM_counts_fit, aes(x = ts_p, y = y_sp), col=4, size =1) + 
  scale_y_log10(limits=c(1e5, 1e8))+ 
  scale_x_log10(limits= c(50, 750), breaks = c(75, 150, 300,  600))+
  labs(title = 'Counts of T1 cells',  y=NULL,  x = 'Host age (days)') 


T2MZP_counts_nlm <- nls(total_counts ~ theta_spline(age.at.S1K, nu, theta0),
                  data =  T2MZP_data, start = list(nu=0.005, theta0=10))
T2MZP_counts_pars <- coef(T2MZP_counts_nlm)

T2MZP_counts_fit <- data.frame(ts_p, "y_sp" = theta_spline(ts_p, 0, T2MZP_counts_pars[2]))

ggplot() + geom_point(data=T2MZP_data, aes(age.at.S1K, total_counts), col=4, size =2) +
  geom_line(data = T2MZP_counts_fit, aes(x = ts_p, y = y_sp), col=4, size =1) + 
  scale_y_log10(limits=c(1e3, 5e4))+ 
  scale_x_log10(limits= c(50, 750), breaks = c(75, 150, 300,  600))+
  labs(title = 'Counts of T2-MZP cells',  y=NULL,  x = 'Host age (days)') 


## fd fit
ts_dpt <- seq(15, 600)

T1_chi_nlm <- nls(fd ~ chi_spline(age.at.S1K - age.at.bmt, chiEst, qEst),
                     data =  T1_data, start = list(chiEst=0.84, qEst=0.1))
T1_chi_pars <- coef(T1_chi_nlm)

fd_fit_T1 <- data.frame(ts_dpt, "y_sp" = chi_spline(ts_dpt, T1_chi_pars[1], T1_chi_pars[2]))
fd_fit_T1 <- data.frame(ts_dpt, "y_sp" = chi_spline(ts_dpt, 0.75, 0.1748))

ggplot() + 
  geom_point(data=T1_data, aes(age.at.S1K-age.at.bmt, fd), col=4, size =2) +
  geom_line(data = fd_fit_T1, aes(x = ts_dpt , y = y_sp), col=4, size =1) + 
  ylim(0,1) + scale_x_log10(limits=c(10, 750)) +
  labs(title = 'Donor fraction in T1 cells',  y=NULL,  x = 'Days post BMT') 


T2_chi_nlm <- nls(fd ~ chi_spline(age.at.S1K - age.at.bmt, chiEst, qEst),
                  data =  T2_data, start = list(chiEst=0.84, qEst=0.1))
T2_chi_pars <- coef(T2_chi_nlm)

fd_fit_T2 <- data.frame(ts_dpt, "y_sp" = chi_spline(ts_dpt, T2_chi_pars[1], T2_chi_pars[2]))
fd_fit_T2 <- data.frame(ts_dpt, "y_sp" = chi_spline(ts_dpt, 0.75, 0.1383))

ggplot() + 
  geom_point(data=T2_data, aes(age.at.S1K-age.at.bmt, fd), col=4, size =2) +
  geom_line(data = fd_fit_T2, aes(x = ts_dpt , y = y_sp), col=4, size =1) + 
  ylim(0,1) + scale_x_log10(limits=c(10, 750)) +
  labs(title = 'Donor fraction in T2 cells',  y=NULL,  x = 'Days post BMT') 

FM_chi_nlm <- nls(fd ~ chi_spline(age.at.S1K - age.at.bmt, chiEst, qEst),
                  data =  FM_data, start = list(chiEst=0.84, qEst=0.1))
FM_chi_pars <- coef(FM_chi_nlm)

fd_fit_FM <- data.frame(ts_dpt, "y_sp" = chi_spline(ts_dpt, FM_chi_pars[1], FM_chi_pars[2]))
fd_fit_FM <- data.frame(ts_dpt, "y_sp" = chi_spline(ts_dpt, 0.8, 0.0372))

ggplot() + 
  geom_point(data=FM_data, aes(age.at.S1K-age.at.bmt, fd), col=4, size =2) +
  geom_line(data = fd_fit_FM, aes(x = ts_dpt , y = y_sp), col=4, size =1) + 
  ylim(0,1) + scale_x_log10(limits=c(10, 750)) +
  labs(title = 'Donor fraction in T2 cells',  y=NULL,  x = 'Days post BMT') 

T2MZP_chi_nlm <- nls(fd ~ chi_spline(age.at.S1K - age.at.BMT, chiEst, qEst),
               data =  T2MZP_data, start = list(chiEst=0.8, qEst=0.01))
T2MZP_chi_pars <- coef(T2MZP_chi_nlm)

fd_fit <- data.frame(ts_dpt, "y_sp" = chi_spline(ts_dpt, T2MZP_chi_pars[1], T2MZP_chi_pars[2]))
fd_fit_mm <- data.frame(ts_dpt, "y_sp" = chi_spline(ts_dpt, 0.75, 0.011))

ggplot() + 
  geom_point(data=T2MZP_data, aes(age.at.S1K-age.at.BMT, fd), col=4, size =2) +
  geom_line(data = fd_fit, aes(x = ts_dpt , y = y_sp), col=4, size =1) + 
  geom_line(data = fd_fit_mm, aes(x = ts_dpt , y = y_sp), col=4, size =1) + 
  ylim(0,1) + scale_x_log10(limits=c(10, 750)) +
  labs(title = 'Donor fraction in T2-MZP cells',  y=NULL,  x = 'Days post BMT') 


## ki67 fits
T2_donor_ki_nlm <- glm(donor_ki67_Tra ~ age.at.S1K,
                          data = T2_data)

T2_ki_pars_donor <- coef(T2_donor_ki_nlm)

T2_host_ki_nlm <- glm(host_ki67_Tra ~ age.at.S1K,
                      data = T2_data)

T2_ki_pars_host <- coef(T2_host_ki_nlm)

T2_ki_fit_donor <- data.frame(ts_p, "y_sp" = ts_p * T2_ki_pars_donor[2] + T2_ki_pars_donor[1])
T2_ki_fit_host <- data.frame(ts_p, "y_sp" =  ts_p * T2_ki_pars_host[2] + T2_ki_pars_host[1])

ggplot() + geom_point(data=T2_data, aes(age.at.S1K, donor_ki67_Tra), col=4, size =2) +
  geom_point(data=T2_data, aes(age.at.S1K, host_ki67_Tra), col=2, size =2) + 
  geom_line(data = T2_ki_fit_donor, aes(x = ts_p, y = y_sp), col=4, size =1) + 
  geom_line(data = T2_ki_fit_host, aes(x = ts_p, y = y_sp), col=2, size =1) + 
  scale_x_log10(limits= c(50, 750), breaks = c(75, 150, 300,  600)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = 'Ki67+ proprtions in T2-MZP cells',  y=NULL,  x = 'Host age (days)') 


FM_donor_ki_nlm <- nls(donor_ki67_source ~ ki_spline(age.at.S1K, eps_0, eps_f, A),
                      data = FM_data, start = list(eps_0 = 0.1, eps_f = 0.03, A=3))

FM_ki_pars_donor <- coef(FM_donor_ki_nlm)

FM_host_ki_nlm <- glm(host_ki67_source ~ age.at.S1K,
                      data = FM_data)

FM_ki_pars_host <- coef(FM_host_ki_nlm)

FM_ki_fit_donor <- data.frame(ts_p, "y_sp" = ki_spline(ts_p, 0.17, 0.1, -70))
FM_ki_fit_host <- data.frame(ts_p, "y_sp" =  ts_p * 0 + 0.17)

ggplot() + geom_point(data=FM_data, aes(age.at.S1K, donor_ki67_source), col=4, size =2) +
  geom_point(data=FM_data, aes(age.at.S1K, host_ki67_source), col=2, size =2) + 
  geom_line(data = FM_ki_fit_donor, aes(x = ts_p, y = y_sp), col=4, size =1) + 
  geom_line(data = FM_ki_fit_host, aes(x = ts_p, y = y_sp), col=2, size =1) + 
  scale_x_log10(limits= c(50, 750), breaks = c(75, 150, 300,  600)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = 'Ki67+ proprtions in FM B cells',  y=NULL,  x = 'Host age (days)') 


T2MZP_donor_ki_nlm <- nls(donor_ki67 ~ ki_spline(age.at.S1K, eps_0, eps_f, A),
                    data = T2MZP_data, start = list(eps_0 = 0.1, eps_f = 0.03, A=3))

T2MZP_ki_pars_donor <- coef(T2MZP_donor_ki_nlm)

T2MZP_host_ki_nlm <- nls(host_ki67 ~ ki_spline(age.at.S1K, eps_0, eps_f, A),
                   data = T2MZP_data, start = list(eps_0 = 0.1, eps_f = 0.03, A=3))

T2MZP_ki_pars_host <- coef(T2MZP_host_ki_nlm)


T2MZP_ki_fit_donor <- data.frame(ts_p, "y_sp" = ki_spline(ts_p, T2MZP_ki_pars_donor[1], T2MZP_ki_pars_donor[2], T2MZP_ki_pars_donor[3]))
T2MZP_ki_fit_donor_mm <- data.frame(ts_p, "y_sp" = ki_spline(ts_p, 0.068, 0.01, 115))
T2MZP_ki_fit_host <- data.frame(ts_p, "y_sp" = ki_spline(ts_p, T2MZP_ki_pars_host[1], 0, 0))

ggplot() + geom_point(data=T2MZP_data, aes(age.at.S1K, donor_ki67), col=4, size =2) +
  geom_point(data=T2MZP_data, aes(age.at.S1K, host_ki67), col=2, size =2) + 
  #geom_line(data = ki_fit_donor, aes(x = ts_p, y = y_sp), col=6, size =1) + 
  geom_line(data = T2MZP_ki_fit_donor_mm, aes(x = ts_p, y = y_sp), col=4, size =1) + 
  geom_line(aes(x = ts_p, y = rep(0.068, length(ts_p))), col=2, size =1) + 
  scale_x_log10(limits= c(50, 750), breaks = c(75, 150, 300,  600)) +
  scale_y_continuous(limits = c(0, 0.46)) +
  labs(title = 'Ki67+ proprtions in T2-MZP cells',  y=NULL,  x = 'Host age (days)') 










