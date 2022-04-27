## clearing the environment and freeing the memory held by previously loaded objects
rm(list = ls()); gc()   

require(tidyverse)
require(rstan)

expose_stan_functions("stan_models/MAP_neutral_cd4.stan")

init_cond <- c(0, 0, exp(10) * 0.2, exp(10) * 0.8) 
parms <- c(0.3, 0.05, 0.004, 4, exp(10), 0.2)
theta_spline(45)

chivec <- sapply(ts_pred1, Chi_T2MZP_timevar)
Chi_T2MZP_timevar(42)

eps_donor_timevar(45)

solve_shm(100, init_cond, parms)

solve_chi(46, 45, init_cond, parms)


# time sequence for predictions specific to age bins within the data
ts_pred1 <- 10^seq(log10(45), log10(750), length.out = 300)
ts_pred2 <- 10^seq(log10(67), log10(750), length.out = 300)
ts_pred3 <- 10^seq(log10(89), log10(750), length.out = 300)
tb_pred1 <- rep(40, 300)
tb_pred2 <- rep(65, 300)
tb_pred3 <- rep(85, 300)

solve_ode_chi(50, 45, init_cond, parms)

# tb_time_pred -- for solving for initial conditions at each tb (here the mean age of each age bin)
tb_time_pred1 <- c(41, Mean_agebin1)
tb_time_pred2 <- c(41, Mean_agebin2)
tb_time_pred3 <- c(41, Mean_agebin3)

num_pred1 = length(ts_pred1)
num_pred2 = length(ts_pred2)
num_pred3 = length(ts_pred3)

chi_spline <- Vectorize(Chi_spline)



init_cond <- c(0, 0, 2e6, 2e7)
parms <- c(psi = 0.2, rho = 0.001, delta = 0.01, Beta = 4, r_psi = 0.001) 
rdata <- c(theta0 =  exp(14.38), nu = 0.002, chiEst = 0.75, qEst = 0.12, eps_donor = 0.99, eps_host = 0.97)

#predictions
ode_df1 <- solve_ode_pred(ts_pred1, init_cond, parms, rdata, tb_pred1, num_pred1, tb_time_pred1)
ode_df2 <- solve_ode_pred(ts_pred2, init_cond, parms, rdata, tb_pred2, num_pred2, tb_time_pred2)
ode_df3 <- solve_ode_pred(ts_pred3, init_cond, parms, rdata, tb_pred3, num_pred3, tb_time_pred3)

stan_pred_df1 <- data.frame("time" = ts_pred1,
                           "y_pred" = matrix(unlist(ode_df1), nrow = length(ode_df1), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4,
         fd = (y_pred.1 + y_pred.2)/ (total_counts * chi_spline(ts_pred1, rdata["chiEst"], rdata["qEst"])),
         host_ki = y_pred.3/(y_pred.3 + y_pred.4),
         donor_ki = y_pred.1/(y_pred.1 + y_pred.2))

stan_pred_df2 <- data.frame("time" = ts_pred2,
                            "y_pred" = matrix(unlist(ode_df2), nrow = length(ode_df2), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4,
         fd = (y_pred.1 + y_pred.2)/ (total_counts * chi_spline(ts_pred2, rdata["chiEst"], rdata["qEst"])),
         host_ki = y_pred.3/(y_pred.3 + y_pred.4),
         donor_ki = y_pred.1/(y_pred.1 + y_pred.2))

stan_pred_df3 <- data.frame("time" = ts_pred3,
                            "y_pred" = matrix(unlist(ode_df3), nrow = length(ode_df3), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4,
         fd = (y_pred.1 + y_pred.2)/ (total_counts * chi_spline(ts_pred3, rdata["chiEst"], rdata["qEst"])),
         host_ki = y_pred.3/(y_pred.3 + y_pred.4),
         donor_ki = y_pred.1/(y_pred.1 + y_pred.2))


ggplot()+
  geom_point(data = stan_pred_df1, aes(x = time, y= total_counts), size = 2) +
  geom_point(data = stan_pred_df2, aes(x = time, y= total_counts), col =2, size = 2) +
  geom_point(data = stan_pred_df3, aes(x = time, y= total_counts), col =3,  size = 2) + scale_y_log10(limits = c(1e6, 1e8))


ggplot()+
  geom_point(data = stan_pred_df1, aes(x = time, y= fd), size = 2) +
  geom_point(data = stan_pred_df2, aes(x = time, y= fd), col =2, size = 2) +
  geom_point(data = stan_pred_df3, aes(x = time, y= fd), col =3,  size = 2) +
  ylim(0, 1)


####
## ODE sol in R
tb = 40;
eps_donor = 0.95;
eps_host = 0.95;

theta_funct <- function(Time){
   exp(9.495) * exp(-0.01 * (Time))
}

integran_func <- function(Time,lambda){
  theta_funct(Time) * exp(lambda * Time)
}

integ_func <- function(Time, lambda){
  integrate(integran_func, 0, Time, lambda=lambda)$value
}

integ_func(100, 0.04)

ode_sol <- function(Time, lambda, y0){
  
  #y0 * exp(-lambda * Time) + (exp(-lambda * Time) * (Theta0/(lambda - nu)) * ((exp(lambda - nu) * Time) - 1))
  y0 * exp(-lambda * Time) + exp(-lambda * Time) * integ_func(Time, lambda)
}

ode_vec <- Vectorize(ode_sol)

sol_time <- seq(40, 300, length.out = 20)

ode_vec(sol_time, 0.04, 1000)


require(rstan)

expose_stan_functions("stan_models/ki67_SHM_T2MZP.stan")
sapply(sol_time, theta_spline)
sapply(sol_time - 45, Chi_T2MZP_timevar)
sapply(sol_time, eps_donor_timevar)

out <- sapply(sol_time, solve_shm, init_cond = c(0, 0, 4000, 1000), parms = c(0.01, 0.05, 0.03, 4, 45))
out

simple_mod <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    
    dy1 <- psi * theta_funct(t) - lambda * y1
    
    dy2 <- (1-psi) * theta_funct(t) - lambda * y2
    
    return(list(c(y1, y2)))
  })
}

params <- c(psi=0,03, lambda = 0.04)
times <- c(sol_time)
state <- c(y1=0, y2=1000)
out <- deSolve::ode(y=state, times=times, func=simple_mod, parms = params)
out



neutral_mod <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    dY1 <- (psi * theta_spline(Time) * Chi_T1(Time - tb) * eps_donor) +
      rho * (2*Y2 + Y1) - (beta  + delta) * Y1
    
    dY2 <- (psi * theta_spline(Time) * Chi_T1(Time - tb) * (1 - eps_donor)) +
      beta * Y1 - (rho + delta) * Y2
    
    dY3 <- (psi * theta_spline(Time) * (1 - Chi_T1(Time - tb)) * eps_host) +
      rho * (2*Y4 + Y3) - (beta  + delta) * Y3
    
    dY2 <- (psi * theta_spline(Time) * (1 - Chi_T1(Time - tb)) * (1 - eps_host)) +
      beta * Y3 - (rho + delta) * Y4
    
    return(list(c(Y1, Y2, Y3, Y4)))
    })
}

pars <- c(psi = 0.3, rho = 0.005, delta = 0.04, beta = 1/4)
yini <- c(Y1 = 1000, Y2 = 1000, Y3 =exp(10) * 0.2, Y4= exp(10) * 0.8 )

solve_time <- seq(40,  300, length.out = 100)

out <- deSolve::ode(yini, solve_time, neutral_mod, pars)
out[1:10, 1:5]

###################################################################################################################################################
###################################################################################################################################################











