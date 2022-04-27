functions{
  // function that describes the changes in source counts
  real theta_spline(real t,
    real nu,               // nu is the rate of decline of source influx assuming that its proportional to the counts of source compartment
    real spline_int) {     // spline_int gives the initial counts of the source compartment

      real theta;
      int ta = 40;         // the earliest age at BMT

      theta = spline_int * exp(-nu * (t-ta));
      return theta;
   }

   // function that describes changes in source chimerism
   real Chi_spline( real t,
     real chiEst,         // chiEst is the level if stabilised chimerism in the source compartment
     real qEst) {         // qEst is the rate with which cimerism chnages in the source compartment

       real chi;

       if (t < 0){
         chi = 0;                       // conditioning the function to adapt to the timepoints before BMT
       } else {
         chi = chiEst * (1 - exp(-qEst * t));
       }

         return chi;
    }


// function that contains the ODE equations to be used in ODE solver
  real[] shm(real t, real[] y, real[] parms, real[] rdata, int[] idata) {
    real dydt[2];             // system of ODEs

    // the array of parameters invloved in ode equations
    real psi    = parms[1];
    real lambda = parms[2];

    // the array of input data (fixed varibaled) invloved in ode equations
    real theta0 = rdata[1];
    real nu     = rdata[2];
    real chiEst = rdata[3];
    real qEst   = rdata[4];

    // for age correction -- 'tb' is the age of host at BMT
     real tb = idata[1];

    dydt[1] = (psi * theta_spline(t, nu, theta0) * Chi_spline(t - tb, chiEst, qEst)) - lambda * y[1];
    dydt[2] = (psi * theta_spline(t, nu, theta0) * (1 - Chi_spline(t - tb, chiEst, qEst))) - lambda * y[2];

    return dydt;
  }

  // function that contains the ODE equations to be used in ODE solver
  real[] shm_pred(real t,  real[] k, real[] parms, real[] rdata, int[] idata) {

   // the array of parameters invloved in ode equations
   real psi    = parms[1];
   real rho    = parms[2];
   real Beta   = parms[3];
   real lambda = parms[4];

   // the array of input data (fixed varibaled) invloved in ode equations
   real theta0    = rdata[1];
   real nu        = rdata[2];
   real chiEst    = rdata[3];
   real qEst      = rdata[4];
   real eps_host  = rdata[5];
   real eps_donor = rdata[6];

  // for age correction -- 'tb' is the age of host at BMT
   real tb = idata[1];

   // the system of ODEs
   real dkdt[4];
   dkdt[1] = (psi * theta_spline(t, nu, theta0) * Chi_spline(t-tb, chiEst, qEst) * eps_donor) + rho * (2 * k[2] + k[1]) - ((1/Beta) + (lambda + rho)) * k[1];
   dkdt[2] = (psi * theta_spline(t, nu, theta0) * Chi_spline(t-tb, chiEst, qEst) * (1-eps_donor)) + (1/Beta) * k[1] - (rho + (lambda + rho)) * k[2];

   dkdt[3] = (psi * theta_spline(t, nu, theta0) * (1-Chi_spline(t-tb, chiEst, qEst)) * eps_host) + rho * (2 * k[4] + k[3]) - ((1/Beta) + (lambda + rho)) * k[3];
   dkdt[4] = (psi * theta_spline(t, nu, theta0) * (1-Chi_spline(t-tb, chiEst, qEst)) * (1-eps_host)) + (1/Beta) * k[3] - (rho + (lambda + rho)) * k[4];

   return dkdt;
 }

   real rho_hat(real lambda_est, real kappa_est, real eps_est, real beta) {

   real rho_est;
   real T_est;

   T_est = 1/beta;

   rho_est = (kappa_est + (lambda_est * T_est) - (lambda_est * T_est * eps_est))/ (T_est * (1 - kappa_est) * 2);

   return rho_est;
   }

   real delta_hat(real lambda_est, real kappa_est, real eps_est, real beta) {

   real delta_est;
   real T_est;

   T_est = 1/beta;

   delta_est = (kappa_est + lambda_est * T_est * (2 - eps_est - kappa_est))/ (T_est * (1 - kappa_est) * 2);

   return delta_est;
   }

   // to calculate total counts (N) at each 'tb'
   // this function solves the ODE for each value of 'tb' using 'ta' as the starting point (t0) and converts the solution into a 1d array
   real[] foreach_init(real tb, real ta, real[] init_cond, real[] parms, real[] rdata, int x_i){
     //first argument to integrate_ode_rk45 must be the name of the ODE function with signature (real, real[], real[], real[], int[])
     return to_array_1d(integrate_ode_rk45(shm, init_cond, ta, rep_array(tb, 1), parms, rdata, rep_array(x_i, 1)));
   }

   // ode solution at each 'tb' -- this serves as initial condition when solving at 't'
   real[,] solve_init(real[] tb_time, real[] init_cond, real[] parms, real[] rdata, int x_i, int num_tb){
    real y_init[num_tb, 4];

    y_init[1] = init_cond;
    for (i in 2:num_tb){
      y_init[i] = foreach_init(tb_time[i], tb_time[1], init_cond, parms, rdata, x_i);
    }

    return y_init;
   }

   // this function solves the ODE for all obsreved timepoints using 'tb' as the strating point (t0) and converts the solution into a 1d array
   real[] foreach_ode(real ts, real t0, real[] init_cond, real[] parms, real[] rdata, int x_i) {
     //first argument to integrate_ode_rk45 must be the name of the ODE function with signature (real, real[], real[], real[], int[])
     return to_array_1d(integrate_ode_rk45(shm, init_cond, t0, rep_array(ts, 1), parms, rdata, rep_array(x_i, 1)));
    }

  // ode solution at each 't' i.e. obsreved timepoints
   real[,] solve_ode(real[] solve_time, real[] init_cond, real[] parms, real[] rdata, int[] tb, int num_index, real[] tb_time, int[] tb_index, int num_tb){
    real y_hat[num_index, 2];       // the array for ODE solution
    real y0[num_tb, 2];             // the array for ode solution of initial conditions
    real init_tb[2];                // the array that defines initial conditions

    //initial conidtion at each tb
    y0 = solve_init(tb_time, init_cond, parms, rdata, 40, num_tb);

    init_tb[1] = init_cond[1];                                           //at tbmt - donor ki67Hi subset size is zero
    init_tb[2] = y0[tb_index[i], 1] + y0[tb_index[i], 2];                                      //at tbmt - donor ki67Lo subset size is zero

    for (i in 1:num_index){
        y_hat[i] = foreach_ode(solve_time[i], tb[i], init_tb, parms, rdata, tb[i]);
    }

    return y_hat;
   }

   // ode solution at each 't_pred' i.e. predictor timepoints
   real[,] solve_ode_pred(real[] solve_time, real[] init_cond, real[] parms, real[] rdata, int[] tb, int num_index, real[] tb_time){
    real y_hat[num_index, 4];
    real y0[2, 4];
    real init_tb[4];

    y0 = solve_init(tb_time, init_cond, parms, rdata, 40, 2);

    init_tb[1] = init_cond[1];                                         //at tbmt - donor ki67Lo subset size is zero
    init_tb[2] = init_cond[1];                                         //at tbmt - donor ki67Lo subset size is zero

    init_tb[3] = y0[2, 1] + y0[2, 3];                                  //at tbmt - all ki67Hi cells would be host
    init_tb[4] = y0[2, 2] + y0[2, 4];                                  //at tbmt - all ki67Hi cells would be host

    y_hat[1] = init_tb;
    for (i in 2:num_index){
      y_hat[i] = foreach_ode(solve_time[i], tb[i], init_tb, parms, rdata, tb[i]);
    }

    return y_hat;
   }

  // function to transform arrays using logit function that have values greater than 1
   real[] logit_boundary_array(real[] x){
     int ndims = size(x);   // dimensions of the array to be transformed
     real answer[ndims];
     real b = 1.5;          // the upper boundary for logit transformation

     // logit function
     for (i in 1: ndims){
       answer[i] = log(x[i]/(b-x[i]));
     }
     return answer;
   }

     // function to transform elements/vectors using logit function that have values greater than 1
    real logit_boundary(real x){
     real answer;
     real b = 1.5;          // the upper boundary for logit transformation

      // logit function
     answer = log(x/(b-x));

     return answer;
   }

   // function to caclulate inverse of logit transformed values
   real expit(real x){
     real answer;

    // expit function
     answer = exp(x)/(1+exp(x));

     return answer;
   }

   // function to caclulate inverse of logit transformed values with a boundary greater than 1
   real expit_boundary(real x){
     real answer;
     real b = 1.5;

     // expit function
     answer =  (b * exp(x)) /( 1+ exp(x));

     return answer;
   }
}

data{
  int<lower  = 1> numObs;                           // number of observations for donor fractions may be different that cell counts
  int<lower  = 1> num_index;
  int<lower  = 1> num_tb;
  real<lower = 0> solve_time[num_index];
  int<lower  = 1> time_index[numObs];
  int<lower  = 1> tb_index[num_index];
  int<lower = 0> ageAtBMT[num_index];
  real<lower =0> tb_time[num_tb];
  int dpBMT[numObs];
  real Nd_0;
  real<lower = 0> counts[numObs];
  real<lower = 0> Nfd[numObs];
  real<lower = 0> ki_donor[numObs];
  real<lower = 0> ki_host[numObs];
  int<lower  = 1> numPred1;
  int<lower  = 1> numPred2;
  int<lower  = 1> numPred3;
  real ts_pred1[numPred1];
  real ts_pred2[numPred2];
  real ts_pred3[numPred3];
  int tb_pred1[numPred1];
  int tb_pred2[numPred2];
  int tb_pred3[numPred3];
  real tb_time_pred1[2];
  real tb_time_pred2[2];
  real tb_time_pred3[2];
  real theta0;
  real nu;
  real chiEst;
  real qEst;
  real eps_donor;
  real eps_host;
  real kappa0;
  }

transformed data{
  real y1[numObs];
  real y2[numObs];
  real rdata[6];
  int tb[num_index];

// fixed variables that match to the data inputs in R file
  rdata[1] = theta0;
  rdata[2] = nu;
  rdata[3] = chiEst;
  rdata[4] = qEst;
  rdata[5] = eps_host;
  rdata[6] = eps_donor;

  for (i in 1:num_index){
    tb[i] = ageAtBMT[i];
  }

  y1 = log(counts);                                // transforming cell counts of donor compartments to feed in to ODEs
  y2 = logit_boundary_array(Nfd);                  // transfored donor fractions normalised to source chimerism to feed in to ODEs
}

parameters{
  // parameters to sample with boundary conditions
  real y0_Log;
  real<lower = 0, upper=1> psi;
  real lambda;

  // stdev within individual datasets to be estimated
  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  }

transformed parameters{
  real k_hat[num_index, 2];  // declaring tge array for ODE solution
  real y1_mean[numObs];      // predictions for dataset1 based on ODE solution
  real y2_mean[numObs];      // predictions for dataset2 based on ODE solution

  real parms[2];              // declaring the array for parameters
  real init_cond[2];          // declaring the array for state variables
  real y0 = exp(y0_Log);      // transformed parameters for better/faster sampling

  // initial conditions and parameters
  // kapp0 = ki67 proportions within host at t0 estimated from the model fit
  // at t0 donor counts i.e Nd_0 = 0 --> Nh_0 = N0
  // also kappa_host_0 = kapp0
  init_cond[1] = Nd_0;
  init_cond[2] = y0;

  parms[1] = psi;
  parms[2] = lambda;

  // solution of the system of ODEs for the predictor values
  k_hat = solve_ode(solve_time, init_cond, parms, rdata, tb, num_index, tb_time, tb_index, num_tb);

  for (i in 1:numObs){
    // total counts
    y1_mean[i] = k_hat[time_index[i], 1] + k_hat[time_index[i], 2];

    // donor fractions normalised with chimerism in the source
    y2_mean[i] = (k_hat[time_index[i], 1] + k_hat[time_index[i], 2])/(y1_mean[i] * Chi_spline(dpBMT[i], chiEst, qEst));
  }
}

model{
  // prior distribution for model parameters
  psi ~ normal(0.5, 0.25);
  y0_Log ~ normal(14, 1);
  lambda ~ normal(0.05, 0.05);

  sigma1 ~ normal(0.6, 0.2);
  sigma2 ~ normal(3, 0.5);

  // model fitting on to data
  y1 ~ normal(log(y1_mean), sigma1);
  y2 ~ normal(logit_boundary_array(y2_mean), sigma2);
}

generated quantities{
  // ODE predictions for 3 different timecourses based on the age at BMT bins
  real k_hat_pred_age1[numPred1, 2];
  real k_hat_pred_age2[numPred2, 2];
  real k_hat_pred_age3[numPred3, 2];


  real k_hat_pred[numPred1, 4];
  // model predictions for the datasets for 3 different timecourses based on the age at BMT bins
  real y1_mean_pred_age1[numPred1];
  real y2_mean_pred_age1[numPred1];
  real y1_mean_pred_age2[numPred2];
  real y2_mean_pred_age2[numPred2];
  real y1_mean_pred_age3[numPred3];
  real y2_mean_pred_age3[numPred3];

  // model predictions for the datasets with stdev estimated from the fits for 3 different timecourses based on the age at BMT bins
  real countspred_age1[numPred1];
  real fdpred_age1[numPred1];
  real countspred_age2[numPred2];
  real fdpred_age2[numPred2];
  real countspred_age3[numPred3];
  real fdpred_age3[numPred3];
  real donor_kiprop_pred1[numPred1];
  real host_kiprop_pred1[numPred1];
  //real donor_kiprop_pred2[numPred2];
  //real host_kiprop_pred2[numPred2];
  //real donor_kiprop_pred3[numPred3];
  //real host_kiprop_pred3[numPred3];
  // log likelihood for individual data sets and combined
  vector[numObs] log_lik;
  vector[numObs] log_lik1;
  vector[numObs] log_lik2;

  real init_ki[4];
  real par_ki[4];
  int idata[0];


  // parameters of interest
  real lambda_inv = 1/lambda;
  real rho = rho_hat(lambda, kappa0, eps_host, 0.25);
  real rho_donor = rho_hat(lambda, kappa0, eps_donor, 0.25);
  real delta_donor = delta_hat(lambda, kappa0, eps_donor, 0.25);
  real rho_host = rho_hat(lambda, kappa0, eps_host, 0.25);
  real delta_host = delta_hat(lambda, kappa0, eps_host, 0.25);


  //ODE solution for different age bins
  k_hat_pred_age1 = solve_ode_pred(ts_pred1, init_cond, parms, rdata, tb_pred1, numPred1, tb_time_pred1);
  k_hat_pred_age2 = solve_ode_pred(ts_pred2, init_cond, parms, rdata, tb_pred2, numPred2, tb_time_pred2);
  k_hat_pred_age3 = solve_ode_pred(ts_pred3, init_cond, parms, rdata, tb_pred3, numPred3, tb_time_pred3);

  // Total cell counts (donor + host) for different age bins
  // age bin1
  for (i in 1:numPred1){
    y1_mean_pred_age1[i] = k_hat_pred_age1[i, 1] + k_hat_pred_age1[i, 2] + k_hat_pred_age1[i, 3] + k_hat_pred_age1[i, 4];
    countspred_age1[i] = exp(normal_rng(log(y1_mean_pred_age1[i]), sigma1));
  }
  // age bin2
  for (i in 1: numPred2){
    y1_mean_pred_age2[i] = k_hat_pred_age2[i, 1] + k_hat_pred_age2[i, 2] + k_hat_pred_age2[i, 3] + k_hat_pred_age2[i, 4];
    countspred_age2[i] = exp(normal_rng(log(y1_mean_pred_age2[i]), sigma1));
  }
  // age bin3
  for (i in 1: numPred3){
    y1_mean_pred_age3[i] = k_hat_pred_age3[i, 1] + k_hat_pred_age3[i, 2] + k_hat_pred_age3[i, 3] + k_hat_pred_age3[i, 4];
    countspred_age3[i] = exp(normal_rng(log(y1_mean_pred_age3[i]), sigma1));
  }

  // Initial conditions for fd and frcations of ki67host and ki67 donor
  y2_mean_pred_age1[1] = Nd_0;
  fdpred_age1[1] = Nd_0;
  y2_mean_pred_age2[1] = Nd_0;
  fdpred_age2[1] = Nd_0;
  y2_mean_pred_age3[1] = Nd_0;
  fdpred_age3[1] = Nd_0;


  for (i in 2:numPred1){
    y2_mean_pred_age1[i] = (k_hat_pred_age1[i, 1] + k_hat_pred_age1[i, 2])/(y1_mean_pred_age1[i] * Chi_spline(ts_pred1[i] - tb_pred1[i], chiEst, qEst));
    fdpred_age1[i] = expit_boundary(normal_rng(logit_boundary(y2_mean_pred_age1[i]), sigma2));
  }

  for (i in 2:numPred2){
    y2_mean_pred_age2[i] = (k_hat_pred_age2[i, 1] + k_hat_pred_age2[i, 2])/(y1_mean_pred_age2[i] * Chi_spline(ts_pred2[i] - tb_pred2[i], chiEst, qEst));
    fdpred_age2[i] = expit_boundary(normal_rng(logit_boundary(y2_mean_pred_age2[i]), sigma2));
  }

  for (i in 2:numPred3){
    y2_mean_pred_age3[i] = (k_hat_pred_age3[i, 1] + k_hat_pred_age3[i, 2])/(y1_mean_pred_age3[i] * Chi_spline(ts_pred3[i] - tb_pred3[i], chiEst, qEst));
    fdpred_age3[i] = expit_boundary(normal_rng(logit_boundary(y2_mean_pred_age3[i]), sigma2));
  }

  // calculating the log predictive accuracy for each point
  for (n in 1:numObs) {
    log_lik1[n] = normal_lpdf(y1[n] | log(y1_mean[n]), sigma1);
    log_lik2[n] = normal_lpdf(y2[n] | logit_boundary(y2_mean[n]), sigma2);
    log_lik[n] = log_lik1[n] + log_lik2[n];
  }


  // Predictions for ki67 host and donor data
  // drawing the timecourse of ki67hi proportions from posterior samples
  init_ki[1] = 0;
  init_ki[2] = 0;                                         // initial conditions for host and donor khi and klo counts
  init_ki[3] = y0 * kappa0;
  init_ki[4] = y0 * (1 - kappa0);

  par_ki[1] = psi;
  par_ki[2] = rho;
  par_ki[3] = 4;
  par_ki[4] = lambda;

  k_hat_pred[1, ] = init_ki;
  k_hat_pred[2:numPred1, ] = integrate_ode_rk45(shm_pred, init_ki, solve_time[1], ts_pred1[2:numPred1], par_ki, rdata, idata);

  for (i in 1:numPred1){
    donor_kiprop_pred1[i] = k_hat_pred[i, 1]/(k_hat_pred[i, 1] + k_hat_pred[i, 2]);
    host_kiprop_pred1[i] = k_hat_pred[i, 3]/(k_hat_pred[i, 3] + k_hat_pred[i, 4]);
  }
}
