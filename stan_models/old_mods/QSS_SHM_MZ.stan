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
    real[] shm(real t,  real[] k, real[] parms, real[] rdata, int[] idata) {

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
    real y_hat[num_index, 4];       // the array for ODE solution
    real y0[num_tb, 4];             // the array for ode solution of initial conditions
    real init_tb[4];                // the array that defines initial conditions

    //initial conidtion at each tb
    y0 = solve_init(tb_time, init_cond, parms, rdata, 40, num_tb);

    init_tb[1] = init_cond[1];                                           //at tbmt - donor ki67Hi subset size is zero
    init_tb[2] = init_cond[2];                                           //at tbmt - donor ki67Lo subset size is zero

    for (i in 1:num_index){
      init_tb[3] = y0[tb_index[i], 1] + y0[tb_index[i], 3];              //at tbmt - all ki67Hi cells would be host
      init_tb[4] = y0[tb_index[i], 2] + y0[tb_index[i], 4];              //at tbmt - all ki67Lo cells would be host
      y_hat[i] = foreach_ode(solve_time[i], tb[i], init_tb, parms, rdata, tb[i]);
    }

    return y_hat;
   }

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
  }

transformed data{
  real y1[numObs];
  real y2[numObs];
  real y3[numObs];
  real y4[numObs];
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
  y3 = logit_boundary_array(ki_donor);             // transforming counts of ki67 positive cells in the  donor compartments to feed in to ODEs
  y4 = logit_boundary_array(ki_host);              // transforming counts of ki67 positive cells in the  host compartments to feed in to ODEs
}

parameters{
  // parameters to sample with boundary conditions
  //real y0_Log;
  real<lower = 0, upper=1> kappa_0;
  real<lower = 0, upper=1> psi;
  real<lower = 0> rho;
  real<lower = 0> lambda;
  real<lower = 0> Beta;

  // stdev within individual datasets to be estimated
  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real<lower = 0> sigma3;
  real<lower = 0> sigma4;
  }

transformed parameters{
  real k_hat[num_index, 4];  // declaring tge array for ODE solution
  real y1_mean[numObs];      // predictions for dataset1 based on ODE solution
  real y2_mean[numObs];      // predictions for dataset2 based on ODE solution
  real y3_mean[numObs];      // predictions for dataset3 based on ODE solution
  real y4_mean[numObs];      // predictions for dataset4 based on ODE solution

  real parms[4];              // declaring the array for parameters
  real init_cond[4];          // declaring the array for state variables

  // steady state assumption allows us to estimate the constant population size
  real y0 = (psi * theta0) / lambda;

  // initial conditions and parameters
  // kapp0 = ki67 proportions within host at t0 estimated from the model fit
  // at t0 donor counts i.e Nd_0 = 0 --> Nh_0 = N0
  // also kappa_host_0 = kapp0
  init_cond[1] = Nd_0;
  init_cond[2] = Nd_0;
  init_cond[3] = y0 * kappa_0;
  init_cond[4] = y0 * (1 - kappa_0);
  parms[1] = psi;
  parms[2] = rho;
  parms[3] = Beta;
  parms[4] = lambda;

  // solution of the system of ODEs for the predictor values
  k_hat = solve_ode(solve_time, init_cond, parms, rdata, tb, num_index, tb_time, tb_index, num_tb);

  for (i in 1:numObs){
    // total counts
    y1_mean[i] = k_hat[time_index[i], 1] + k_hat[time_index[i], 2] + k_hat[time_index[i], 3] + k_hat[time_index[i], 4];

    // donor fractions normalised with chimerism in the source
    y2_mean[i] = (k_hat[time_index[i], 1] + k_hat[time_index[i], 2])/(y1_mean[i] * Chi_spline(dpBMT[i], chiEst, qEst));

    // fractions of ki67 positive cells in the donor compartment
    y3_mean[i] = (k_hat[time_index[i], 1])/(k_hat[time_index[i], 1] + k_hat[time_index[i], 2]);

    // fractions of ki67 positive cells in the host compartment
    y4_mean[i] = (k_hat[time_index[i], 3])/(k_hat[time_index[i], 3] + k_hat[time_index[i], 4]);
  }
}

model{
  // prior distribution for model parameters
  psi ~ normal(0.5, 0.25);
  kappa_0 ~  normal(0.1, 0.15);
  rho ~ normal(0.01, 0.1);
  lambda ~ normal(0.05, 0.05);
  Beta ~ normal(4, 1.2);

  sigma1 ~ normal(0.6, 0.2);
  sigma2 ~ normal(3, 0.5);
  sigma3 ~ normal(0.3, 0.2);
  sigma4 ~ normal(0.9, 2);

  // model fitting on to data
  y1 ~ normal(log(y1_mean), sigma1);
  y2 ~ normal(logit_boundary_array(y2_mean), sigma2);
  y3 ~ normal(logit_boundary_array(y3_mean), sigma3);
  y4 ~ normal(logit_boundary_array(y4_mean), sigma4);
}

generated quantities{
  // ODE predictions for 3 different timecourses based on the age at BMT bins
  real k_hat_pred_age1[numPred1, 4];
  real k_hat_pred_age2[numPred2, 4];
  real k_hat_pred_age3[numPred3, 4];
  // model predictions for the datasets for 3 different timecourses based on the age at BMT bins
  real y1_mean_pred_age1[numPred1];
  real y2_mean_pred_age1[numPred1];
  real y1_mean_pred_age2[numPred2];
  real y2_mean_pred_age2[numPred2];
  real y1_mean_pred_age3[numPred3];
  real y2_mean_pred_age3[numPred3];
  real y3_mean_pred1[numPred1];
  real y4_mean_pred1[numPred1];
  real y3_mean_pred2[numPred2];
  real y4_mean_pred2[numPred2];
  real y3_mean_pred3[numPred3];
  real y4_mean_pred3[numPred3];
  // model predictions for the datasets with stdev estimated from the fits for 3 different timecourses based on the age at BMT bins
  real countspred_age1[numPred1];
  real fdpred_age1[numPred1];
  real countspred_age2[numPred2];
  real fdpred_age2[numPred2];
  real countspred_age3[numPred3];
  real fdpred_age3[numPred3];
  real donor_kiprop_pred1[numPred1];
  real host_kiprop_pred1[numPred1];
  real donor_kiprop_pred2[numPred2];
  real host_kiprop_pred2[numPred2];
  real donor_kiprop_pred3[numPred3];
  real host_kiprop_pred3[numPred3];
  // log likelihood for individual data sets and combined
  vector[numObs] log_lik;
  vector[numObs] log_lik1;
  vector[numObs] log_lik2;
  vector[numObs] log_lik3;
  vector[numObs] log_lik4;

  // parameters of interest
  real lambda_inv = 1/lambda;
  real rho_inv = 1/rho;
  real delta = lambda + rho;
  real delta_inv = 1/delta;

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

  // at t0 the ki67 proportions in donor of target = ki67 proportions in source donor compartemnt
  y3_mean_pred1[1] = eps_donor;
  donor_kiprop_pred1[1] = expit(normal_rng(logit_boundary(y3_mean_pred1[1]), sigma3));
  y3_mean_pred2[1] = eps_donor;
  donor_kiprop_pred2[1] = expit(normal_rng(logit_boundary(y3_mean_pred2[1]), sigma3));
  y3_mean_pred3[1] = eps_donor;
  donor_kiprop_pred3[1] = expit(normal_rng(logit_boundary(y3_mean_pred3[1]), sigma3));

  // ki67 proportions at t0 estimated from the model fit
  y4_mean_pred1[1] = kappa_0;
  host_kiprop_pred1[1] = expit(normal_rng(logit_boundary(y4_mean_pred1[1]), sigma4));
  y4_mean_pred2[1] = kappa_0;
  host_kiprop_pred2[1] = expit(normal_rng(logit_boundary(y4_mean_pred2[1]), sigma4));
  y4_mean_pred3[1] = kappa_0;
  host_kiprop_pred3[1] = expit(normal_rng(logit_boundary(y4_mean_pred3[1]), sigma4));

  for (i in 2:numPred1){
    y2_mean_pred_age1[i] = (k_hat_pred_age1[i, 1] + k_hat_pred_age1[i, 2])/(y1_mean_pred_age1[i] * Chi_spline(ts_pred1[i] - tb_pred1[i], chiEst, qEst));
    fdpred_age1[i] = expit_boundary(normal_rng(logit_boundary(y2_mean_pred_age1[i]), sigma2));

    y3_mean_pred1[i] = k_hat_pred_age1[i, 1]/(k_hat_pred_age1[i, 1] + k_hat_pred_age1[i, 2]);
    donor_kiprop_pred1[i] =  expit_boundary(normal_rng(logit_boundary(y3_mean_pred1[i]), sigma3));

    y4_mean_pred1[i] = k_hat_pred_age1[i, 3]/(k_hat_pred_age1[i, 3] + k_hat_pred_age1[i, 4]);
    host_kiprop_pred1[i] = expit_boundary(normal_rng(logit_boundary(y4_mean_pred1[i]), sigma4));
  }

  for (i in 2:numPred2){
    y2_mean_pred_age2[i] = (k_hat_pred_age2[i, 1] + k_hat_pred_age2[i, 2])/(y1_mean_pred_age2[i] * Chi_spline(ts_pred2[i] - tb_pred2[i], chiEst, qEst));
    fdpred_age2[i] = expit_boundary(normal_rng(logit_boundary(y2_mean_pred_age2[i]), sigma2));

    y3_mean_pred2[i] = k_hat_pred_age2[i, 1]/(k_hat_pred_age2[i, 1] + k_hat_pred_age2[i, 2]);
    donor_kiprop_pred2[i] =  expit_boundary(normal_rng(logit_boundary(y3_mean_pred2[i]), sigma3));

    y4_mean_pred2[i] = k_hat_pred_age2[i, 3]/(k_hat_pred_age2[i, 3] + k_hat_pred_age2[i, 4]);
    host_kiprop_pred2[i] = expit_boundary(normal_rng(logit_boundary(y4_mean_pred2[i]), sigma4));
  }

  for (i in 2:numPred3){
    y2_mean_pred_age3[i] = (k_hat_pred_age3[i, 1] + k_hat_pred_age3[i, 2])/(y1_mean_pred_age3[i] * Chi_spline(ts_pred3[i] - tb_pred3[i], chiEst, qEst));
    fdpred_age3[i] = expit_boundary(normal_rng(logit_boundary(y2_mean_pred_age3[i]), sigma2));

    y3_mean_pred3[i] = k_hat_pred_age3[i, 1]/(k_hat_pred_age3[i, 1] + k_hat_pred_age3[i, 2]);
    donor_kiprop_pred3[i] =  expit_boundary(normal_rng(logit_boundary(y3_mean_pred3[i]), sigma3));

    y4_mean_pred3[i] = k_hat_pred_age3[i, 3]/(k_hat_pred_age3[i, 3] + k_hat_pred_age3[i, 4]);
    host_kiprop_pred3[i] = expit_boundary(normal_rng(logit_boundary(y4_mean_pred3[i]), sigma4));
  }


  // calculating the log predictive accuracy for each point
  for (n in 1:numObs) {
    log_lik1[n] = normal_lpdf(y1[n] | log(y1_mean[n]), sigma1);
    log_lik2[n] = normal_lpdf(y2[n] | logit_boundary(y2_mean[n]), sigma2);
    log_lik3[n] = normal_lpdf(y3[n] | logit_boundary(y3_mean[n]), sigma3);
    log_lik4[n] = normal_lpdf(y4[n] | logit_boundary(y4_mean[n]), sigma4);
    log_lik[n] = log_lik1[n] + log_lik2[n] + log_lik3[n] + log_lik4[n];
  }
}
