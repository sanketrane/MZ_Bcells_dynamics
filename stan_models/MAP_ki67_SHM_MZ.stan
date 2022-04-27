functions{
  // function that describes the changes in source counts
  real theta_spline(real time){

      real nu = 1.86e-03; real theta0 = 14.36;

      real theta;
      int t0 = 40;         // the earliest age at BMT

      // nu is the rate of decline of source influx assuming that its proportional to the counts of source compartment
      //theta0 gives the initial counts of the source compartment
      theta = exp(theta0) * exp(-nu * (time-t0));
      return theta;
   }

   real[] thet_s(real[] time){
     int numdim = size(time);
     real answer[numdim];

     for (i in 1:numdim){
       answer[i] = theta_spline(time[i]);
     }

     return answer;
   }

   // function that describes changes in source chimerism
   real Chi_spline( real time){

      real chiEst = 0.76; real qEst = 0.094;

       real chi;
       // chiEst is the level if stabilised chimerism in the source compartment
       // qEst is the rate with which cimerism chnages in the source compartment
       if (time < 0){
         chi = 0;                       // conditioning the function to adapt to the timepoints before BMT
       } else {
         chi = chiEst * (1 - exp(-qEst * time));
       }

         return chi;
    }

    // function that contains the ODE equations to be used in ODE solver
    real[] shm(real time,  real[] k, real[] parms, real[] rdata, int[] idata) {

     // the array of parameters invloved in ode equations
     real psi    = parms[1];
     real rho    = parms[2];
     real lambda = parms[3];
     real Beta   = parms[4];

     real eps_donor = 0.99; real eps_host = 0.95;

     // age of BMT in each recipient
     real ageAtBMT = parms[5];

     // the system of ODEs
     real dkdt[4];

     dkdt[1] = (psi * theta_spline(time) * Chi_spline(time-ageAtBMT) * eps_donor) + rho * (2 * k[2] + k[1]) - ((1/Beta) + (lambda + rho)) * k[1];
     dkdt[2] = (psi * theta_spline(time) * Chi_spline(time-ageAtBMT) * (1-eps_donor)) + (1/Beta) * k[1] - (rho + (lambda + rho)) * k[2];

     dkdt[3] = (psi * theta_spline(time) * (1-Chi_spline(time-ageAtBMT)) * eps_host) + rho * (2 * k[4] + k[3]) - ((1/Beta) + (lambda + rho)) * k[3];
     dkdt[4] = (psi * theta_spline(time) * (1-Chi_spline(time-ageAtBMT)) * (1-eps_host)) + (1/Beta) * k[3] - (rho + (lambda + rho)) * k[4];

     return dkdt;
   }

   real[] solve_shm(real solve_time, real[] init_cond, real[] parms) {
     // solves the ode for each timepoint from t0
     return to_array_1d(integrate_ode_rk45(shm, init_cond, 40.0, rep_array(solve_time, 1), parms, {0.0}, {0}));
    }

   real[] solve_chi(real solve_time,            // time point of observation
     real ageAtBMT,
     real[] init_cond,
     real[] parms){

       real y_solve[4];
       real params[5];

       real y0[4];
       real init_tb[4];                         // init conditions at the mean age of BMT for the group

       params[1:4] = parms[1:4];
       params[5] = ageAtBMT;                                           // age at BMT

       // y0 is the solution for the initial conditions for the age at BMT
       if (ageAtBMT > 40.0){
         y0 = solve_shm(ageAtBMT, init_cond, params);
       } else {
         y0 = init_cond;
       }

       // init conditions at the BMT
       init_tb[1] = 0;                                           //at tbmt - # donor is zero
       init_tb[2] = 0;                                           //at tbmt - # donor is zero
       init_tb[3] = y0[1] + y0[3];                               //at tbmt - all ki67Hi cells are host
       init_tb[4] = y0[2] + y0[4];                               //at tbmt - all ki67Lo cells are host

       y_solve = to_array_1d(integrate_ode_rk45(shm, init_tb, ageAtBMT, rep_array(solve_time, 1), params, {0.0}, {0}));

       return y_solve;
     }

    real[,] solve_ode_chi(real[] solve_time,
     real[] ageAtBMT,
     real[] init_cond,
     real[] parms){

       int numdim = size(solve_time);
       real y_solve[numdim, 4];

       for (i in 1:numdim) {
         y_solve[i] = solve_chi(solve_time[i], ageAtBMT[i], init_cond, parms);
       }

       return y_solve;
     }

    vector math_reduce(vector global_params, vector local_params, real[] x_r, int[] x_i){
      // data for each shard
      int n = size(x_i); // n = 1
      real solve_time = x_r[1];
      int ageAtBMT = x_i[1];                          // time zero -- for chimeras age at BMT
      real tb_time = ageAtBMT/1.0;

      //params
      real N0 = global_params[5];           // cell counts at time zero
      real kappa_0 = global_params[6];            // Ki proportos at time zero
      real init_cond[4];

      // ODE solution -- predictions for the observed timecourse
      real chi_solve[4];
      real chi_counts_mean;
      real host_counts_mean;
      real donor_counts_mean;
      real host_ki_mean;
      real donor_ki_mean;

      vector[4*n] y_mean_stacked;

      // initial conditions
      init_cond[1] = 0.0;
      init_cond[2] = 0.0;
      init_cond[3] = kappa_0 * N0;
      init_cond[4] = (1 - kappa_0) * N0;

      // each shard receives a unique datpoint
      // ODE solution for chimera dataset -- x_r = data time and x_i = time at BMT

      chi_solve = solve_chi(solve_time, tb_time, init_cond, to_array_1d(global_params));
      chi_counts_mean = chi_solve[1] + chi_solve[2] + chi_solve[3] + chi_solve[4];
      donor_counts_mean = chi_solve[1] + chi_solve[2];
      host_counts_mean = chi_solve[3] + chi_solve[4];

      donor_ki_mean = chi_solve[1]/donor_counts_mean;
      host_ki_mean = chi_solve[3]/host_counts_mean;
      y_mean_stacked[1] = chi_counts_mean;
      y_mean_stacked[2] = donor_counts_mean/(chi_counts_mean * Chi_spline(solve_time - tb_time));
      y_mean_stacked[3] = donor_ki_mean;
      y_mean_stacked[4] = host_ki_mean;

      return y_mean_stacked;
    }

    // functions for transformation of fractions in (0,a), where a >=1
    real logit_inverse(real x){
       real ans;

       ans = exp(x)/(1+exp(x));

       return ans;
     }

   // functions for transformation of fractions in (0,a), where a >=1
   real[] asinsqrt_array(real[] x){
     int ndims = size(x);
     real answer[ndims];
     real a = 1.2;

     for (i in 1: ndims){
       answer[i] = asin(sqrt(x[i])/sqrt(a));
     }
     return answer;
   }

   real asinsqrt_real(real x){
     real a = 1.2;

     real answer = asin(sqrt(x)/sqrt(a));
     return answer;
   }

   real asinsqrt_inv(real x){
     real a = 1.2;

     real answer = a * (sin(x))^2;
     return answer;
   }
}

data{
  int n_shards;
  int<lower  = 1> numObs;                           // number of observations for donor fractions may be different that cell counts
  //int<lower  = 1> num_index;
  int<lower  = 1> numPred;
  int<lower  = 1> time_index[numObs];
  real<lower = 0> solve_time[n_shards];
  int ageAtBMT[n_shards];
  real<lower = 0> counts[numObs];
  real<lower = 0> Nfd[numObs];
  real<lower = 0> ki_donor[numObs];
  real<lower = 0> ki_host[numObs];
  real ts_pred1[numPred];
  real ts_pred2[numPred];
  real ts_pred3[numPred];
  real tb_pred1[numPred];
  real tb_pred2[numPred];
  real tb_pred3[numPred];
  }

transformed data{
  int x_i[n_shards, 1];         // each shard gets a single data point
  real x_r[n_shards, 1];        // each shard gets a single data point

  // empty set of per shard params
  vector[0] local_params[n_shards];  // shard specific params --  useful for hierarchical modelling

  // data split into shards
  for (s in 1:n_shards){
   x_i[s, 1] = ageAtBMT[s];                             // age at BMT split
   x_r[s, 1] = solve_time[s];                     // time split
  }
}

parameters{
  // parameters to sample with boundary conditions
  real<lower = 11, upper=17> N0_Log;
  real<lower = 0, upper=1> kappa_0;
  real<lower = 0> psi;
  real<lower = 0> delta;
  real<lower = 0> rho;
  real<lower = 0> Beta;

  // stdev within individual datasets to be estimated
  real<lower=0> sigma_chi_counts;
  real<lower=0> sigma_Nfd;
  real<lower=0> sigma_donor_ki;
  real<lower=0> sigma_host_ki;
}

transformed parameters{
  vector[6] global_params;
  vector[n_shards] y1_mean;                  // PDE prediction for counts from chimera data
  vector[n_shards] y2_mean;                  // PDE prediction for Nfd from chimera data
  vector[n_shards] y3_mean;                  // PDE prediction for ki proportions in donor compartment from chimera data
  vector[n_shards] y4_mean;                  // PDE prediction for ki proportions in host compartment from chimera data
  vector[(4*n_shards)] y_mean_stacked;       // compliled output across all nodes

  real total_counts_mean[numObs];
  real Nfd_mean[numObs];
  real donor_ki_mean[numObs];
  real host_ki_mean[numObs];

  real N0 = exp(N0_Log);

  global_params[1] = psi;
  global_params[2] = delta;
  global_params[3] = rho;
  global_params[4] = Beta;
  global_params[5] = N0;
  global_params[6] = kappa_0;

  // combining the output from all the shards
  y_mean_stacked = map_rect(math_reduce, global_params, local_params, x_r, x_i);

  for (i in 1:n_shards){
    y1_mean[i] = y_mean_stacked[4*i - 3];
    y2_mean[i] = y_mean_stacked[4*i - 2];
    y3_mean[i] = y_mean_stacked[4*i - 1];
    y4_mean[i] = y_mean_stacked[4*i];
  }

  total_counts_mean = to_array_1d(y1_mean[time_index]);
  Nfd_mean = to_array_1d(y2_mean[time_index]);
  donor_ki_mean = to_array_1d(y3_mean[time_index]);
  host_ki_mean = to_array_1d(y4_mean[time_index]);
}

model{
  // prior distribution for model parameters
  psi ~ normal(0.5, 0.25);
  N0_Log ~ normal(14, 1);
  kappa_0 ~  normal(0.1, 0.15);
  rho ~ normal(0.005, 0.2);
  delta ~ normal(0.04, 0.2);
  Beta ~ normal(3.5, 0.8);

  sigma_chi_counts ~ normal(0, 2);
  sigma_Nfd ~ normal(0, 2);
  sigma_donor_ki ~ normal(0, 2);
  sigma_host_ki ~ normal(0, 2);

  // model fitting on to data
  log(counts) ~ normal(log(total_counts_mean), sigma_chi_counts);
  asinsqrt_array(Nfd) ~ normal(asinsqrt_array(Nfd_mean), sigma_Nfd);
  asinsqrt_array(ki_donor) ~ normal(asinsqrt_array(donor_ki_mean), sigma_donor_ki);
  asinsqrt_array(ki_host) ~ normal(asinsqrt_array(host_ki_mean), sigma_host_ki);
}

generated quantities{
  // ODE predictions for 3 different timecourses based on the age at BMT bins
  real y_chi_pred1[numPred, 4];
  real y_chi_pred2[numPred, 4];
  real y_chi_pred3[numPred, 4];

  // model predictions for the datasets for 3 different timecourses based on the age at BMT bins
  real y1_mean_pred1[numPred]; real y2_mean_pred1[numPred]; real y3_mean_pred1[numPred]; real y4_mean_pred1[numPred];
  real y1_mean_pred2[numPred]; real y2_mean_pred2[numPred]; real y3_mean_pred2[numPred]; real y4_mean_pred2[numPred];
  real y1_mean_pred3[numPred]; real y2_mean_pred3[numPred]; real y3_mean_pred3[numPred]; real y4_mean_pred3[numPred];

  // model predictions for the datasets with stdev estimated from the fits for 3 different timecourses based on the age at BMT bins
  real chicounts_pred1[numPred]; real Nfd_pred1[numPred]; real donorki_pred1[numPred]; real hostki_pred1[numPred];
  real chicounts_pred2[numPred]; real Nfd_pred2[numPred]; real donorki_pred2[numPred]; real hostki_pred2[numPred];
  real chicounts_pred3[numPred]; real Nfd_pred3[numPred]; real donorki_pred3[numPred]; real hostki_pred3[numPred];

  // log likelihood for individual data sets and combined
  vector[numObs] log_lik1; vector[numObs] log_lik2; vector[numObs] log_lik3; vector[numObs] log_lik4;

  real eps_donor = 0.99; real eps_host = 0.95;

  // initial conditions
  real init_cond[4];
  init_cond[1] = 0.0;
  init_cond[2] = 0.0;
  init_cond[3] = kappa_0 * N0;
  init_cond[4] = (1 - kappa_0) * N0;

  //ODE solution for different age bins
  y_chi_pred1 = solve_ode_chi(ts_pred1, tb_pred1, init_cond, to_array_1d(global_params));
  y_chi_pred2 = solve_ode_chi(ts_pred2, tb_pred2, init_cond, to_array_1d(global_params));
  y_chi_pred3 = solve_ode_chi(ts_pred3, tb_pred3, init_cond, to_array_1d(global_params));

  for (i in 1:numPred){
    y1_mean_pred1[i] = y_chi_pred1[i, 1] + y_chi_pred1[i, 2] + y_chi_pred1[i, 3] + y_chi_pred1[i, 4];
    y1_mean_pred2[i] = y_chi_pred2[i, 1] + y_chi_pred2[i, 2] + y_chi_pred2[i, 3] + y_chi_pred2[i, 4];
    y1_mean_pred3[i] = y_chi_pred3[i, 1] + y_chi_pred3[i, 2] + y_chi_pred3[i, 3] + y_chi_pred3[i, 4];

    y2_mean_pred1[i] = (y_chi_pred1[i, 1] + y_chi_pred1[i, 2])/(y3_mean_pred1[i] * Chi_spline(ts_pred1[i] - 45));
    y2_mean_pred2[i] = (y_chi_pred2[i, 1] + y_chi_pred2[i, 2])/(y3_mean_pred2[i] * Chi_spline(ts_pred2[i] - 67));
    y2_mean_pred3[i] = (y_chi_pred3[i, 1] + y_chi_pred3[i, 2])/(y3_mean_pred3[i] * Chi_spline(ts_pred3[i] - 89));

    y3_mean_pred1[i] = y_chi_pred1[i, 1]/(y_chi_pred1[i, 1] + y_chi_pred1[i, 2]);
    y3_mean_pred2[i] = y_chi_pred2[i, 1]/(y_chi_pred2[i, 1] + y_chi_pred2[i, 2]);
    y3_mean_pred3[i] = y_chi_pred3[i, 1]/(y_chi_pred3[i, 1] + y_chi_pred3[i, 2]);

    y4_mean_pred1[i] = y_chi_pred1[i, 3]/(y_chi_pred1[i, 3] + y_chi_pred1[i, 4]);
    y4_mean_pred2[i] = y_chi_pred2[i, 3]/(y_chi_pred2[i, 3] + y_chi_pred2[i, 4]);
    y4_mean_pred3[i] = y_chi_pred3[i, 3]/(y_chi_pred3[i, 3] + y_chi_pred3[i, 4]);

    chicounts_pred1[i] = exp(normal_rng(log(y1_mean_pred1[i]), sigma_chi_counts));
    chicounts_pred2[i] = exp(normal_rng(log(y1_mean_pred2[i]), sigma_chi_counts));
    chicounts_pred3[i] = exp(normal_rng(log(y1_mean_pred3[i]), sigma_chi_counts));

    Nfd_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y2_mean_pred1[i]), sigma_Nfd));
    Nfd_pred2[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y2_mean_pred2[i]), sigma_Nfd));
    Nfd_pred3[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y2_mean_pred3[i]), sigma_Nfd));

    donorki_pred1[i] =asinsqrt_inv(normal_rng(asinsqrt_real(y3_mean_pred1[i]), sigma_donor_ki));
    donorki_pred2[i] =asinsqrt_inv(normal_rng(asinsqrt_real(y3_mean_pred2[i]), sigma_donor_ki));
    donorki_pred3[i] =asinsqrt_inv(normal_rng(asinsqrt_real(y3_mean_pred3[i]), sigma_donor_ki));

    hostki_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y4_mean_pred1[i]), sigma_host_ki));
    hostki_pred2[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y4_mean_pred2[i]), sigma_host_ki));
    hostki_pred3[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y4_mean_pred3[i]), sigma_host_ki));
  }


  // calculating the log predictive accuracy for each point
  for (n in 1:numObs) {
    log_lik1[n] = normal_lpdf(counts[n] | log(y1_mean[n]), sigma_chi_counts);
    log_lik2[n] = normal_lpdf(Nfd[n] | asinsqrt_real(y2_mean[n]), sigma_Nfd);
    log_lik3[n] = normal_lpdf(ki_donor[n] | asinsqrt_real(y3_mean[n]), sigma_donor_ki);
    log_lik4[n] = normal_lpdf(ki_host[n] | asinsqrt_real(y4_mean[n]), sigma_host_ki);
  }
}
