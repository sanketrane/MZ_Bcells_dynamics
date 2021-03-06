functions{
  real theta_spline(real t,
    real nu,               // nu is the rate of decline of source influx assuming that its proportional to the counts of source compartment
    real theta0) {     // theta0 gives the initial counts of the source compartment

      real theta;
      int ta = 40;         // the earliest age at BMT

      theta = theta0 * exp(-nu * (t-ta));
      return theta;
   }

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

  real[] khm(real t,  real[] k, real[] parms, real[] rdata, int[] idata) {
    real psi        = parms[1];
    real f_fast     = parms[2];
    real rhoFast    = parms[3];
    real rhoSlow    = parms[4];
    real Beta       = parms[5];
    real lambdaFast = parms[6];
    real lambdaSlow = parms[7];

    real theta0    = rdata[1];
    real nu        = rdata[2];
    real chiEst    = rdata[3];
    real qEst      = rdata[4];
    real eps_donor = rdata[5];
    real eps_host  = rdata[6];

    real tb = idata[1];

    real dkdt[8];

    dkdt[1] = (psi * f_fast * theta_spline(t, nu, theta0) * Chi_spline(t-tb, chiEst, qEst) * eps_donor) + rhoFast * (2 * k[2] + k[1]) - ((1/Beta) + (lambdaFast+rhoFast)) * k[1];
    dkdt[2] = (psi * f_fast * theta_spline(t, nu, theta0) * Chi_spline(t-tb, chiEst, qEst) * (1-eps_donor)) + (1/Beta) * k[1] - (rhoFast + (lambdaFast+rhoFast)) * k[2];
    dkdt[3] = (psi * (1-f_fast) * theta_spline(t, nu, theta0) * Chi_spline(t-tb, chiEst, qEst) * eps_donor) + rhoSlow * (2 * k[4] + k[3]) - ((1/Beta) + (lambdaSlow+rhoSlow)) * k[3];
    dkdt[4] = (psi * (1-f_fast) * theta_spline(t, nu, theta0) * Chi_spline(t-tb, chiEst, qEst) * (1-eps_donor)) + (1/Beta) * k[3] - (rhoSlow + (lambdaSlow+rhoSlow)) * k[4];

    dkdt[5] = (psi * f_fast * theta_spline(t, nu, theta0) * (1-Chi_spline(t-tb, chiEst, qEst)) * eps_host) + rhoFast * (2 * k[6] + k[5]) - ((1/Beta) + (lambdaFast+rhoFast)) * k[5];
    dkdt[6] = (psi * f_fast * theta_spline(t, nu, theta0) * (1-Chi_spline(t-tb, chiEst, qEst)) * (1-eps_host)) + (1/Beta) * k[5] - (rhoFast + (lambdaFast+rhoFast)) * k[6];
    dkdt[7] = (psi * (1-f_fast) * theta_spline(t, nu, theta0) * (1-Chi_spline(t-tb, chiEst, qEst)) * eps_host) + rhoSlow * (2 * k[8] + k[7]) - ((1/Beta) + (lambdaSlow+rhoSlow)) * k[7];
    dkdt[8] = (psi * (1-f_fast) * theta_spline(t, nu, theta0) * (1-Chi_spline(t-tb, chiEst, qEst)) * (1-eps_host)) + (1/Beta) * k[7] - (rhoSlow + (lambdaSlow+rhoSlow)) * k[8];


    return dkdt;
  }

  // to calculate total counts (N) at each 'tb'
  // this function solves the ODE for each value of 'tb' using 'ta' as the starting point (t0) and converts the solution into a 1d array
  real[] foreach_init(real tb, real ta, real[] init_cond, real[] parms, real[] rdata, int x_i){
    //first argument to integrate_ode_rk45 must be the name of the ODE function with signature (real, real[], real[], real[], int[])
    return to_array_1d(integrate_ode_rk45(khm, init_cond, ta, rep_array(tb, 1), parms, rdata, rep_array(x_i, 1)));
  }

  // ode solution at each 'tb' -- this serves as initial condition when solving at 't'
  real[,] solve_init(real[] tb_time, real[] init_cond, real[] parms, real[] rdata, int x_i, int num_tb){
   real y_init[num_tb, 8];

   y_init[1] = init_cond;
   for (i in 2:num_tb){
     y_init[i] = foreach_init(tb_time[i], tb_time[1], init_cond, parms, rdata, x_i);
   }

   return y_init;
  }

  // this function solves the ODE for all obsreved timepoints using 'tb' as the strating point (t0) and converts the solution into a 1d array
  real[] foreach_ode(real ts, real t0, real[] init_cond, real[] parms, real[] rdata, int x_i) {
    //first argument to integrate_ode_rk45 must be the name of the ODE function with signature (real, real[], real[], real[], int[])
    return to_array_1d(integrate_ode_rk45(khm, init_cond, t0, rep_array(ts, 1), parms, rdata, rep_array(x_i, 1)));
   }

 // ode solution at each 't' i.e. obsreved timepoints
  real[,] solve_ode(real[] solve_time, real[] init_cond, real[] parms, real[] rdata, int[] tb, int num_index, real[] tb_time, int[] tb_index, int num_tb){
   real y_hat[num_index, 8];       // the array for ODE solution
   real y0[num_tb, 8];             // the array for ode solution of initial conditions
   real init_tb[8];                // the array that defines initial conditions

   //initial conidtion at each tb
   y0 = solve_init(tb_time, init_cond, parms, rdata, 40, num_tb);

   init_tb[1] = init_cond[1];                                          //at tbmt - donor ki67Hi subset size is zero
   init_tb[2] = init_cond[2];                                          //at tbmt - donor ki67Hi subset size is zero
   init_tb[3] = init_cond[3];                                          //at tbmt - donor ki67Hi subset size is zero
   init_tb[4] = init_cond[4];                                          //at tbmt - donor ki67Hi subset size is zero

   for (i in 1:num_index){
     init_tb[5] = y0[tb_index[i], 1] + y0[tb_index[i], 5];              //at tbmt - all ki67Hi cells would be host
     init_tb[6] = y0[tb_index[i], 2] + y0[tb_index[i], 6];              //at tbmt - all ki67Lo cells would be host
     init_tb[7] = y0[tb_index[i], 3] + y0[tb_index[i], 7];              //at tbmt - all ki67Hi cells would be host
     init_tb[8] = y0[tb_index[i], 4] + y0[tb_index[i], 8];              //at tbmt - all ki67Lo cells would be host
     y_hat[i] = foreach_ode(solve_time[i], tb[i], init_tb, parms, rdata, tb[i]);
   }

   return y_hat;
  }

  real[,] solve_ode_pred(real[] solve_time, real[] init_cond, real[] parms, real[] rdata, int[] tb, int num_index, real[] tb_time){
   real y_hat[num_index, 8];
   real y0[2, 8];
   real init_tb[8];

   y0 = solve_init(tb_time, init_cond, parms, rdata, 40, 2);

   init_tb[1] = init_cond[1];                                          //at tbmt - donor ki67Hi subset size is zero
   init_tb[2] = init_cond[2];                                          //at tbmt - donor ki67Hi subset size is zero
   init_tb[3] = init_cond[3];                                          //at tbmt - donor ki67Hi subset size is zero
   init_tb[4] = init_cond[4];                                          //at tbmt - donor ki67Hi subset size is zero
   init_tb[5] = y0[2, 1] + y0[2, 5];                                   //at tbmt - all ki67Hi cells would be host
   init_tb[6] = y0[2, 2] + y0[2, 6];                                   //at tbmt - all ki67Lo cells would be host
   init_tb[7] = y0[2, 3] + y0[2, 7];                                   //at tbmt - all ki67Hi cells would be host
   init_tb[8] = y0[2, 4] + y0[2, 8];                                   //at tbmt - all ki67Lo cells would be host

   y_hat[1] = init_tb;
   for (i in 2:num_index){
     y_hat[i] = foreach_ode(solve_time[i], tb[i], init_tb, parms, rdata, tb[i]);
   }

   return y_hat;
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
  y2 =  asinsqrt_array(Nfd);                  // transfored donor fractions normalised to source chimerism to feed in to ODEs
  y3 =  asinsqrt_array(ki_donor);             // transforming counts of ki67 positive cells in the  donor compartments to feed in to ODEs
  y4 =  asinsqrt_array(ki_host);              // transforming counts of ki67 positive cells in the  host compartments to feed in to ODEs
}

parameters{
  real y0_Log;
  real<lower = 0, upper = 1> kappaF_0;
  real<lower = 0, upper = 1> kappaS_0;
  real<lower = 0, upper = 1> psi;
  real<lower = 0, upper = 1> f_fast;
  real rhoFast_Log;
  real rhoSlow_Log;
  real lambdaFast_Log;
  real<upper = lambdaFast_Log>  lambdaSlow_Log;
  real<lower = 0, upper = 1>  alpha;
  real<lower = 0> Beta;

  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real<lower = 0> sigma3;
  real<lower = 0> sigma4;
  }

  transformed parameters{
    real k_hat[num_index, 8];
    real y1_mean[numObs];
    real y2_mean[numObs];
    real y3_mean[numObs];
    real y4_mean[numObs];

    real parms[7];
    real init_cond[8];

    real y0         = exp(y0_Log);                   // transformed parameters for better/faster sampling
    real rhoFast    = exp(rhoFast_Log);        // transformed parameters for better/faster sampling
    real rhoSlow    = exp(rhoSlow_Log);        // transformed parameters for better/faster sampling
    real lambdaFast = exp(lambdaFast_Log);        // transformed parameters for better/faster sampling
    real lambdaSlow = exp(lambdaSlow_Log);        // transformed parameters for better/faster sampling


    init_cond[1] = Nd_0;
    init_cond[2] = Nd_0;
    init_cond[3] = Nd_0;
    init_cond[4] = Nd_0;

    init_cond[5] = y0 * alpha * kappaF_0;
    init_cond[6] = y0 * alpha * (1- kappaF_0);
    init_cond[7] = y0 * (1-alpha) * kappaS_0;
    init_cond[8] = y0 * (1-alpha) * (1- kappaS_0);

    parms[1] = psi;
    parms[2] = f_fast;
    parms[3] = rhoFast;
    parms[4] = rhoSlow;
    parms[5] = Beta;
    parms[6] = lambdaFast;
    parms[7] = lambdaSlow;

    k_hat = solve_ode(solve_time, init_cond, parms, rdata, tb, num_index, tb_time, tb_index, num_tb);

    for (i in 1:numObs){
      // total counts
      y1_mean[i] = k_hat[time_index[i], 1] + k_hat[time_index[i], 2] + k_hat[time_index[i], 3] + k_hat[time_index[i], 4] +
       k_hat[time_index[i], 5] + k_hat[time_index[i], 6] + k_hat[time_index[i], 7] + k_hat[time_index[i], 8];

      // donor fractions normalised with chimerism in the source
      y2_mean[i] = (k_hat[time_index[i], 1] + k_hat[time_index[i], 2] + k_hat[time_index[i], 3] + k_hat[time_index[i], 4])/
      (y1_mean[i] * Chi_spline(dpBMT[i], chiEst, qEst));

      // fractions of ki67 positive cells in the donor compartment
      y3_mean[i] = (k_hat[time_index[i], 1] + k_hat[time_index[i], 3])/
      (k_hat[time_index[i], 1] + k_hat[time_index[i], 2] + k_hat[time_index[i], 3] + k_hat[time_index[i], 4]);

      // fractions of ki67 positive cells in the host compartment
      y4_mean[i] = (k_hat[time_index[i], 5] + k_hat[time_index[i], 7])/
      (k_hat[time_index[i], 5] + k_hat[time_index[i], 6] + k_hat[time_index[i], 7] + k_hat[time_index[i], 8]);
    }
}

model{
  // prior distribution for model parameters
  psi ~ normal(0.5, 0.25);
  f_fast ~ uniform(0, 1);
  rhoSlow_Log ~ normal(-4, 1.2);
  rhoFast_Log ~ normal(-4, 1.2);
  lambdaFast_Log ~ normal(-3, 1.2);
  lambdaSlow_Log ~ normal(-4, 1.2);
  alpha ~ uniform(0, 1);
  Beta ~ normal(3.5, 0.8);

  y0_Log ~ normal(14, 1);
  kappaF_0 ~ normal(0.1, 0.2);
  kappaS_0 ~ normal(0.1, 0.2);

  sigma1 ~ normal(0.6, 0.2);
  sigma2 ~ normal(3, 0.5);
  sigma3 ~ normal(0.3, 0.2);
  sigma4 ~ normal(0.9, 2);

  y1 ~ normal(log(y1_mean), sigma1);
  y2 ~ normal(asinsqrt_array(y2_mean), sigma2);
  y3 ~ normal(asinsqrt_array(y3_mean), sigma3);
  y4 ~ normal(asinsqrt_array(y4_mean), sigma4);
}

generated quantities{
  // ODE predictions for 3 different timecourses based on the age at BMT bins
  real k_hat_pred_age1[numPred1, 8];
  real k_hat_pred_age2[numPred2, 8];
  real k_hat_pred_age3[numPred3, 8];
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
  vector[numObs] log_lik1;
  vector[numObs] log_lik2;
  vector[numObs] log_lik3;
  vector[numObs] log_lik4;

  // parameters of interest
  real lambdaFast_inv = 1/lambdaFast;
  real lambdaSlow_inv = 1/lambdaSlow;
  real rhoFast_inv = 1/rhoFast;
  real rhoSlow_inv = 1/rhoFast;
  real deltaFast = lambdaFast + rhoFast;
  real deltaFast_inv = 1/deltaFast;
  real deltaSlow = lambdaSlow + rhoFast;
  real deltaSlow_inv = 1/deltaSlow;

  //ODE solution for different age bins
  k_hat_pred_age1 = solve_ode_pred(ts_pred1, init_cond, parms, rdata, tb_pred1, numPred1, tb_time_pred1);
  k_hat_pred_age2 = solve_ode_pred(ts_pred2, init_cond, parms, rdata, tb_pred2, numPred2, tb_time_pred2);
  k_hat_pred_age3 = solve_ode_pred(ts_pred3, init_cond, parms, rdata, tb_pred3, numPred3, tb_time_pred3);

  // Total cell counts (donor + host)
  for (i in 1:numPred1){
    y1_mean_pred_age1[i] = k_hat_pred_age1[i, 1] + k_hat_pred_age1[i, 2] + k_hat_pred_age1[i, 3] + k_hat_pred_age1[i, 4] +
     k_hat_pred_age1[i, 5] + k_hat_pred_age1[i, 6] + k_hat_pred_age1[i, 7] + k_hat_pred_age1[i, 8];
    countspred_age1[i] = exp(normal_rng(log(y1_mean_pred_age1[i]), sigma1));
  }

  for (i in 1: numPred2){
    y1_mean_pred_age2[i] = k_hat_pred_age2[i, 1] + k_hat_pred_age2[i, 2] + k_hat_pred_age2[i, 3] + k_hat_pred_age2[i, 4] +
     k_hat_pred_age2[i, 5] + k_hat_pred_age2[i, 6] + k_hat_pred_age2[i, 7] + k_hat_pred_age2[i, 8];
    countspred_age2[i] = exp(normal_rng(log(y1_mean_pred_age2[i]), sigma1));
  }

  for (i in 1: numPred3){
    y1_mean_pred_age3[i] = k_hat_pred_age3[i, 1] + k_hat_pred_age3[i, 2] + k_hat_pred_age3[i, 3] + k_hat_pred_age3[i, 4] +
     k_hat_pred_age3[i, 5] + k_hat_pred_age3[i, 6] + k_hat_pred_age3[i, 7] + k_hat_pred_age3[i, 8];
    countspred_age3[i] = exp(normal_rng(log(y1_mean_pred_age3[i]), sigma1));
  }

  // Initial conditions for fd and frcations of ki67host and ki67 donor
  y2_mean_pred_age1[1] = Nd_0;
  fdpred_age1[1] = Nd_0;
  y2_mean_pred_age2[1] = Nd_0;
  fdpred_age2[1] = Nd_0;
  y2_mean_pred_age3[1] = Nd_0;
  fdpred_age3[1] = Nd_0;

  y3_mean_pred1[1] = eps_donor;
  donor_kiprop_pred1[1] =  asinsqrt_inv(normal_rng( asinsqrt_real(y3_mean_pred1[1]), sigma3));
  y3_mean_pred2[1] = eps_donor;
  donor_kiprop_pred2[1] =  asinsqrt_inv(normal_rng( asinsqrt_real(y3_mean_pred2[1]), sigma3));
  y3_mean_pred3[1] = eps_donor;
  donor_kiprop_pred3[1] =  asinsqrt_inv(normal_rng( asinsqrt_real(y3_mean_pred3[1]), sigma3));

  y4_mean_pred1[1] = (y0 * alpha * kappaF_0 + y0 * (1-alpha) * kappaS_0)/(y0);
  host_kiprop_pred1[1] =  asinsqrt_inv(normal_rng( asinsqrt_real(y4_mean_pred1[1]), sigma4));
  y4_mean_pred2[1] = (y0 * alpha * kappaF_0 + y0 * (1-alpha) * kappaS_0)/(y0);
  host_kiprop_pred2[1] =  asinsqrt_inv(normal_rng( asinsqrt_real(y4_mean_pred2[1]), sigma4));
  y4_mean_pred3[1] = (y0 * alpha * kappaF_0 + y0 * (1-alpha) * kappaS_0)/(y0);
  host_kiprop_pred3[1] =  asinsqrt_inv(normal_rng( asinsqrt_real(y4_mean_pred3[1]), sigma4));

  for (i in 2:numPred1){
    y2_mean_pred_age1[i] = (k_hat_pred_age1[i, 1] + k_hat_pred_age1[i, 2] + k_hat_pred_age1[i, 3] + k_hat_pred_age1[i, 4])/(y1_mean_pred_age1[i] * Chi_spline(ts_pred1[i] - tb_pred1[i], chiEst, qEst));
    fdpred_age1[i] =  asinsqrt_inv(normal_rng( asinsqrt_real(y2_mean_pred_age1[i]), sigma2));

    y3_mean_pred1[i] = (k_hat_pred_age1[i, 1] + k_hat_pred_age1[i, 3])/(k_hat_pred_age1[i, 1] + k_hat_pred_age1[i, 2] + k_hat_pred_age1[i, 3] + k_hat_pred_age1[i, 4]);
    donor_kiprop_pred1[i] =   asinsqrt_inv(normal_rng( asinsqrt_real(y3_mean_pred1[i]), sigma3));

    y4_mean_pred1[i] = (k_hat_pred_age1[i, 5] + k_hat_pred_age1[i, 7])/(k_hat_pred_age1[i, 5] + k_hat_pred_age1[i, 6] + k_hat_pred_age1[i, 7] + k_hat_pred_age1[i, 8]);
    host_kiprop_pred1[i] =  asinsqrt_inv(normal_rng( asinsqrt_real(y4_mean_pred1[i]), sigma4));
  }

  for (i in 2:numPred2){
    y2_mean_pred_age2[i] = (k_hat_pred_age2[i, 1] + k_hat_pred_age2[i, 2] + k_hat_pred_age2[i, 3] + k_hat_pred_age2[i, 4])/(y1_mean_pred_age2[i] * Chi_spline(ts_pred2[i] - tb_pred2[i], chiEst, qEst));
    fdpred_age2[i] =  asinsqrt_inv(normal_rng( asinsqrt_real(y2_mean_pred_age2[i]), sigma2));

    y3_mean_pred2[i] = (k_hat_pred_age2[i, 1] + k_hat_pred_age2[i, 3])/(k_hat_pred_age2[i, 1] + k_hat_pred_age2[i, 2] + k_hat_pred_age2[i, 3] + k_hat_pred_age2[i, 4]);
    donor_kiprop_pred2[i] =   asinsqrt_inv(normal_rng( asinsqrt_real(y3_mean_pred2[i]), sigma3));

    y4_mean_pred2[i] = (k_hat_pred_age2[i, 5] + k_hat_pred_age2[i, 7])/(k_hat_pred_age2[i, 5] + k_hat_pred_age2[i, 6] + k_hat_pred_age2[i, 7] + k_hat_pred_age2[i, 8]);
    host_kiprop_pred2[i] =  asinsqrt_inv(normal_rng( asinsqrt_real(y4_mean_pred2[i]), sigma4));
  }

  for (i in 2:numPred3){
    y2_mean_pred_age3[i] = (k_hat_pred_age3[i, 1] + k_hat_pred_age3[i, 2] + k_hat_pred_age3[i, 3] + k_hat_pred_age3[i, 4])/(y1_mean_pred_age3[i] * Chi_spline(ts_pred3[i] - tb_pred3[i], chiEst, qEst));
    fdpred_age3[i] =  asinsqrt_inv(normal_rng( asinsqrt_real(y2_mean_pred_age3[i]), sigma2));

    y3_mean_pred3[i] = (k_hat_pred_age3[i, 1] + k_hat_pred_age3[i, 3])/(k_hat_pred_age3[i, 1] + k_hat_pred_age3[i, 2] + k_hat_pred_age3[i, 3] + k_hat_pred_age3[i, 4]);
    donor_kiprop_pred3[i] =   asinsqrt_inv(normal_rng( asinsqrt_real(y3_mean_pred3[i]), sigma3));

    y4_mean_pred3[i] = (k_hat_pred_age3[i, 5] + k_hat_pred_age3[i, 7])/(k_hat_pred_age3[i, 5] + k_hat_pred_age3[i, 6] + k_hat_pred_age3[i, 7] + k_hat_pred_age3[i, 8]);
    host_kiprop_pred3[i] =  asinsqrt_inv(normal_rng( asinsqrt_real(y4_mean_pred3[i]), sigma4));
  }


  // calculating the log predictive accuracy for each point
  for (n in 1:numObs) {
    log_lik1[n] = normal_lpdf(y1[n] | log(y1_mean[n]), sigma1);
    log_lik2[n] = normal_lpdf(y2[n] |  asinsqrt_real(y2_mean[n]), sigma2);
    log_lik3[n] = normal_lpdf(y3[n] |  asinsqrt_real(y3_mean[n]), sigma3);
    log_lik4[n] = normal_lpdf(y4[n] |  asinsqrt_real(y4_mean[n]), sigma4);
  }
}
