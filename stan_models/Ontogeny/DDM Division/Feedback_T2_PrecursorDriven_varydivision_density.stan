functions {
  //Function that describes changes in the source chimerism
  real chi_source(real time) {
    real chi_est = 0.83;
    //real nu = 0.1257;
    real nu = 0.08;
    //real chi_0 = 0.273067;
    real t0 = 10; // t0 is the minimun days post bmt 
    real chi_s;
    
    if ((time - t0) < 0) {
      chi_s = 0; // conditioning the function to adapt to the timepoints before BMT
    } else {
      chi_s = chi_est * (1 - exp(-nu * (time - t0)));
  
    }
    
    return chi_s;
  }

    //Function that describes changes in the Ki67+ donor cells 
  real eps_donor(real time) {

    real m = -1.9e-4;
    real c = 0.97;
    real t0 = 40; // t0 is the minimun days post bmt 
    real e_d;
    
    if ((time - t0) < 0) {
      e_d = 0; // conditioning the function to adapt to the timepoints before BMT
    } else {
      e_d = m*(time-t0) + c;
  
    }

    return e_d;
  }

    //Function that describes changes in the Ki67+ host cells 
  real eps_host(real time) {

    real m = -3.8e-5;
    real c = 0.88;
    real t0 = 40; // t0 is the minimun days post bmt 
    real e_h;
    
    if ((time - t0) < 0) {
      e_h = 0; // conditioning the function to adapt to the timepoints before BMT
    } else {
      e_h = m*(time-t0) + c;
  
    }

    return e_h;
  }

  //Function that describes changes in the source T2 counts over time
  real P(real time) {
    real P_est = 14.1;
    real n = 1.8;
    real nu = 0.25;
    real t0 = 10; // t0 is the minimun days post bmt 
    real P;

    P = exp(P_est) * (1 + (time - t0)^n * exp(-nu * (time - t0)));
    return P;}
   
   //Varying the division rate with the MZ counts
   real b_m(real M, real b, real K) {
   real r_m = b * (1 - (M/K));
   //return fmax(0.0, r_m);
   return r_m;
    }
  
  //Function that contain the ODE equations to be used in the ODE solver
  vector Ode_sys_ont(real time, //time
                 vector M, //state,
                 array[] real parms //parameters
                 ) {
    //real M_stable = parms[3]; // influx rate
    real phi = parms[1]; // division rate
    real b = parms[2]; // division rate
    real d = parms[3]; // death rate
    real K = parms[4]; // Population size at which the death rate had doubled
  

    // real theta_spline = 1308516;
    // real lamda = delta - rho; //net loss rate
    //real delta = ((phi* P(time))/M_stable) + rho; 
    //real phi = (lamda * M_stable)/ theta_spline; // influx rate at time t
  
    // ODE for calulating the change in MZ counts over time
    vector[1] dMdt;
    dMdt[1] = phi * P(time) + b_m(M[1], b, K) * M[1] - d * M[1];
    return dMdt;
  }
 
  real ode_pred_ont(real solve_time, array[] real parms) {
   vector[1] M0 = [12058]'; //MZ counts at t=11
    real t0 = 11; //time at t=11
      vector[1] sol = ode_rk45(Ode_sys_ont, M0, t0, rep_array(solve_time, 1), parms)[1];
      return sol[1];
  }
  
 //Function that contain the ODE equations to be used in the ODE solver
  vector Ode_sys(real time, //time
                 vector X, //state,
                 array[] real parms //parameters
                 ) {
    // parameters
    real phi = parms[1]; //loss rate
    real b = parms[2]; //division rate
    real d = parms[3]; // death rate
    real K = parms[4]; // Population size at which the death rate had doubled
    real ageATBMT = parms[5]; //age at BMT
    
    // fixed parameters
    real P = 1308516; //mean of the precursor counts 
    real beta = 3.5;

    //real M_bar= X[1] + X[2] + X[3] + X[4]; // total counts at steady state
    //real delta = ((phi* P)/M_bar) + rho; //influx rate calculated from the loss and division rate
    //real phi = ((delta - rho) * M_bar) / P; // influx rate at time t

    // the system of ODEs
    vector[4] dXdt;
    //donor Ki67+ cells
    dXdt[1] = phi * eps_donor(time) * chi_source(time - ageATBMT) * P
              + b_m((X[1] + X[2]+ X[3] + X[4]), b, K) * (X[1] + 2 * X[2]) - ((1 / beta) + d) * X[1];
    //donor Ki67- cells
    dXdt[2] = phi * (1 - eps_donor(time)) * chi_source(time - ageATBMT) * P
              + (1 / beta) * X[1] - (b_m((X[1] + X[2]+ X[3] + X[4]), b, K) + d) * X[2];
    //host Ki67+ cells
    dXdt[3] = phi * eps_host(time) * (1 - chi_source(time - ageATBMT)) * P
              + b_m((X[1] + X[2]+ X[3] + X[4]), b, K) * (X[3] + 2 * X[4]) - ((1 / beta) + d) * X[3];
    //host Ki67- cells
    dXdt[4] = phi * (1 - eps_host(time)) * (1 - chi_source(time - ageATBMT))* P
             + (1 / beta) * X[3] - (b_m((X[1] + X[2]+ X[3] + X[4]), b, K) + d) * X[4];
    return dXdt;
  }

 //we define a solver to generate initial conditions at each age at BMT, using the earliest age at BMT as t_0
  vector initial_cond_generator(vector initial_cond, real solve_ageatBMT,
                                 array[] real parms) {
    array[1] vector[4] solution;
    array[5] real par_tb; // parameters for the intial_cond_generator
    par_tb[1 : 4] = parms[1 : 4];
    par_tb[5] = 40.0; // use the earliest age at BMT 
    solution = ode_rk45(Ode_sys, initial_cond, 40, rep_array(solve_ageatBMT, 1), par_tb); // ode_rk5 is the ode solver which give an array vector of 4
    return to_vector(solution[1]);
  }
  
  // we define a solver to generate the predicted values at each time point using initial conditions generated by the initial_cond_generator
  vector prediction_generator(real solve_time, real ageatBMT, array[] real parms){
  vector[4] initial_cond;
  vector[4] initial_tbmt; // solution of the initial_cond_generator
  vector[4] y0; // initial conditions for the solve time evaluation
  array[1] vector[4] solution;
  vector[1] N0; 
  
  array[5] real par_tb; // parameters for the intial_cond_generator
  par_tb[1 : 4] = parms[1 : 4];
  par_tb[5] = ageatBMT;
  
  N0[1] = ode_pred_ont(40, parms[1:4]); // initial conditions at t=10
  // initial conditions at time t=11
  initial_cond[1] = 0; //MZ counts at t=11
  initial_cond[2] = 0; // donor Ki67- cells
  initial_cond[3] = N0[1] * parms[5]; // host Ki67+ cells
  initial_cond[4] = N0[1] * (1 - parms[5]); // host Ki67- cells
  
  initial_tbmt = initial_cond_generator(initial_cond, ageatBMT, par_tb);
  // For every time point(each mouse), we set the first two initial conditions (for donor compartment) to be zero because age at bmt(t0), donor counts are zero.
  y0[1] = 0;
  y0[2] = 0;
  y0[3] = initial_tbmt[1] + initial_tbmt[3];
  y0[4] = initial_tbmt[2] + initial_tbmt[4];
  
  solution = ode_rk45(Ode_sys, y0, ageatBMT, rep_array(solve_time, 1), par_tb);
  return to_vector(solution[1]);
  }
  
  // define a function to solve the ODEs for each solve time using corresponding age at BMT
  vector ode_pred(array[] real pred_time, real ageatBMT,array[] real parms) {
  vector[4] initial_cond;
  int num_time = size(pred_time);
  vector[1] N0; 
  vector[4 * num_time] solution;
  
  vector[4] initial_tbmt; // solution of the initial_cond_generator
  vector[4] y0; // initial conditions for the solve time evaluation
  array[5] real par_tb; // parameters for the intial_cond_generator
  
  par_tb[1 : 4] = parms[1 : 4];
  par_tb[5] = ageatBMT;

  
  N0[1] = ode_pred_ont(40, parms[1 : 4]); // initial conditions at t=10
  // initial conditions at time t=11
  initial_cond[1] = 0; //MZ counts at t=11
  initial_cond[2] = 0; // donor Ki67- cells
  initial_cond[3] = N0[1] * parms[5]; // host Ki67+ cells
  initial_cond[4] = N0[1] * (1 - parms[5]); // host Ki67- cells
  
  initial_tbmt = initial_cond_generator(initial_cond, ageatBMT, par_tb);
  
  // For every time point(each mouse), we set the first two initial conditions (for donor compartment) to be zero because age at bmt(t0), donor counts are zero.
  y0[1] = 0;
  y0[2] = 0;
  y0[3] = initial_tbmt[1] + initial_tbmt[3];
  y0[4] = initial_tbmt[2] + initial_tbmt[4];
  
  for (i in 1 : num_time) {
  solution[4 * i - 3 : 4 * i] = ode_rk45(Ode_sys, y0, ageatBMT, rep_array(pred_time[i], 1),
                                par_tb)[1];                                                     
  }
  
  return solution;
  }
  
  // define a function to parallelize the ode_solver
  vector map_rect_helper_func(vector eta, vector theta, array[] real x_r,
   array[] int x_i) {
  array[5] real parms = to_array_1d(eta[1 : 5]);
  real solve_time = x_r[1];
  real ageatBMT = x_r[2];
  
  vector[4] solution = prediction_generator(solve_time, ageatBMT, parms);
  return solution;
  }
  
  array[] vector ode_solver_parallel( data array[,] real solve_time_ageatBMT,
          array[] real parms) {
  int num_time = size(solve_time_ageatBMT[ : , 1]);
  vector[5] eta = to_vector(parms);
  array[num_time] vector[0] theta;
  
  array[num_time, 0] int x_i;
  
  vector[4 * num_time] solution_concat = map_rect(map_rect_helper_func,
                         eta, theta, solve_time_ageatBMT, x_i);
  
  array[num_time] vector[4] solution;
  
  for (i in 1 : num_time) {
  solution[i,  : ] = solution_concat[(4 * (i - 1) + 1) : (4 * i)];
  }
  
  return solution;
  }
  
}
data {
  int<lower=1> numobs; // number of observations in the dataset
  int<lower=1> numobsNfd;
  int<lower=1> numsolve; // number of unique time points for ODE solver
  int<lower=1> numsolveNfd; // number of unique time points for ODE solver
  int<lower=1> numpred;
  int<lower=1> numpred1; // number of predicted time points 
  int<lower=1> numpred2; // number of predicted time points
  array[numobs] real data_time; // observed Nfd dataset
  array[numobs] int time_index_map; // time map for dataset
  array[numobsNfd] int time_index_map_Nfd; // time map for dataset
  array[numsolve] real solve_time; // solve time for ode solver for Ki67 dataset
  array[numsolveNfd] real solve_time_Nfd; // solve time for ode solver for Nfd dataset
  array[numsolveNfd] real solve_ageatBMT;
  array[numobsNfd] real dayspbmt; // observed data set 
  //array[numsolve] real ageatbmt;
  array[numobs] real totalMZ;
  array[numobsNfd] real Nfd;
  array[numobsNfd] real ki_donor;
  array[numobsNfd] real ki_host;
  array[numpred] real time_pred_solve; // predicted time points
  array[numpred1] real time_pred1_solve; // predicted time points
  array[numpred2] real time_pred2_solve;
  real time_pred1_bmt;
  real time_pred2_bmt;
}

transformed data {
  array[numobs] real y1;
  array[numobsNfd] real y2;
  array[numobsNfd] real y3;
  array[numobsNfd] real y4;

  array[numsolveNfd, 2] real solve_time_ageatBMT;
  solve_time_ageatBMT[ : , 1] = solve_time_Nfd;
  solve_time_ageatBMT[ : , 2] = solve_ageatBMT;

  y1 = log(totalMZ); // transforming cell counts of donor compartments to feed in to ODEs
  y2 = Nfd; // transfored donor fractions normalised to source chimerism to feed in to ODEs
  y3 = ki_donor; // transforming counts of ki67 positive cells in the donor compartments to feed in to ODEs
  y4 = ki_host; // transforming counts of ki67 positive cells in the host compartments to feed in to ODEs
}
parameters {
  //real M0_Log;
  //real<lower=0> M_stable; // steady state MZ counts
  real<lower=0> phi; // influx rate
  real<lower=0> d; 
  real<lower=0> b;
  real<lower=0, upper=15.3> K_Log; // Population size at which the death rate had doubled
  real<lower=0, upper=1> kappa_0;
  vector<lower=0>[4] sigma; // measurement error
 
 }

transformed parameters {
  // solution of the system of ODEs for the predictor values
  array[numsolveNfd] vector[4] k_hat; // declaring the array for ODE solution
  array[numobs] real totalMz_mean; // predictions for dataset1 based on ODE solution
  array[numobsNfd] real Nfd_mean; // predictions for dataset2 based on ODE solution
  array[numobsNfd] real Kid_mean; // predictions for dataset3 based on ODE solution
  array[numobsNfd] real Kih_mean; // predictions for dataset4 based on ODE solution

  array[numsolve] real MZ_counts_mean;
  array[numsolveNfd] real donor_fractions_mean;
  array[numsolveNfd] real host_ki_mean;
  array[numsolveNfd] real donor_ki_mean;

  // initial conditions and parameters
  real K = exp(K_Log);
  array[4] real parms_ont = {phi, b, d, K};
  array[5] real parms = {phi, b, d, K, 0.16}; // parameters to be passed to the ODE solver
  

// solution of the system of ODEs for the predictor values
  for (i in 1 : numsolve) {
   MZ_counts_mean[i] = ode_pred_ont(solve_time[i], parms_ont);

  }

  k_hat = ode_solver_parallel(solve_time_ageatBMT, parms);

  for (i in 1 : numsolveNfd) {

    // donor fractions normalised with chimerism in the source
    donor_fractions_mean[i] = (k_hat[i, 1] + k_hat[i, 2]) / ( k_hat[i, 1] + k_hat[i, 2] + k_hat[i, 3] + k_hat[i, 4]) ;

    // fractions of ki67 positive cells in the donor compartment
    donor_ki_mean[i] = k_hat[i, 1] / (k_hat[i, 1] + k_hat[i, 2]);

    // fractions of ki67 positive cells in the host compartment
    host_ki_mean[i] = k_hat[i, 3] / (k_hat[i, 3] + k_hat[i, 4]);

    
  }

  for (i in 1 : numobs) {
    // total counts
    totalMz_mean[i] = MZ_counts_mean[time_index_map[i]];
  }

  for (i in 1 : numobsNfd) {
   // donor fractions normalised with chimerism in the source
    Nfd_mean[i] = donor_fractions_mean[time_index_map_Nfd[i]] / chi_source(dayspbmt[i]);
    //Nfd_mean[i] = donor_fractions_mean[time_index_map[i]] / 0.79;
 
    // fractions of ki67 positive cells in the donor compartment
    Kid_mean[i] = donor_ki_mean[time_index_map_Nfd[i]];

    // fractions of ki67 positive cells in the host compartment
    Kih_mean[i] = host_ki_mean[time_index_map_Nfd[i]];

  }
}
model {
 // prior distributions for the parameters
  K_Log ~ normal(14, 2);
  b ~ normal(0.1, 0.25);
  d ~ normal(0.035, 0.0034);
  phi ~ normal(0.019, 0.0026);

 // prior distributions for the measurement error
 sigma[1] ~ normal(0.6, 0.5);
 sigma[2] ~ normal(3, 1);
 sigma[3] ~ normal(0.5, 0.5);
 sigma[4] ~ normal(0.9, 1);

  // likelihood of the observed data given the model prediction 
  y1 ~ normal(log(totalMz_mean), sigma[1]);
  y2 ~ normal(Nfd_mean, sigma[2]);
  y3 ~ normal(Kid_mean, sigma[3]);
  y4 ~ normal(Kih_mean, sigma[4]);
}
generated quantities {

array[numpred] real totalMz_pred;
  for (i in 1:numpred) {
    totalMz_pred[i] = ode_pred_ont(time_pred_solve[i], parms_ont);
  }
  
  vector[4 * numpred1] k_hat_pred1 = ode_pred(time_pred1_solve, time_pred1_bmt, parms);
  vector[4 * numpred2] k_hat_pred2 = ode_pred(time_pred2_solve, time_pred2_bmt, parms);

  // Predictions for the predicted time points
  array[numpred1] real MZ_counts_pred1;
  array[numpred1] real Nfd_pred1;
  array[numpred1] real host_ki_pred1;
  array[numpred1] real donor_ki_pred1;

  array[numpred2] real MZ_counts_pred2;
  array[numpred2] real Nfd_pred2;
  array[numpred2] real host_ki_pred2;
  array[numpred2] real donor_ki_pred2;

  // Errors for the predicted time points
  // array[numpred1] real MZ_counts_pred1_err;
  // array[numpred1] real Nfd_pred1_err;
  // array[numpred1] real host_ki_pred1_err;
  // array[numpred1] real donor_ki_pred1_err;

  // array[numpred2] real MZ_counts_pred2_err;
  // array[numpred2] real Nfd_pred2_err;
  // array[numpred2] real host_ki_pred2_err;
  // array[numpred2] real donor_ki_pred2_err;

  
  for (i in 1 : numpred1) {
    MZ_counts_pred1[i] = k_hat_pred1[4 * i - 3] + k_hat_pred1[4 * i - 2] + k_hat_pred1[4 * i - 1] + k_hat_pred1[4 * i];
    //MZ_counts_pred1_err[i] = exp(normal_rng(MZ_counts_pred1[i], sigma[1]));

    Nfd_pred1[i] = (k_hat_pred1[4 * i - 3] + k_hat_pred1[4 * i - 2]) / (MZ_counts_pred1[i] * chi_source(time_pred1_solve[i] - time_pred1_bmt));
    //Nfd_pred1_err[i] = normal_rng(Nfd_pred1[i], sigma[2]);
    //Nfd_pred1[i] = (k_hat_pred1[4 * i - 3] + k_hat_pred1[4 * i - 2]) / (MZ_counts_pred1[i] * 0.79);

    donor_ki_pred1[i] = k_hat_pred1[4 * i - 3] / (k_hat_pred1[4 * i - 3] + k_hat_pred1[4 * i - 2]);
    //donor_ki_pred1_err[i] = normal_rng(donor_ki_pred1[i], sigma[3]);

    host_ki_pred1[i] = k_hat_pred1[4 * i - 1] / (k_hat_pred1[4 * i - 1] + k_hat_pred1[4 * i]);
   // host_ki_pred1_err[i] = normal_rng(host_ki_pred1[i], sigma[4]);
  }

  for (i in 1 : numpred2) {
    MZ_counts_pred2[i] = k_hat_pred2[4 * i - 3] + k_hat_pred2[4 * i - 2] + k_hat_pred2[4 * i - 1] + k_hat_pred2[4 * i];
    //MZ_counts_pred2_err[i] = exp(normal_rng(MZ_counts_pred2[i], sigma[1]));

    Nfd_pred2[i] = (k_hat_pred2[4 * i - 3] + k_hat_pred2[4 * i - 2]) / (MZ_counts_pred2[i] * chi_source(time_pred2_solve[i] - time_pred2_bmt));
    //Nfd_pred2_err[i] = normal_rng(Nfd_pred2[i], sigma[2]);

    donor_ki_pred2[i] = k_hat_pred2[4 * i - 3] / (k_hat_pred2[4 * i - 3] + k_hat_pred2[4 * i - 2]);
    //donor_ki_pred2_err[i] = normal_rng(donor_ki_pred2[i], sigma[3]);
    
    host_ki_pred2[i] = k_hat_pred2[4 * i - 1] / (k_hat_pred2[4 * i - 1] + k_hat_pred2[4 * i]);
    //host_ki_pred2_err[i] = normal_rng(host_ki_pred2[i], sigma[4]);
  }

  // log likelihood of the observed data given the model prediction
  vector[numobs] log_lik_Mz;
  vector[numobsNfd] log_lik_Nfd;
  vector[numobsNfd] log_lik_Kid;
  vector[numobsNfd] log_lik_Kih;
 
  // calculating the log predictive accuracy for each point
  for (i in 1 : numobs) {
    log_lik_Mz[i] = normal_lpdf(y1[i] | log(totalMz_mean[i]), sigma[1]);
}
  for (i in 1 : numobsNfd) {
    log_lik_Nfd[i] = normal_lpdf(y2[i] | Nfd_mean[i], sigma[2]);
    log_lik_Kid[i] = normal_lpdf(y3[i] | Kid_mean[i], sigma[3]);
    log_lik_Kih[i] = normal_lpdf(y4[i] | Kih_mean[i], sigma[4]);
  }

// Calculate residuals for the observed data
  array[numobs] real residuals_Mz;
  array[numobsNfd] real residuals_Nfd;
  array[numobsNfd] real residuals_Kid;
  array[numobsNfd] real residuals_Kih;

  for (i in 1 : numobs) {
    residuals_Mz[i] = y1[i] - log(totalMz_mean[i]);
}
  for (i in 1 : numobsNfd) {
    residuals_Nfd[i] = y2[i] - Nfd_mean[i];
    residuals_Kid[i] = y3[i] - Kid_mean[i];
    residuals_Kih[i] = y4[i] - Kih_mean[i];
  }
}