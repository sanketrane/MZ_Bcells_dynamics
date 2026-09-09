functions {
   /* ... function declarations and definitions ... */
   vector emp_func(vector time, real theta_0, real nu, real t0) {
    int nsize= size(time);
    vector[nsize] y_pred;
    for (i in 1:nsize){
      y_pred[i] = exp(theta_0) / (1+exp(-nu * (time[i]-t0)));
    }
    return y_pred;
   }
}

data {
  int<lower=0> numobs;   // number of observations 
  int<lower=0> num_pred; // number of observations for prediction
  vector[numobs] numsolve;    // observed time point for fitting
  vector[numobs] numFOcounts; // T1 counts for observed data  
  vector[num_pred] time_pred; // time points for prediction 
}

parameters {
  real theta_0;  // intercept
  real<lower=0> nu;  // coefficients
  real<lower=0> t0;
  real<lower=0, upper=1> sig; // standard deviation
}
transformed parameters {
  // real theta_0_Log = exp(theta_0);
  vector[numobs] y1;
  y1 = log(numFOcounts);
  vector[numobs] y_pred = emp_func(numsolve, theta_0, nu, t0);
}

model {
  theta_0 ~ lognormal(14, 3);
  nu ~ normal(0.5, 1); 
  t0 ~ normal(30, 2); 
  sig ~ normal(0.5, 0.5);//normal prior on interval 0, 1
  y1 ~ normal(log(y_pred), sig); 
}

generated quantities {
  vector [num_pred] FO_counts_pred;
  FO_counts_pred = emp_func(time_pred, theta_0, nu,t0);
  
}

