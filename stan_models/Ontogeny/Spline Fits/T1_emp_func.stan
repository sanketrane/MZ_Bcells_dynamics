functions {
   /* ... function declarations and definitions ... */
   vector emp_func(vector time, real theta_0, real nu, real n) {
    int nsize= size(time);
    vector[nsize] y_pred;
    for (i in 1:nsize){
      y_pred[i] = exp(theta_0) * (1+ (time[i]-10)^n * exp(-nu * (time[i]-10)));
    }
    return y_pred;
   }
}

data {
  int<lower=0> numobs;   // number of observations 
  int<lower=0> num_pred; // number of observations for prediction
  vector[numobs] numsolve;    // observed time point for fitting
  vector[numobs] numT1counts; // T1 counts for observed data  
  vector[num_pred] time_pred; // time points for prediction 
}

parameters {
  real <lower=13, upper=15> theta_0;  // intercept
  real nu;  // coefficients
  real n;
  real<lower=0, upper=1> sig; // standard deviation
}
transformed parameters {
  // real theta_0_Log = exp(theta_0);
  vector[numobs] y1;
  y1 = log(numT1counts);
  vector[numobs] y_pred = emp_func(numsolve, theta_0, nu, n);
}

model {
  theta_0 ~ lognormal(14, 1);
  nu ~ normal(0.5, 1); 
  n ~ normal(2, 0.5); 
  sig ~ normal(0.5, 0.5);//normal prior on interval 0, 1
  y1 ~ normal(log(y_pred), sig); 
}

generated quantities {
  vector [num_pred] T1_counts_pred;
  T1_counts_pred = emp_func(time_pred, theta_0, nu,n);
  
}

