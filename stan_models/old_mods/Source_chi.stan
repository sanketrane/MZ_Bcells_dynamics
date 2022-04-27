functions{
  real Chi_spline( real t,
      real chiEst,
      real qEst) {

        real chi;

        chi = chiEst * (1 - exp(-qEst * t));
        return chi;
    }
}

data{
int<lower = 1> numObs;
int<lower = 0>  Time[numObs];
int<lower  = 1> numPred;
int<lower  = 0> ts_pred[numPred];
real<lower = 0> chi_source[numObs];
}

transformed data{
real y[numObs];

for (i in 1:numObs){
y[i] = chi_source[i];
}
}

parameters{
  real<lower = 0, upper = 1> chiEst;
  real qEstlog;
  real<lower = 0> sigma_chi;
}

transformed parameters{
 real ymean[numObs];
 real qEst;

 qEst = exp(qEstlog);

 for(i in 1:numObs){
 ymean[i] = Chi_spline(Time[i], chiEst, qEst);
 }
}

model{
  chiEst ~ uniform(0, 1);
  qEstlog ~ normal(-2, 2);
  sigma_chi ~ normal(0, 1);

  y ~ normal(ymean, sigma_chi);
}

generated quantities{
  real y_chipred[numPred];
  real chipred[numPred];

  for(i in 1:numPred){
  y_chipred[i] = Chi_spline(ts_pred[i], chiEst, qEst);
  }

for(i in 1:numPred){
    chipred[i] = normal_rng(y_chipred[i], sigma_chi);
    }
}
