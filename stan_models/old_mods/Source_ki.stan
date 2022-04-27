functions{
  real ki_spline( real t,
      real ki_0,
      real ki_slope) {

        real ki_frac;

        ki_frac = ki_0 * exp(-ki_slope * t);
        return ki_frac;
    }
}

data{
int<lower = 1> numObs;
int<lower = 0>  Time[numObs];
int<lower  = 1> numPred;
int<lower  = 0> ts_pred[numPred];
real<lower = 0> ki_host[numObs];
real<lower = 0> ki_donor[numObs];
}

transformed data{
real y1[numObs];
real y2[numObs];

for (i in 1:numObs){
  y1[i] = ki_host[i];
  y2[i] = ki_donor[i];
 }
}

parameters{
  real<lower = 0, upper = 1> ki_0_host;
  real<lower = 0, upper = 1> ki_0_donor;
  real ki_slope_Log;
  //real ki_slope_DLog;
  real<lower = 0> sigma_ki_host;
  real<lower = 0> sigma_ki_donor;
}

transformed parameters{
 real y1mean[numObs];
 real y2mean[numObs];
 real ki_slope_host = exp(ki_slope_Log);
 real ki_slope_donor = exp(ki_slope_Log);

 for(i in 1:numObs){
   y1mean[i] = ki_spline(Time[i], ki_0_host, ki_slope_host);
   y2mean[i] = ki_spline(Time[i], ki_0_donor, ki_slope_donor);
 }
}

model{
  ki_0_host ~ normal(0.8, 0.1);
  ki_0_donor ~ normal(0.99, 0.1);
  ki_slope_Log ~ normal(-4, 2);
  //ki_slope_DLog ~ normal(-4, 2);

  sigma_ki_host ~ normal(0, 2);
  sigma_ki_donor ~ normal(0, 2);

  y1 ~ normal(y1mean, sigma_ki_host);
  y2 ~ normal(y2mean, sigma_ki_donor);
}

generated quantities{
  real y_kihost_pred[numPred];
  real ki_host_pred[numPred];
  real y_kidonor_pred[numPred];
  real ki_donor_pred[numPred];


  for(i in 1:numPred){
  y_kihost_pred[i] = ki_spline(ts_pred[i], ki_0_host, ki_slope_host);
  y_kidonor_pred[i] = ki_spline(ts_pred[i], ki_0_donor, ki_slope_donor);
  }

for(i in 1:numPred){
  ki_donor_pred[i] = normal_rng(y_kidonor_pred[i], sigma_ki_donor);
  ki_host_pred[i] = normal_rng(y_kihost_pred[i], sigma_ki_host);
    }
}
