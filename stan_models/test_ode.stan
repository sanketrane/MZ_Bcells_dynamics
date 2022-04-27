functions{
  // function that describes the changes in source counts
  real theta_func(real time){

      real nu = 0.01; real theta0 = 9.495;

      real theta;
      int t0 = 40;         // the earliest age at BMT

      // nu is the rate of decline of source influx
      //theta0 gives the initial counts of the source compartment
      theta = exp(theta0) * exp(-nu * (time-t0));
      return theta;
   }

  // function that describes changes in source chimerism
  real Chi_T1( real time){

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
   real[] ode(real time,  real[] k, real[] parms, real[] rdata, int[] idata) {

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

     dkdt[1] = (psi * theta_func(time) * Chi_T1(time-ageAtBMT) * eps_donor) + rho * (2 * k[2] + k[1]) - ((1/Beta) + (lambda + rho)) * k[1];
     dkdt[2] = (psi * theta_func(time) * Chi_T1(time-ageAtBMT) * (1-eps_donor)) + (1/Beta) * k[1] - (rho + (lambda + rho)) * k[2];

     dkdt[3] = (psi * theta_func(time) * (1-Chi_T1(time-ageAtBMT)) * eps_host) + rho * (2 * k[4] + k[3]) - ((1/Beta) + (lambda + rho)) * k[3];
     dkdt[4] = (psi * theta_func(time) * (1-Chi_T1(time-ageAtBMT)) * (1-eps_host)) + (1/Beta) * k[3] - (rho + (lambda + rho)) * k[4];

     return dkdt;
  }

  real[,] solve_shm(real[] solve_time, real[] init_cond, real[] parms) {
    // solves the ode for each timepoint from t0
    return integrate_ode_rk45(ode, init_cond, 40.0, solve_time, parms, {0.0}, {0});
   }
}
