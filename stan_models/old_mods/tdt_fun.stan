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

       // y0 is the solution for the initial conditions for the age at BMT
       if (ageAtBMT > 40.0){
         y0 = solve_shm(ageAtBMT, init_cond, parms);
       } else {
         y0 = init_cond;
       }

       // init conditions at the BMT
       init_tb[1] = 0;                                           //at tbmt - # donor is zero
       init_tb[2] = 0;                                           //at tbmt - # donor is zero
       init_tb[3] = y0[1] + y0[3];                               //at tbmt - all ki67Hi cells are host
       init_tb[4] = y0[2] + y0[4];                               //at tbmt - all ki67Lo cells are host

       params[1:4] = parms[1:4];
       params[5] = ageAtBMT;                                           // age at BMT

       y_solve = to_array_1d(integrate_ode_rk45(shm, init_tb, ageAtBMT, rep_array(solve_time, 1), params, {0.0}, {0}));

       return y_solve;
     }

   real[] solve_ode_chi(real[] solve_time,
     real[] ageAtBMT,
     real[] init_cond,
     real[] parms){

       int numdim = size(solve_time);
       real y_solve[4*numdim];
       real chi_counts_mean;
       real host_counts_mean;
       real donor_counts_mean;
       real host_ki_mean;
       real donor_ki_mean;

       for (i in 1:numdim) {
         y_solve[4*i-3:4*i] = solve_chi(solve_time[i], ageAtBMT[i], init_cond, parms);

       }

       return y_solve;
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
