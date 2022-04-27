functions{
    real theta_spline(real t, real spline_int) {
      real fit1;

      real theta_max;
      real t_peak;
      real s1;
      real s2;
      real s3;
      real alph;

      theta_max = 1.5E6;
      alph      = 0.5378948;
      t_peak    = 50.08298;
      s1        = 0.04673187;
      s2        = 3.754336E-3;
      s3        = 3.532062E-2;

      fit1 = theta_max * (1 - exp(-s1 * t)) * (1 - int_step(t - t_peak)) +
             (theta_max * (1 - exp(-s1 * t_peak)) * (alph * exp(-s2 * (t - t_peak)) + (1 - alph) * exp(-s3 * (t - t_peak)))) * int_step(t - t_peak);

      fit1 = spline_int * fit1;

      return fit1;
    }

    real ki_fit(real t) {
      real fit2;
      fit2 = exp(-0.03337899 * (t + 2.92554110)) + 0.13732103;
      return fit2;
    }

  real[] shm(real t,
    real[] y,
    real[] parms,
    real[] rdata,
    int[] idata) {

    real theta0;
    real lambda;
    real div_nai;
    real theta_lo;
    real ki;
    real dydt[2];

    theta0 = parms[1];
    lambda = parms[2];
    div_nai = parms[3];
    theta_lo = parms[4];
    ki       = parms[5];

    // ki hi
    dydt[1] = theta_spline(t, theta0) * ki_fit(t) + 2 * div_nai * y[2] + div_nai * y[1] -
              (ki + lambda) * y[1];
    // ki lo
    dydt[2] = theta_spline(t, theta_lo) * (1 - ki_fit(t)) + ki * y[1] - div_nai * y[2] - lambda * y[2];
    return dydt;
}


   // to calculate total counts (N) at each 'tb'
   // this function solves the ODE for each value of 'tb' using 'ta' as the starting point (t0) and converts the solution into a 1d array
   real[] foreach_init(real t, t0, real[] init_cond, real[] parms, real[] rdata){
     //first argument to integrate_ode_rk45 must be the name of the ODE function with signature (real, real[], real[], real[], int[])
     int idata[0];

     return to_array_1d(integrate_ode_rk45(shm, init_cond, t0, rep_array(tb, 1), parms, rdata, idata));
   }

   // ode solution at each 'tb' -- this serves as initial condition when solving at 't'
   real[,] solve_init(real[] solve_time, real[] init_cond, real[] parms, real[] rdata, int num_pred{
    real y_hat[num_pred, 4];

    y_hat[1] = init_cond;
    for (i in 2:num_pred){
      y_hat[i] = foreach_init(solve_time[i], solve_time[1], init_cond, parms, rdata);
    }

    return y_hat;
   }
 }
