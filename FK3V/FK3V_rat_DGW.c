/******************************************************************************
Fenton Karma Rat Model by Dominic Whittaker

Please cite the following if the code is used:

Investigation of the role of myocyte orientations in cardiac arrhythmia using image-based models
Whittaker, D. G., Benson, A. P., Teh, I., Schneider, J. E., Colman, M. A.
Biophysical Journal (2019), 117, 1-13 (in press).
https://doi.org/10.1016/j.bpj.2019.09.041
******************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

#define FIRST_ACT 25
#define SECOND_ACT 75

const double dt = 0.01;
const double dx = 0.1;
// const double BCL = 200;
const double stim_dur = 1.0;
const double stim_amp = -0.2;
const int paceAP = 3;
const int NX = 100;
const int APD_rec_site = 50;

double V_0 = -82;
double V_fi = 35;
double V_thresh = -70.0;

double BCL = 200;

int heav (double x);

int main() {

  // Files
  FILE *output_V;
  FILE *output_CV;
  FILE *output_APD; 
  output_V = fopen("FK3V_V_1D_LV.dat", "w");
  output_CV = fopen("FK3V_CV_1D_LV.dat", "w");
  output_APD = fopen("FK3V_APD_1D_LV.dat", "w");

  // Variables
  int m, p, q, r;
  double tau_v_minus, tau_so, CV;
  double dw_dt, dv_dt, dd_dt, du_dt;
  double v[NX], w[NX], d[NX], u[NX], tempu[NX], V[NX], VP[NX], lap[NX];
  double I_fi[NX], I_so[NX], I_si[NX], I_ion[NX], I_stim, cell_stim;
  double t, t_end, stim_time;
  int i, j, k, N;
  long int n;
  bool APD_flag;

  int N_BCL = BCL / dt;
  int N_stim = stim_dur / dt;

  // Restitution protocol parameters
  int S2[] = {80,100,120,140,160,180,200}; // Use this to output APD restitution
  // int S2[] = {200};
  int N_S2;
  int size = sizeof(S2)/sizeof(S2[0]);
  int count = 0;

  // Factors from optimisation to Benoist APD restitution data
  double f[3] = {1.087619, 2.024560, 0.979126};
  double g[3] = {0.998880, 1.146296, 0.960981};

  // Restitution protocol loop
  for(i = 0; i < size; i++) {

    t = 0.0;
    n = 0;
    N = -1;
    N_S2 = S2[i] / dt;
    stim_time = 0.0;
    t_end = paceAP * BCL + 1000.0;

    // APD & CV variables
    double V_max = -100.0, V_low = -80.26, V_90;
    double APD, APD_end, APD_start = 0.0, obj_APD;
    double CV_start = 0.0;

    // Initialize state variables
    for (j = 0; j < NX; j++) {
      v[j] = 1.0;
      w[j] = 1.0;
      d[j] = 0.0;
      u[j] = 0.0;
      tempu[j] = 0.0;
      V[j] = 0.0;
      VP[j] = 0.0;
    }

    // See original paper for different FK variants
    double u_c          = 0.13;
    double u_v          = 0.04;
    double tau_v1_minus = g[0]*f[0]*50;//19.6;
    double tau_v2_minus = g[1]*f[1]*200;//1250;
    double tau_v_plus   = 3.33;
    double tau_w_minus  = 41;
    double tau_w_plus   = 870;
    double u_csi        = 0.85;
    double tau_d        = 0.1;//0.25;
    double tau_o        = 12.5;
    double tau_si       = g[2]*f[2]*45;//29;
    double xk           = 10;
    double tau_r        = 33.33; //+20% BASE; -20% RV
    double D            = 0.1555;

    // Time loop
    while (t < t_end) {

      // Stimlulus Current
      if(n == stim_time){
        N++;
        printf("Stim at time %f\n", t);
      }
      if(n - stim_time <= N_stim && n - stim_time > 0)
        I_stim = stim_amp;
      else 
        I_stim = 0.0;

      // Cable loop
      for (j = 0; j < NX; j++) {

        cell_stim = 0.0;
        if (j < 10)
          cell_stim = I_stim;

        tau_v_minus = heav(u[j] - u_v) * tau_v1_minus + heav(u_v - u[j]) * tau_v2_minus;
        
        // State Variables
        dv_dt = (heav(u_c - u[j]) * (1.0 - v[j])) / tau_v_minus - (heav(u[j] - u_c) * v[j]) / tau_v_plus;
        dw_dt = (heav(u_c - u[j]) * (1.0 - w[j])) / tau_w_minus - (heav(u[j] - u_c) * w[j]) / tau_w_plus;

        v[j] += dt * dv_dt;  //fast gate
        w[j] += dt * dw_dt;  //slow gate
                
        // Currents
        I_fi[j] = (-v[j] / tau_d) * heav(u[j] - u_c) * (1.0 - u[j]) * (u[j] - u_c);
        I_so[j] = (u[j] / tau_o) * heav(u_c - u[j]) + (heav(u[j] - u_c)) / tau_r;
        I_si[j] = (-w[j] / (2.0 * tau_si)) * (1.0 + tanh(xk * (u[j] - u_csi)));

        I_ion[j] = - (I_fi[j] + I_so[j] + I_si[j] + cell_stim);

      } // End cable loop

      int cc, y;

        // Boundary conditions!
        for (y = 0; y < NX; y++)  {
          if(y==0) lap[y] = (u[y + 1] - u[y]);
          else if(y==(NX-1)) lap[y] = (u[y - 1] - u[y]);
          else lap[y] = (u[y + 1] + u[y - 1] - 2.0 * u[y]);
        }

        for (y = 0; y < NX; y++)  {
          u[y] += dt*(I_ion[y] + ((D / (dx * dx) * lap[y])));
          V[y] = V_0 + u[y]*(V_fi - V_0); // Transform voltages into real units
        }
      
      // CV calculation
      if (V[FIRST_ACT]>=-40. && VP[FIRST_ACT]<-40.) CV_start = t;
      if (V[SECOND_ACT]>=-40. && VP[SECOND_ACT]<-40.) {
        CV = 100 * (SECOND_ACT - FIRST_ACT) * dx / (t - CV_start);
      }
      
      // APD90 calculation
      if (n == stim_time - 1){
        V_max = -100.0;
        APD_flag = true;
      }
      if (V[APD_rec_site] > V_max) {
        V_max = V[APD_rec_site];
        V_90 = V_max - (V_max - V_0)*0.90;
        if(V[APD_rec_site] > V_thresh && APD_flag) {
          APD_start = t;
          APD_flag = false;
        }
      }
      double baseline = V_90;
      if (V[APD_rec_site] <= baseline && VP[APD_rec_site] > baseline && n - stim_time > 50) {
        if(N < paceAP)
          stim_time = (N + 1) * BCL / dt;
        else if(N == paceAP){
          stim_time = (N * BCL + S2[i])/ dt;
        }
        APD = t - APD_start;
        APD_end = t;
        count++;
      }

      // Update voltage!
      for (k = 0; k < NX; k++) {
        VP[k] = V[k];
      }

      if (n % 100 == 0 && S2[i] == 200) { // Only plot for S2 of 200 ms
      // if (n % 100 == 0) {
        for (j = 0; j < NX; j++) {
          fprintf(output_V, "%f\t", u[j]); // Output voltages along 1D strand
        }
        fprintf(output_V, "\n");
      }

      n++;
      t = n * dt;

    } // End time loop

    printf("APD = %0.1f ms\n", APD);
    printf("CV = %0.1f cm/s\n", CV);
    fprintf(output_APD, "%f\t%d\t%f\n", S2[i] - APD, S2[i], APD); // Plot using either 1:3 or 2:3
    fprintf(output_CV,"\t%f\t%f\n", S2[i] - APD, CV); 
    
  } // End restitution loop

  fclose(output_V);
  fclose(output_CV);
  fclose(output_APD);

  return 0;

} // End main         

int heav (double x) {
  if (x >= 0) 
    return 1;
  else
    return 0;
}