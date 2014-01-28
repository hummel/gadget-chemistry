#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "allvars.h"
#ifdef JH_HEATING
#include "proto.h"
#endif

#define pi 3.1415927
#define c 2.99792458e10
#define h_nu 6.6262e-27
#define h_eV 4.13567e-15
#define k_B 1.3806e-16
#define pc 3.085678e18

#ifdef JH_HEATING
void initialize_heat_ion_rates()
{
  int i;
  double z, J0;
  J0 = All.xrbIntensity;

#ifdef JH_VARIABLE_HEATING
  /* Variable background ramping up from high z */
  z = 1.0 / (All.Time) - 1;
  i = 0;
  do
    {
      J0 = All.Jxr[i];
      i++;
    }
  while(All.Jz[i] > z);
  J0 = J0 * All.xrbIntensity;
#endif /* JH_VARIABLE_HEATING */

  calculate_heat_ion_rates(0, J0);
  calculate_heat_ion_rates(1, J0);
  calculate_heat_ion_rates(2, J0);
  All.heat_ion[6] = 0.0; // No LW Background!
  
  if(ThisTask==0)
    {
      printf("\nz: %lg  J0: %lg\n", z, J0);
      printf("Ionization Rates:\n   HI = %lg \n   HeI = %lg \n   HeII = %lg \n",
	     All.heat_ion[3], All.heat_ion[4], All.heat_ion[5]);
      printf("Heating Rates:\n   HI = %lg \n   HeI = %lg \n   HeII = %lg \n",
	     All.heat_ion[0], All.heat_ion[1], All.heat_ion[2]);
    } 
}

void calculate_heat_ion_rates(int rad_type, double J_0)
{
  int i = 0;
  int N_i_steps = 100;
  double A_0 = 6.3e-18;
  double nu_min, nu_max, nu_ion, nu_0;
  double ion_rate = 0, heat_rate = 0;
  double Z, logvmin, logvmax, Freq, epsilon, sigma;
  double F_nu, Freq_start, Freq_end, DFreq;

  // Energy integration limits in eV
  double E_0 = 1.0e3;
  double E_min = 1.0e3;
  double E_max = 10.0e3;
  // Atomic Number
  double Z_HI = 1.0;
  double Z_HeI = 0.89;
  double Z_HeII = 2.0;
  // Ionization Thresholds in Hz
  double nu_ion_HI = 3.3e15;
  double nu_ion_HeI = 5.95e15;
  double nu_ion_HeII = 1.32e16;
  
  
  if(rad_type == 0)
    {
      Z = Z_HI;
      nu_ion = nu_ion_HI;
    }
  if(rad_type == 1)
    {
      Z = Z_HeI;
      nu_ion = nu_ion_HeI;
    }
  if(rad_type == 2)
    {
      Z = Z_HeII;
      nu_ion = nu_ion_HeII;
    }
  nu_0 = E_0 / h_eV;
  nu_min = E_min / h_eV;
  nu_max = E_max / h_eV;
  logvmin = log10(nu_min);
  logvmax = log10(nu_max);
  if(nu_min < nu_ion)
    {
      logvmin = log10(nu_ion);
    }
  for(i = 0; i < N_i_steps; i++)
    {
      Freq = (logvmax - logvmin) / (double)(N_i_steps) * (i + 0.5) + logvmin;
      Freq_start = (logvmax - logvmin) / (double)(N_i_steps) * (i) + logvmin;
      Freq_end   = (logvmax - logvmin) / (double)(N_i_steps) * (i + 1.0) + logvmin;
      Freq       = pow(10, Freq);
      Freq_start = pow(10, Freq_start);
      Freq_end   = pow(10, Freq_end);
      DFreq      = Freq_end - Freq_start;
      
      F_nu = 4 * pi * J_0 * pow( Freq/nu_0, -1.5);
      epsilon = sqrt(Freq / nu_ion - 1.0);
      sigma = A_0 / pow(Z,2.0) * pow(nu_ion/Freq,4.0) 
	* exp(4.0 - (4.0*atan(epsilon) / epsilon)) / (1.0-exp(-2.0*PI / epsilon));
      ion_rate += F_nu * sigma / (h_nu * Freq) * DFreq;
      heat_rate += F_nu * sigma * ( 1.0 - nu_ion / Freq ) * DFreq;
    }
  if(rad_type == 0)
    {
      All.heat_ion[3] = ion_rate; // HI ion
      All.heat_ion[0] = heat_rate; // HI heat
    }
  if(rad_type == 1)
    {
      All.heat_ion[4] = ion_rate; // HeI ion
      All.heat_ion[1] = heat_rate; // HeI heat
    }
  if(rad_type == 2)
    {
      All.heat_ion[5] = ion_rate; // HeII ion
      All.heat_ion[2] = heat_rate; // HeII heat
    }
}
#endif /* JH_HEATING */

#ifdef RAYTRACE_TG
double heat_ion_rates(int rad_type, double L3, double T3)
  {
    int i = 0;
    int N_steps = 10000;
    int flag_sun = 0;
    double A_0 = 6.3e-18;
    double sigma = 2.0*pow(pi,5.0)*pow(k_B,4.0)/15.0/pow(h_nu,3.0)/pow(c,2.0);

    double L0 = pow(10.0, 5.568)*3.827e33;
    double L1 = pow(10.0, 6.095)*3.827e33;
    double L2 = pow(10.0, 6.574)*3.827e33;
    double T0 = pow(10.0, 4.922);
    double T1 = pow(10.0, 4.975);
    double T2 = pow(10.0, 4.999);

    double L = 0;
    double T_eff = 0;
    double prefactor = 0;
    double Z_HI = 1.0;
    double Z_HeI = 0.89;
    double Z_HeII = 2.0;
    double k_LW = 1.1e8;
    double nu_L_LW = 2.7e15;
    double nu_L_HI = 3.3e15;
    double nu_L_HeI = 5.95e15;
    double nu_L_HeII = 1.32e16;
    double nu_max = 1.0e20;
    double heat_HI = 0.0;
    double heat_HeI = 0.0;
    double heat_HeII = 0.0;
    double ion_LW = 0.0;
    double ion_HI = 0.0;
    double ion_HeI = 0.0;
    double ion_HeII = 0.0;
    double nu_LW[N_steps];
    double nu_HI[N_steps];
    double nu_HeI[N_steps];
    double nu_HeII[N_steps];
    double dnu_LW[N_steps];
    double dnu_HI[N_steps];
    double dnu_HeI[N_steps];
    double dnu_HeII[N_steps];
    double epsilon_HI[N_steps];
    double epsilon_HeI[N_steps];
    double epsilon_HeII[N_steps];
    double sigma_HI[N_steps];
    double sigma_HeI[N_steps];
    double sigma_HeII[N_steps];
    double B_LW[N_steps];
    double B_HI[N_steps];
    double B_HeI[N_steps];
    double B_HeII[N_steps];
    double x, y;

    int num=7;
    double heat_ion_vec[num];
 
    if(All.ray_flag_sun == 0)
      {
        L = L0;
        T_eff = T0;
      }

    if(All.ray_flag_sun == 1)
      {
        L = L1;
        T_eff = T1;
      }

    if(All.ray_flag_sun == 2)
      {
        L = L2;
        T_eff = T2;
      }

    if(All.ray_flag_sun == 3)
      {
        L = L3;
        T_eff = T3;
      }

    if(L > 0)
     prefactor = L/4.0/sigma/pow(T_eff,4.0)/pow(pc,2.0);
    else
     prefactor = 0;    

    for(i = 0; i < N_steps; i++)
      {
        nu_LW[i] = nu_L_LW*pow(nu_L_HI / nu_L_LW, (i + 0.5) / (double)N_steps);
        nu_HI[i] = nu_L_HI*pow(nu_L_HeII / nu_L_HI, (i + 0.5) / (double)N_steps);
        nu_HeI[i] = nu_L_HeI*pow(nu_L_HeII / nu_L_HeI, (i + 0.5) / (double)N_steps);
        nu_HeII[i] = nu_L_HeII*pow(nu_max / nu_L_HeII, (i + 0.5) / (double)N_steps);

        dnu_LW[i] = nu_L_LW*pow(nu_L_HI / nu_L_LW, (i + 1) / (double)N_steps) - nu_L_LW*pow(nu_L_HI / nu_L_LW, i / (double)N_steps);
        dnu_HI[i] = nu_L_HI*pow(nu_L_HeII / nu_L_HI, (i + 1) / (double)N_steps) - nu_L_HI*pow(nu_L_HeII / nu_L_HI, i / (double)N_steps);
        dnu_HeI[i] = nu_L_HeI*pow(nu_L_HeII / nu_L_HeI, (i + 1) / (double)N_steps) - nu_L_HeI*pow(nu_L_HeII / nu_L_HeI, i / (double)N_steps);
        dnu_HeII[i] = nu_L_HeII*pow(nu_max / nu_L_HeII, (i + 1) / (double)N_steps) - nu_L_HeII*pow(nu_max / nu_L_HeII, i / (double)N_steps);

        epsilon_HI[i] = sqrt(nu_HI[i]/nu_L_HI-1.0);
        epsilon_HeI[i] = sqrt(nu_HeI[i]/nu_L_HeI-1.0);
        epsilon_HeII[i] = sqrt(nu_HeII[i]/nu_L_HeII-1.0);

        x= (nu_HeI[i]/3.286e15) - .4434;
        y= pow(x,2) + 4.563;

        sigma_HI[i] = A_0/pow(Z_HI,2.0)*pow(nu_L_HI/nu_HI[i],4.0)*exp(4.0-4.0*atan(epsilon_HI[i])/epsilon_HI[i])/(1.0-exp(-2.0*pi/epsilon_HI[i]));
        //sigma_HeI[i] = A_0/pow(Z_HeI,2.0)*pow(nu_L_HeI/nu_HeI[i],4.0)*exp(4.0-4.0*atan(epsilon_HeI[i])/epsilon_HeI[i])/(1.0-exp(-2.0*pi/epsilon_HeI[i]));
        sigma_HeI[i] = 9.492e-16*(pow(x-1,2) + 4.158)*pow(y,-1.953)*pow(1 + .825*pow(y,0.25),-3.188); 
        sigma_HeII[i] = A_0/pow(Z_HeII,2.0)*pow(nu_L_HeII/nu_HeII[i],4.0)*exp(4.0-4.0*atan(epsilon_HeII[i])/epsilon_HeII[i])/(1.0-exp(-2.0*pi/epsilon_HeII[i]));

        B_LW[i] = 2.0*h_nu*pow(nu_LW[i],3.0)/pow(c,2.0)/(exp(h_nu*nu_LW[i]/k_B/T_eff)-1.0);
        B_HI[i] = 2.0*h_nu*pow(nu_HI[i],3.0)/pow(c,2.0)/(exp(h_nu*nu_HI[i]/k_B/T_eff)-1.0);
        B_HeI[i] = 2.0*h_nu*pow(nu_HeI[i],3.0)/pow(c,2.0)/(exp(h_nu*nu_HeI[i]/k_B/T_eff)-1.0);
        B_HeII[i] = 2.0*h_nu*pow(nu_HeII[i],3.0)/pow(c,2.0)/(exp(h_nu*nu_HeII[i]/k_B/T_eff)-1.0);

        heat_HI += prefactor*B_HI[i]*sigma_HI[i]*(1.0-nu_L_HI/nu_HI[i])*dnu_HI[i];
        heat_HeI += prefactor*B_HeI[i]*sigma_HeI[i]*(1.0-nu_L_HeI/nu_HeI[i])*dnu_HeI[i];
        heat_HeII += prefactor*B_HeII[i]*sigma_HeII[i]*(1.0-nu_L_HeII/nu_HeII[i])*dnu_HeII[i];

        //ion_LW += prefactor*k_LW*B_LW[i]*dnu_LW[i];
        ion_HI += prefactor*B_HI[i]*sigma_HI[i]/h_nu/nu_HI[i]*dnu_HI[i];
        ion_HeI += prefactor*B_HeI[i]*sigma_HeI[i]/h_nu/nu_HeI[i]*dnu_HeI[i];
        ion_HeII += prefactor*B_HeII[i]*sigma_HeII[i]/h_nu/nu_HeII[i]*dnu_HeII[i];
      }
   
    ion_LW = prefactor*k_LW*B_LW[i-4830];

    if(All.NumCurrentTiStep == 1)
      printf("sigma_H0 = %lg, sigma_He0 = %lg \n", sigma_HI[i-100], sigma_HeI[i-100]);

/*
    printf("heat_HI = %g\n", heat_HI);
    printf("heat_HeI = %g\n", heat_HeI);
    printf("heat_HeII = %g\n", heat_HeII);

    printf("ion_LW = %g\n", ion_LW);
    printf("ion_HI = %g\n", ion_HI);
    printf("ion_HeI = %g\n", ion_HeI);
    printf("ion_HeII = %g\n", ion_HeII);
*/
    heat_ion_vec[0] = heat_HI;
    heat_ion_vec[1] = heat_HeI;
    heat_ion_vec[2] = heat_HeII;
    heat_ion_vec[3] = ion_HI;
    heat_ion_vec[4] = ion_HeI;
    heat_ion_vec[5] = ion_HeII;
    heat_ion_vec[6] = ion_LW;

   if(rad_type == 0) 
     return heat_HI;

   else if(rad_type == 1)
     return heat_HeI;
    
   else if(rad_type == 2)
     return heat_HeII;

   else if(rad_type == 3)
     return ion_HI;

   else if(rad_type == 4)
     return ion_HeI;

   else if(rad_type == 5)
     return ion_HeII;

   else if(rad_type == 6)
     return ion_LW;

   else
   {
    printf("rad_type not properly set!\n");
    exit(0);
   }
}

#endif
