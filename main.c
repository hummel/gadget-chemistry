/* this remains in this file. Please add this file to the list of 
 * files for compilation in the Makefile
 */

/* follow the evolution of a single particle at given overdensity delta 
 * adapted from bg_init_cool.c: DriverTestCoolThermalEvolution()
 * Andreas Pawlik, UT Austin, 2010
 */
#include "omp.h"
#include "allvars.h"
#include "proto.h"
#include "f2c.h"
#include "chemcool_consts.h"

#define pi 3.1415927
#define m_H 1.6726e-24
#define k_B 1.3806e-16
#define G 6.672e-8
#define X 0.76
#define yr 3.1536e7
#define unit_length 3.085678e21
#define unit_mass 1.989e43
#define unit_energy 1.989e53

int load_data();
int N_part;
int ThisTask, NTask;

int main(int argc, char **argv)
{
    FILE *CoolTime;
    char dir[500];
    char dir2[500];
    char buf[500];
    char buf2[500]; 
    int  n, i, nsp, snapshot, ipar[NIPAR];
    double mu, h2frac, u, density, n_H;
    double energy, temp, divv, dl, dtcool;
    double rpar[NRPAR], y[NSPEC], ydot[NSPEC];
    double abundances[TRAC_NUM];
    double t_start;

    if(argc < 2)
      {
	printf("Parameters are missing.\n");
	printf("Call with <ParameterFile> [<SnapshotNumber>]\n");
	exit(0);
      }

    strcpy(ParameterFile, argv[1]);

    if(argc >= 3)
      snapshot = atoi(argv[2]);
    else
      {
	printf("Parameters are missing.\n");
	printf("Call with <ParameterFile> [<SnapshotNumber>]\n");
	exit(0);
      }


    NTask = 1;
    read_parameter_file(ParameterFile);
    allocate_commbuffers();
    set_units();

    COOLR.phi_pah = 1.0;
    COOLI.iflag_ad = 3;
    COOLI.iflag_mn = 4;
    COOLI.iflag_3bh2a = 1;
    COOLI.iflag_3bh2b = 1;
    COOLI.iflag_h3pra = 1;

    /* Initialize Cooling Functions */
    chemcool_init();

#ifdef JH_HEATING
    initialize_heat_ion_rates();

    for(i=0; i<=6; i++)
      {
	COOLR.heat_ion[i] = All.heat_ion[i];
      }
#endif


    sprintf(buf, "%s/snapshot_%03d", All.InitCondFile, snapshot);
    read_ic(buf);
    //N_gas = load_data();
    printf("processing...\n");
    nsp = NSPEC;
    t_start = -1.; /* Nothing in rate_eq depends on t_start, so it doesn't
		      matter what value we give it */

    sprintf(buf, "%s%s_%04d.dat", All.OutputDir, All.SnapshotFileBase, snapshot);
    CoolTime=fopen(buf,"w");
    #pragma omp parallel for private(n,i,y,density,h2frac,mu,u,temp,energy, \
				     n_H,dl,divv,rpar,ipar,ydot,dtcool)
    for(n = 0; n < N_gas; n++)
      {
	if(P[n].Mass < 2e-12)
	  {
	    for (i=0; i<TRAC_NUM; i++)
	      {
		y[i] = SphP[n].TracAbund[i];
	      }
	    density = SphP[n].Density * unit_mass / pow(unit_length, 3.0) 
	      * pow(All.HubbleParam, 2.0) / pow(All.Time, 3.0);
	    /*
	       The value stored on the disk is particle internal energy, not entropy.
	       This is updated by Gadget eventually, but we're not letting it get that far here.
	       In the event that SphP[n].Entropy ACTUALLY contains the entropy, this should be:
	       u  = SphP[n].Entropy * pow(SphP[n].Density, SphP[n].Gamma) / (SphP[n].Gamma - 1.0);
	    */
	    u = SphP[n].Entropy;
	    u = u * unit_energy / unit_mass; /* convert internal energy to cgs units */
	    energy = density * u; /* Convert from mass specific energy to volume specific energy. */
	    y[ITMP] = energy;
	    
	    n_H = density * X/m_H;
	    dl = SphP[n].Hsml * unit_length;
	    divv = 0.0;

	    h2frac = y[IH2] * 2.0;
	    mu = 1.0 / ((0.24/4.0) + ((1.0-h2frac)*0.76) + (h2frac*.76/2.0));
	    temp = mu * m_H / k_B * (SphP[n].Gamma-1.0) * u;
	    //printf("n=%d n_H=%g density=%g, temp=%g internal energy=%g\n",
	    //n, n_H, density, temp, u, h2frac);

	    rpar[0] = n_H; // hydrogen number density
	    rpar[1] = dl; // smoothing length
	    rpar[2] = divv; // divergence ofthe velocity field
	    ipar[0] = 0; // ???
	    
	    RATE_EQ(&nsp, &t_start, y, ydot, rpar, ipar);
	    
	    if (ydot[ITMP] == 0.0)
	      {
		/* Cooling time is formally infinite. Since we can't return infinity,
		   however, we make do with a very big number: 10^20 seconds. */
		dtcool = 1e20;
	      }
	    else
	      {
		/* We assume that the energy is non-zero */
		dtcool = y[ITMP] / ydot[ITMP];
	      }
	    if(dtcool < 0) 
	      dtcool *= -1; /* make sure timestep is not negative */
	    
	    fprintf(CoolTime,"%15.11g %8d %15.6g %15.11g %15.11g\n", 
		    All.Time, P[n].ID, n_H, temp, dtcool);
	  }
      }
    
    fclose(CoolTime);
    free(P);
    free(SphP);
    
    printf("done!\n");
    return 0;
}

int load_data(void)
{
  All.MaxPart = 4;
  All.MaxPartSph = 4;
  allocate_memory();

  P[0].Pos[0] = 50.26185814;
  P[0].Pos[1] = 50.24031015;
  P[0].Pos[2] = 49.61333174;
  P[0].Mass = 1.0339542807005775e-12;
  P[0].ID = 4011904;
  SphP[0].Gamma = 1.4523031457733506;
  SphP[0].Entropy = 2.75628e-6;
  SphP[0].Density = 15032.003627719097;
  SphP[0].Hsml = 8.5362500665953259e-06;
  SphP[0].TracAbund[0] = 3.40758378e-001;
  SphP[0].TracAbund[1] = 1.47742345e-010;
  SphP[0].TracAbund[2] = 2.97435944e-015;
  SphP[0].TracAbund[3] = 1.88810849e-005;
  SphP[0].TracAbund[4] = 8.68327877e-080;
  SphP[0].TracAbund[5] = 5.64397803e-197;

  P[1].Pos[0] = 50.26185814;
  P[1].Pos[1] = 50.24031015;
  P[1].Pos[2] = 49.61333174;
  P[1].Mass = 1.0339542807005775e-12;
  P[1].ID = 3691688;
  SphP[1].Gamma = 1.6664895908462534;
  SphP[1].Entropy = 26998505.46;
  SphP[1].Density = 5.5454793718537619e-05;
  SphP[1].Hsml = 0.0058222455615665521;
  SphP[1].TracAbund[0] = 3.88116604e-004;
  SphP[1].TracAbund[1] = 7.24870648e-006;
  SphP[1].TracAbund[2] = 1.60593500e-010;
  SphP[1].TracAbund[3] = 5.52034755e-008;
  SphP[1].TracAbund[4] = 1.87174183e-085;
  SphP[1].TracAbund[5] = 1.28703248e-239;

  P[2].Pos[0] = 50.26185814;
  P[2].Pos[1] = 50.24031015;
  P[2].Pos[2] = 49.61333174;
  P[2].Mass = 1.0339542807005775e-12;
  P[2].ID = 3678762;
  SphP[2].Gamma = 1.6662185785540509;
  SphP[2].Entropy = 0.37089861271490865;
  SphP[2].Density = 3.766345493332401;
  SphP[2].Hsml = 0.00014175038541804796;
  SphP[2].TracAbund[0] = 8.40674666e-004;
  SphP[2].TracAbund[1] = 1.20046114e-008;
  SphP[2].TracAbund[2] = 2.79837883e-013;
  SphP[2].TracAbund[3] = 4.91481098e-008;
  SphP[2].TracAbund[4] = 1.31852790e-116;
  SphP[2].TracAbund[5] = 8.55831086e-225;

  P[3].Pos[0] = 50.26185814;
  P[3].Pos[1] = 50.24031015;
  P[3].Pos[2] = 49.61333174;
  P[3].Mass = 1.0339542807005775e-12;
  P[3].ID = 3322633;
  SphP[3].Gamma = 1.6665151416569162;
  SphP[3].Entropy = 58579920.298570916;
  SphP[3].Density = 3.7917358113706579e-05;
  SphP[3].Hsml = 0.0064916396213253085;
  SphP[3].TracAbund[0] = 3.20500793e-004;
  SphP[3].TracAbund[1] = 8.31632314e-006;
  SphP[3].TracAbund[2] = 1.87643312e-010;
  SphP[3].TracAbund[3] = 3.67416796e-008;
  SphP[3].TracAbund[4] = 1.35279413e-088;
  SphP[3].TracAbund[5] = 8.74635849e-268;

  All.Time = 0.038405460020081036;
  All.HubbleParam = 0.7;
  return(4);
  }

/*! returns the maximum of two double
 */
double dmax(double x, double y)
{
  if(x > y)
    return x;
  else
    return y;
}
