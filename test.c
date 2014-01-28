/* this remains in this file. Please add this file to the list of 
 * files for compilation in the Makefile
 */

/* follow the evolution of a single particle at given overdensity delta 
 * adapted from bg_init_cool.c: DriverTestCoolThermalEvolution()
 * Andreas Pawlik, UT Austin, 2010
 */
#include "allvars.h"
#include "proto.h"
#include "f2c.h"
#include "chemcool_consts.h"

#define pi 3.1415927
#define m_H 1.6726e-24
#define k_B 1.3806e-16
#define Hubble 3.2407789e-18
#define G 6.672e-8
#define X 0.76
#define yr 3.1536e7
#define gamma 5.0/3.0
#define unit_length 3.085678e21
#define unit_mass 1.989e43
#define unit_energy 1.989e53

int load_data();
double apTestChemCool(double z);
int N_part;

double Time, zred, Hparam;
int ThisTask, NTask;

int main(int argc, char **argv)
{

    double value;
    char dir[500];
    char dir2[500];
    char buf[500];
    char buf2[500]; 
    int i, j, n, m;
    double min, max;
    double time_in_Myr, MeanWeight, u, elec;
    int files;
    FILE *outfile;
 
    files = 1;
    COOLR.phi_pah = 1.0;
    COOLI.iflag_ad = 3;
    COOLI.iflag_mn = 4;
    COOLI.iflag_3bh2a = 1;
    COOLI.iflag_3bh2b = 1;
    COOLI.iflag_h3pra = 1;
    chemcool_init();

    sprintf(dir, "/scratch/cerberus/d4/jhummel/stampede/vanilla/snapshot");
    sprintf(dir2, "./");


    n = 1900;
    sprintf(buf, "%s_%d", dir, n);
    printf("reading %d...\n", n);
    N_gas=load_data();
    printf("processing %d...\n", n);
    double Pos[3][N_gas], Vel[3][N_gas], Temp[N_gas], nh[N_gas], ne[N_gas];
    double TracAbund[N_gas][6];


    sprintf(buf, "%s_%d.cooling1", dir2, n);
    sprintf(buf2, "%s%d_cooling", dir2, n);

    for(n=0; n < N_gas; n++)
      {
	if(P[n].ID < 0)
	  continue;
	nh[n] = SphP[n].Density*unit_mass/pow(unit_length,3.0)*pow(Hparam,2.0)/pow(Time,3.0)*X/m_H;
	elec = SphP[n].TracAbund[1] + SphP[n].TracAbund[2] + SphP[n].TracAbund[4] + 2.0*SphP[n].TracAbund[5];
	MeanWeight = 4.0 / (3.*0.76+1.0+4.*0.76*elec) * m_H;
	/* convert internal energy to cgs units */
	u  = SphP[n].Entropy * pow(SphP[n].Density, SphP[n].Gamma) / (SphP[n].Gamma - 1.0);
	u = u * unit_energy/ unit_mass;
	Temp[n]= MeanWeight/k_B * (SphP[n].Gamma-1.0) * u;
     }
    printf("nh=%g temp=%g efrac=%g u=%g\n",nh[0], Temp[0], elec, u);

     for(n = 0; n < N_gas; n++)
       {
	 Pos[0][n] = P[n].Pos[0];
	 Pos[1][n] = P[n].Pos[1];
	 Pos[2][n] = P[n].Pos[2];
	 Vel[0][n] = P[n].Vel[0];
	 Vel[1][n] = P[n].Vel[1];
	 Vel[2][n] = P[n].Vel[2];
	 ne[n] = SphP[n].TracAbund[1] + SphP[n].TracAbund[2] + SphP[n].TracAbund[4] + 2.0*SphP[n].TracAbund[5];
	 TracAbund[n][0] = SphP[n].TracAbund[0];
	 TracAbund[n][1] = SphP[n].TracAbund[1];
	 TracAbund[n][2] = SphP[n].TracAbund[2];
	 TracAbund[n][3] = SphP[n].TracAbund[3];
	 TracAbund[n][4] = SphP[n].TracAbund[4];
	 TracAbund[n][5] = SphP[n].TracAbund[5];
       }
     
     double temp;
     double dl, divv;
     int part_id, this_task, cur_ti_step;
     double abundances[TRAC_NUM];
     //double H2_rates[30];
     double n_H, ntot;
     double n_H_old, temp_old;
     double abhp, abh2, abhd, abdp, abhep, abhepp, abe, abhI, abheI, abdI;
     double dtime_in_s, this_dtime_in_s;
     double gamma1;
     
     outfile = fopen(buf2, "w");
     for(n = 0; n < N_gas; n++)
       {
	 abundances[IH2] = TracAbund[n][0]; //All.InitMolHydroAbund;
	 abundances[IHP] = TracAbund[n][1]; //All.InitHPlusAbund;
	 abundances[IDP] = TracAbund[n][2]; //All.InitDIIAbund;
	 abundances[IHD] = TracAbund[n][3]; //All.InitHDAbund;
	 abundances[IHEP] = TracAbund[n][4]; //All.InitHeIIAbund;
	 abundances[IHEPP] = TracAbund[n][5]; //All.InitHeIIIAbund;
	 n_H = nh[n];
	 temp = Temp[n];
	 n_H_old = n_H;
	 temp_old = temp;
 
	 abhp = abundances[IHP];
	 abh2 = abundances[IH2];
	 abhd = abundances[IHD];
	 abdp = abundances[IDP];
	 abhep = abundances[IHEP];
	 abhepp = abundances[IHEPP];
	 abe = ne[n];
	 
	 
	 part_id = 1;
	 this_task = 1;
	 cur_ti_step = 1;
	 dtime_in_s = 1.0 * yr;
	 this_dtime_in_s = dtime_in_s;
	 dl = 2.6238e+18;
	 divv = 1.2064e-13;
	 gamma1 = (5. + 5. * ABHE - 3. * abh2 + 5. * abe) / (3. + 3. * ABHE - abh2 + 3. * abe);
	 ntot = (1. + ABHE - abh2 + abe) * n_H;
	 //u = temp / ((gamma1 - 1.0) * ntot * BOLTZMANN);
	 COOLR.redshift = zred;
	 printf("n=%d n_H=%g, temp=%g, u=%g H2=%g HP=%g DP=%g\n", n, n_H, temp, u, abundances[IH2], abundances[IHP], abundances[IDP]);

	 //EVOLVE_ABUNDANCES(&this_dtime_in_s, &dl, &n_H, &divv, &u, abundances, &cur_ti_step, &this_task, &part_id);
	 value = GetCoolTime(SphP, 0);
	 
	 //fprintf(outfile,"%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t \n", COOLR.H2_rates[0], COOLR.H2_rates[1], COOLR.H2_rates[2], COOLR.H2_rates[3], COOLR.H2_rates[4], COOLR.H2_rates[5], COOLR.H2_rates[6], COOLR.H2_rates[7], COOLR.H2_rates[8], COOLR.H2_rates[9], COOLR.H2_rates[10], COOLR.H2_rates[11], COOLR.H2_rates[12], COOLR.H2_rates[13], COOLR.H2_rates[14], COOLR.H2_rates[15], COOLR.H2_rates[16], COOLR.H2_rates[17], COOLR.H2_rates[18], COOLR.H2_rates[19], COOLR.H2_rates[20], COOLR.H2_rates[21], COOLR.H2_rates[22], COOLR.H2_rates[23], COOLR.H2_rates[24],COOLR.H2_rates[25],  COOLR.H2_rates[26], COOLR.H2_rates[27], COOLR.H2_rates[28], n_H_old, temp_old);
	 
	 
       }
     
     fclose(outfile);
     //free(P);
     //exit(0);
     
     
     
     
     // value = apTestChemCool(10.0);
     /*    
     outfile = fopen(buf2, "w");
     
     fwrite(&N_gas,sizeof(int),1,outfile);
     fwrite(&header.time,sizeof(double),1,outfile);
     fwrite(&header.redshift,sizeof(double),1,outfile);
     fwrite(&time_in_Myr,sizeof(double),1,outfile);
     fwrite(&center_x,sizeof(double),1,outfile);
     fwrite(&center_y,sizeof(double),1,outfile);
     fwrite(&center_z,sizeof(double),1,outfile);
     fwrite(Pos,sizeof(Pos),1,outfile);
     fwrite(Temp,sizeof(Temp),1,outfile);
     fwrite(nh,sizeof(nh),1,outfile);
     fwrite(ne,sizeof(ne),1,outfile);
     
     fclose(outfile);
     */
     free(P);
     free(SphP);
     
     printf("done2!\n");
     return 0;
}


int load_data(void)
  {
    All.MaxPart = 1;
    All.MaxPartSph = 1;
    allocate_memory();

    P[0].Pos[0] = 50.26185814;
    P[0].Pos[1] = 50.24031015;
    P[0].Pos[2] = 49.61333174;
    P[0].ID = 3691688;
    P[0].Mass = 1.0339542807005775e-12;
    SphP[0].Gamma = 1.6664895908462534;
    SphP[0].Entropy = 26998505.46;
    SphP[0].Density = 5.5454793718537619e-05;
    SphP[0].Hsml = 0.0058222455615665521;
    SphP[0].TracAbund[0] = 3.88116604e-004;
    SphP[0].TracAbund[1] = 7.24870648e-006;
    SphP[0].TracAbund[2] = 1.60593500e-010;
    SphP[0].TracAbund[3] = 5.52034755e-008;
    SphP[0].TracAbund[4] = 1.87174183e-085;
    SphP[0].TracAbund[5] = 1.28703248e-239;
    Time= 0.038405460020081036;
    zred= 25.037964379989997;
    Hparam = 0.69999999999999996;
    return(1);
  }

