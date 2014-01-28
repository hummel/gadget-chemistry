#ifdef CHEMCOOL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"
#include "f2c.h"

/*initialization of chemcool*/
void chemcool_init(void)
{
  if(ThisTask==0)
    {
      printf("initialize cooling functions...\n");
      fflush(stdout);
    }

  COOLINMO();
  CHEMINMO();
  INIT_TOLERANCES();
  LOAD_H2_TABLE();

  if(ThisTask==0)
    {
      printf("initialization of cooling functions finished.\n");
      fflush(stdout);
    }

}

double initial_electron_fraction(void)
{
  double abe, abhp, abcp, absip, abop, abdp, abhep, abhcop;
  double abhepp, abmgp, absipp, abch3p;

  abe = 0.0;

  switch(All.ChemistryNetwork) {
  case 1:
    abhp   = All.InitHPlusAbund;
    abdp   = All.InitDIIAbund;
    abhep  = All.InitHeIIAbund;
    abhepp = All.InitHeIIIAbund;
    abe    = abhp + abdp + abhep + 2.0 * abhepp;
    break;
  case 2:
    abhp   = All.InitHPlusAbund;
    abdp   = All.InitDIIAbund;
    abhep  = All.InitHeIIAbund;
    abhepp = All.InitHeIIIAbund;
    abcp   = All.InitCIIAbund;
    absip  = All.InitSiIIAbund;
    abop   = All.InitOIIAbund;
    abe    = abhp + abdp + abhep + 2.0 * abhepp;
    abe    += abcp + absip + abop;
    break;
  case 3:
    abhp   = All.InitHPlusAbund;
    abdp   = All.InitDIIAbund;
    abhep  = All.InitHeIIAbund;
    abhepp = All.InitHeIIIAbund;
    abcp   = All.InitCIIAbund;
    absip  = All.InitSiIIAbund;
    abop   = All.InitOIIAbund;
    abhcop = All.InitHCOPlusAbund;
    absipp = All.InitSiIIIAbund;
    abch3p = All.InitCH3PlusAbund;
    abmgp  = All.InitMgPlusAbund; /* Set to zero currently */
    abe    = abhp + abdp + abhep + 2.0 * abhepp;
    abe    += abcp + absip + abop;
    abe    += abhcop + 2.0 * absipp + abch3p + abmgp;
    break;
 case 4:
 case 5:
    abhp  = All.InitHPlusAbund;
    abe   = abhp;
    break;
 case 7:
    abhp   = All.InitHPlusAbund;
    abhep  = All.InitHeIIAbund;
    abcp   = All.InitCIIAbund;
    abop   = All.InitOIIAbund;
    abhcop = All.InitHCOPlusAbund;
    abch3p = All.InitCH3PlusAbund;
    abe    = abhp + abhep + abcp + abop + abhcop + abch3p;
    break;
 case 8:
 case 10:
   abe = 0.0;
   break;
 default:
   break;
  }
  return abe;
}

double compute_electron_fraction(FLOAT abundances[NSPEC])
{
  double abe;

  switch(All.ChemistryNetwork) {
  case 1:
    abe = abundances[IHP] + abundances[IDP] + abundances[IHEP] 
        + 2.0 * abundances[IHEPP];
    break;
  case 2:
    abe = abundances[IHP] + abundances[IDP] + abundances[IHEP]
        + 2.0 * abundances[IHEPP] + abundances[IC] + abundances[IO]  
        + abundances[ISi];
    break;
  case 3:
    abe = abundances[IHP] + abundances[IDP] + abundances[IHEP] 
        + 2.0 * abundances[IHEPP] +  abundances[IC] + abundances[IO]  
        + abundances[ISi]  + abundances[IHCOP]
        + 2.0 * abundances[ISIPP] + abundances[ICH3P];
    break;
  case 4:
  case 5:
    abe = abundances[IHP];
    break;
  case 7:
    abe = abundances[IHP] + abundances[IHEP]
        +  abundances[IC] + abundances[IO]  + abundances[IHCOP]
        + abundances[ICH3P];
    break;
  case 8:
  case 10:
    abe = 0.0;
    break;
  default:
    abe = 0.0; 
    break;
  }
  return abe;
}

double compute_initial_gamma(void)
{
  double abh2, abe, gamma;

  /* Simple estimate of starting gamma; this should be sufficiently
   * accurate to get us going, provided that the initial H2 fraction
   * is small
   */
  abh2 = All.InitMolHydroAbund;
  abe  = initial_electron_fraction();
  gamma = (5.0 + 5.0 * ABHE - 3.0 * abh2 + 5.0 * abe) / 
          (3.0 + 3.0 * ABHE - abh2 + 3.0 * abe);
  return gamma;
}

/* Computes initial molecular weight in units of PROTONMASS */
double compute_initial_molecular_weight(void)
{
  double abe, abh2, abheI, abheII, abheIII, abhI, abhp, mu;

  /* Ignore minor corrections due to heavy elements, deuterium */
  switch(All.ChemistryNetwork) {
  case 1:
  case 2:
  case 3:
    abhp = All.InitHPlusAbund;
    abh2 = All.InitMolHydroAbund;
    abhI = 1.0 - abhp - 2.0 * abh2;

    abheIII = All.InitHeIIIAbund;
    abheII  = All.InitHeIIAbund;
    abheI   = ABHE - abheII - abheIII;

    abe = abhp + abheII + 2.0 * abheIII;
    break;
  case 7:
    abhp = All.InitHPlusAbund;
    abh2 = All.InitMolHydroAbund;
    abhI = 1.0 - abhp - 2.0 * abh2;

    abheIII = 0.0;
    abheII  = All.InitHeIIAbund;
    abheI   = ABHE - abheII;

    abe = abhp + abheII;
    break;
 case 4:
 case 5:
    abhp = All.InitHPlusAbund;
    abh2 = All.InitMolHydroAbund;
    abhI = 1.0 - abhp - 2.0 * abh2;

    abheIII = 0.0;
    abheII  = 0.0;
    abheI   = ABHE;

    abe = abhp;
    break;
  case 8:
  case 10:
    abhp = 0.0;
    abe  = 0.0;
    abh2 = All.InitMolHydroAbund;
    abhI = 1.0 - 2.0 * abh2;
    abheIII = 0.0;
    abheII  = 0.0;
    abheI   = ABHE;
    break;
  default:
    abe = abhp = abhI = abh2 = abheI = abheII = abheIII = 0.0;
    break;
  }

  mu = abhI + abhp + 2.0 * abh2 + 4.0 * (abheI + abheII + abheIII);
  mu /= abhI + abhp + abh2 + abheI + abheII + abheIII + abe;

  return mu;
}
#endif /* CHEMCOOL */
