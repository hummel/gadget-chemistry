#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <sys/types.h>
//#include <unistd.h>

#include "allvars.h"
#include "proto.h"
#include "f2c.h"


/*! \file begrun.c
 *  \brief initial set-up of a simulation run
 *
 *  This file contains various functions to initialize a simulation run. In
 *  particular, the parameterfile is read in and parsed, the initial
 *  conditions or restart files are read, and global variables are
 *  initialized to their proper values.
 */




/*! Computes conversion factors between internal code units and the
 *  cgs-system.
 */
void set_units(void)
{
  double meanweight;
#ifdef CHEMCOOL
  double gamm1;
#endif

  All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitTime_in_Megayears = All.UnitTime_in_s / SEC_PER_MEGAYEAR;

  if(All.GravityConstantInternal == 0)
    All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2);
  else
    All.G = All.GravityConstantInternal;

  All.UnitDensity_in_cgs = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
  All.UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
  All.UnitCoolingRate_in_cgs = All.UnitPressure_in_cgs / All.UnitTime_in_s;
  All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

  /* convert some physical input parameters to internal units */

  All.Hubble = HUBBLE * All.UnitTime_in_s;

  if(ThisTask == 0)
    {
      printf("\nHubble (internal units) = %g\n", All.Hubble);
      printf("G (internal units) = %g\n", All.G);
      printf("UnitMass_in_g = %g \n", All.UnitMass_in_g);
      printf("UnitTime_in_s = %g \n", All.UnitTime_in_s);
      printf("UnitVelocity_in_cm_per_s = %g \n", All.UnitVelocity_in_cm_per_s);
      printf("UnitDensity_in_cgs = %g \n", All.UnitDensity_in_cgs);
      printf("UnitEnergy_in_cgs = %g \n", All.UnitEnergy_in_cgs);
      printf("\n");
    }

#ifdef CHEMCOOL
  meanweight = compute_initial_molecular_weight();
  gamm1      = compute_initial_gamma() - 1.0;
  All.MinEgySpec = 1 / meanweight * (1.0 / gamm1) * (BOLTZMANN / PROTONMASS) * All.MinGasTemp;
  All.MinEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
#else /* CHEMCOOL */
  meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: we assume neutral gas here */
#ifdef ISOTHERM_EQS
  All.MinEgySpec = 0;
#else /* ISOTHERM_EQS */
#ifdef POLYTROPE
  All.MinEgySpec = 0;
#else /* POLYTROPE */
  All.MinEgySpec = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.MinGasTemp;
  All.MinEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
#endif /* POLYTROPE */
#endif /* ISOTHERM_EQS */
#endif /* CHEMCOOL */

}




/*! This function parses the parameterfile in a simple way.  Each paramater
 *  is defined by a keyword (`tag'), and can be either of type double, int,
 *  or character string.  The routine makes sure that each parameter
 *  appears exactly once in the parameterfile, otherwise error messages are
 *  produced that complain about the missing parameters.
 */
void read_parameter_file(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd, *fdout;
  char buf[200], buf1[200], buf2[200], buf3[400];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int  errorFlag = 0;


  if(sizeof(long long) != 8)
    {
      if(ThisTask == 0)
	printf("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
      exit(0);
    }

  if(sizeof(int) != 4)
    {
      if(ThisTask == 0)
	printf("\nType `int' is not 32 bit on this platform. Stopping.\n\n");
      exit(0);
    }

  if(sizeof(float) != 4)
    {
      if(ThisTask == 0)
	printf("\nType `float' is not 32 bit on this platform. Stopping.\n\n");
      exit(0);
    }

  if(sizeof(double) != 8)
    {
      if(ThisTask == 0)
	printf("\nType `double' is not 64 bit on this platform. Stopping.\n\n");
      exit(0);
    }


  if(ThisTask == 0)		/* read parameter file on process 0 */
    {
      nt = 0;

      strcpy(tag[nt], "InitCondFile");
      addr[nt] = All.InitCondFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputDir");
      addr[nt] = All.OutputDir;
      id[nt++] = STRING;

      strcpy(tag[nt], "SnapshotFileBase");
      addr[nt] = All.SnapshotFileBase;
      id[nt++] = STRING;

      strcpy(tag[nt], "EnergyFile");
      addr[nt] = All.EnergyFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "CpuFile");
      addr[nt] = All.CpuFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "InfoFile");
      addr[nt] = All.InfoFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "TimingsFile");
      addr[nt] = All.TimingsFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "RestartFile");
      addr[nt] = All.RestartFile;
      id[nt++] = STRING;

      /* Ionizing Background */
#if defined(XRAY_BACKGROUND) || defined(COSMIC_RAY_BACKGROUND)
      strcpy(tag[nt], "HeatFile");
      addr[nt] = All.HeatFile;
      id[nt++] = STRING;
#endif

      /* X-ray background intensity */
#ifdef XRAY_BACKGROUND
      strcpy(tag[nt], "xrbIntensity");
      addr[nt] = &All.xrbIntensity;
      id[nt++] = DOUBLE;

#ifdef XRAY_VARIABLE_HEATING
      strcpy(tag[nt], "xrbFile");
      addr[nt] = All.xrbFile;
      id[nt++] = STRING;

#endif /* XRAY_VARIABLE_HEATING */
#endif /* XRAY_BACKGROUND */

      /* Cosmic ray background intensity */
#ifdef COSMIC_RAY_BACKGROUND
      strcpy(tag[nt], "crbIntensity");
      addr[nt] = &All.crbIntensity;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "CRheatPerInteraction");
      addr[nt] = &All.CR_heat;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "CRspectrum_min");
      addr[nt] = &All.CR_spectrum_min;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "CRspectrum_max");
      addr[nt] = &All.CR_spectrum_max;
      id[nt++] = DOUBLE;

#ifdef CR_VARIABLE_HEATING
      strcpy(tag[nt], "crbFile");
      addr[nt] = All.crbFile;
      id[nt++] = STRING;

#endif /* CR_VARIABLE_HEATING */
#endif /* COSMIC_RAY_BACKGROUND */


      /*SINK*/
      strcpy(tag[nt], "SinkFile");
      addr[nt] = All.SinkFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "ResubmitCommand");
      addr[nt] = All.ResubmitCommand;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListFilename");
      addr[nt] = All.OutputListFilename;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListOn");
      addr[nt] = &All.OutputListOn;
      id[nt++] = INT;

      strcpy(tag[nt], "Omega0");
      addr[nt] = &All.Omega0;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "OmegaBaryon");
      addr[nt] = &All.OmegaBaryon;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "OmegaLambda");
      addr[nt] = &All.OmegaLambda;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "HubbleParam");
      addr[nt] = &All.HubbleParam;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BoxSize");
      addr[nt] = &All.BoxSize;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "PeriodicBoundariesOn");
      addr[nt] = &All.PeriodicBoundariesOn;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeOfFirstSnapshot");
      addr[nt] = &All.TimeOfFirstSnapshot;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "CpuTimeBetRestartFile");
      addr[nt] = &All.CpuTimeBetRestartFile;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeBetStatistics");
      addr[nt] = &All.TimeBetStatistics;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeBegin");
      addr[nt] = &All.TimeBegin;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeMax");
      addr[nt] = &All.TimeMax;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeBetSnapshot");
      addr[nt] = &All.TimeBetSnapshot;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
      addr[nt] = &All.UnitVelocity_in_cm_per_s;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "UnitLength_in_cm");
      addr[nt] = &All.UnitLength_in_cm;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "UnitMass_in_g");
      addr[nt] = &All.UnitMass_in_g;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TreeDomainUpdateFrequency");
      addr[nt] = &All.TreeDomainUpdateFrequency;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ErrTolIntAccuracy");
      addr[nt] = &All.ErrTolIntAccuracy;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ErrTolTheta");
      addr[nt] = &All.ErrTolTheta;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ErrTolForceAcc");
      addr[nt] = &All.ErrTolForceAcc;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MinGasHsmlFractional");
      addr[nt] = &All.MinGasHsmlFractional;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MaxSizeTimestep");
      addr[nt] = &All.MaxSizeTimestep;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MinSizeTimestep");
      addr[nt] = &All.MinSizeTimestep;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MaxRMSDisplacementFac");
      addr[nt] = &All.MaxRMSDisplacementFac;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ArtBulkViscConst");
      addr[nt] = &All.ArtBulkViscConst;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "CourantFac");
      addr[nt] = &All.CourantFac;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "DesNumNgb");
      addr[nt] = &All.DesNumNgb;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MaxNumNgbDeviation");
      addr[nt] = &All.MaxNumNgbDeviation;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ComovingIntegrationOn");
      addr[nt] = &All.ComovingIntegrationOn;
      id[nt++] = INT;

      strcpy(tag[nt], "ICFormat");
      addr[nt] = &All.ICFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "SnapFormat");
      addr[nt] = &All.SnapFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesPerSnapshot");
      addr[nt] = &All.NumFilesPerSnapshot;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesWrittenInParallel");
      addr[nt] = &All.NumFilesWrittenInParallel;
      id[nt++] = INT;

      strcpy(tag[nt], "ResubmitOn");
      addr[nt] = &All.ResubmitOn;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfTimestepCriterion");
      addr[nt] = &All.TypeOfTimestepCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfOpeningCriterion");
      addr[nt] = &All.TypeOfOpeningCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeLimitCPU");
      addr[nt] = &All.TimeLimitCPU;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningHalo");
      addr[nt] = &All.SofteningHalo;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningDisk");
      addr[nt] = &All.SofteningDisk;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningBulge");
      addr[nt] = &All.SofteningBulge;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningGas");
      addr[nt] = &All.SofteningGas;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningStars");
      addr[nt] = &All.SofteningStars;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningBndry");
      addr[nt] = &All.SofteningBndry;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningHaloMaxPhys");
      addr[nt] = &All.SofteningHaloMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningDiskMaxPhys");
      addr[nt] = &All.SofteningDiskMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningBulgeMaxPhys");
      addr[nt] = &All.SofteningBulgeMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningGasMaxPhys");
      addr[nt] = &All.SofteningGasMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningStarsMaxPhys");
      addr[nt] = &All.SofteningStarsMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningBndryMaxPhys");
      addr[nt] = &All.SofteningBndryMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BufferSize");
      addr[nt] = &All.BufferSize;
      id[nt++] = INT;

      strcpy(tag[nt], "PartAllocFactor");
      addr[nt] = &All.PartAllocFactor;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TreeAllocFactor");
      addr[nt] = &All.TreeAllocFactor;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "GravityConstantInternal");
      addr[nt] = &All.GravityConstantInternal;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "InitGasTemp");
      addr[nt] = &All.InitGasTemp;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MinGasTemp");
      addr[nt] = &All.MinGasTemp;
      id[nt++] = DOUBLE;

#ifdef CHEMCOOL
      strcpy(tag[nt],"H2RefDustEff"); 
      addr[nt]=&All.H2RefDustEff;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"OxyAbund"); 
      addr[nt]=&All.OxyAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"CarbAbund"); 
      addr[nt]=&All.CarbAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"SiAbund"); 
      addr[nt]=&All.SiAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"DeutAbund"); 
      addr[nt]=&All.DeutAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"MgAbund"); 
      addr[nt]=&All.MgAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"UVField"); 
      addr[nt]=&All.UVField;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"PhiPAH"); 
      addr[nt]=&All.PhiPAH;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitDustTemp"); 
      addr[nt]=&All.InitDustTemp;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"DustToGasRatio"); 
      addr[nt]=&All.DustToGasRatio;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"AVConversionFactor"); 
      addr[nt]=&All.AVConversionFactor;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"CosmicRayIonRate"); 
      addr[nt]=&All.CosmicRayIonRate;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitRedshift"); 
      addr[nt]=&All.InitRedshift;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"ExternalDustExtinction"); 
      addr[nt]=&All.ExternalDustExtinction;
      id[nt++]=DOUBLE;

      strcpy(tag[nt], "PhotochemApprox");
      addr[nt] = &All.PhotochemApprox;
      id[nt++] = INT;

      strcpy(tag[nt], "ChemistryNetwork");
      addr[nt] = &All.ChemistryNetwork;
      id[nt++] = INT;

      strcpy(tag[nt], "ADRateFlag");
      addr[nt] = &All.ADRateFlag;
      id[nt++] = INT;

      strcpy(tag[nt], "MNRateFlag");
      addr[nt] = &All.MNRateFlag;
      id[nt++] = INT;

      strcpy(tag[nt], "AtomicFlag");
      addr[nt] = &All.AtomicFlag;
      id[nt++] = INT;

      strcpy(tag[nt], "ThreeBodyFlagA");
      addr[nt] = &All.ThreeBodyFlagA;
      id[nt++] = INT;

      strcpy(tag[nt], "ThreeBodyFlagB");
      addr[nt] = &All.ThreeBodyFlagB;
      id[nt++] = INT;

      strcpy(tag[nt], "H3PlusRateFlag");
      addr[nt] = &All.H3PlusRateFlag;
      id[nt++] = INT;

      strcpy(tag[nt],"InitMolHydroAbund"); 
      addr[nt]=&All.InitMolHydroAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitHPlusAbund"); 
      addr[nt]=&All.InitHPlusAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitDIIAbund"); 
      addr[nt]=&All.InitDIIAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitHDAbund"); 
      addr[nt]=&All.InitHDAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitHeIIAbund"); 
      addr[nt]=&All.InitHeIIAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitHeIIIAbund"); 
      addr[nt]=&All.InitHeIIIAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitCIIAbund"); 
      addr[nt]=&All.InitCIIAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitSiIIAbund"); 
      addr[nt]=&All.InitSiIIAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitOIIAbund"); 
      addr[nt]=&All.InitOIIAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitCOAbund"); 
      addr[nt]=&All.InitCOAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitC2Abund"); 
      addr[nt]=&All.InitC2Abund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitOHAbund"); 
      addr[nt]=&All.InitOHAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitH2OAbund"); 
      addr[nt]=&All.InitH2OAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitO2Abund"); 
      addr[nt]=&All.InitO2Abund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitHCOPlusAbund"); 
      addr[nt]=&All.InitHCOPlusAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitCHAbund"); 
      addr[nt]=&All.InitCHAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitCH2Abund"); 
      addr[nt]=&All.InitCH2Abund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitSiIIIAbund"); 
      addr[nt]=&All.InitSiIIIAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitCH3PlusAbund"); 
      addr[nt]=&All.InitCH3PlusAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitMgPlusAbund"); 
      addr[nt]=&All.InitMgPlusAbund;
      id[nt++]=DOUBLE;
#endif /* CHEMCOOL */

#ifdef POLYTROPE
      strcpy(tag[nt], "WhichEOS");
      addr[nt] = &All.WhichEOS;
      id[nt++] = INT;

      strcpy(tag[nt], "EOSFullTableSize");
      addr[nt] = &All.EOSFullTableSize;
      id[nt++] = INT;
#endif /* POLYTROPE */

      /*akj*/

      strcpy(tag[nt],"TurbulenceOn");
      addr[nt]=&All.TurbulenceOn;
      id[nt++]=INT;

#ifdef TURBULENCE

      strcpy(tag[nt],"DrvTimestep");
      addr[nt]=&All.DrvTimestep;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"VelDispConst");
      addr[nt]=&All.VelDispConst;
      id[nt++]=INT;
 
      strcpy(tag[nt],"LDrv");
      addr[nt]=&All.LDrv;
      id[nt++]=DOUBLE;
 
      strcpy(tag[nt],"MachNumber");
      addr[nt]=&All.MachNumber;
      id[nt++]=DOUBLE;
 
      strcpy(tag[nt],"StartDriving");
      addr[nt]=&All.StartDriving;
      id[nt++]=DOUBLE;
 
      strcpy(tag[nt],"DrvIndx");
      addr[nt]=&All.DrvIndx;
      id[nt++]=INT;
 
      strcpy(tag[nt],"kMin");
      addr[nt]=&All.kMin;
      id[nt++]=INT;
 
      strcpy(tag[nt],"kMax");
      addr[nt]=&All.kMax;
      id[nt++]=INT;
 
      strcpy(tag[nt],"Seed0");
      addr[nt]=&All.Seed0;
      id[nt++]=INT;
 
      strcpy(tag[nt],"Seed1");
      addr[nt]=&All.Seed1;
      id[nt++]=INT;
 
      strcpy(tag[nt],"Seed2");
      addr[nt]=&All.Seed2;
      id[nt++]=INT;

      strcpy(tag[nt],"MWeight");
      addr[nt]=&All.MWeight;
      id[nt++]=DOUBLE;
#endif 

       /* SINK: add parameters for sinks */
      strcpy(tag[nt], "HSinkCreate");
      addr[nt] = &All.HSinkCreate;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "RInner");
      addr[nt] = &All.RInner;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ROuter");
      addr[nt] = &All.ROuter;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SinkCriticalDens");
      addr[nt] = &All.SinkCriticalDens;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SinkCriticalRedshift");
      addr[nt] = &All.SinkCriticalRedshift;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MaxNumSinks");
      addr[nt] = &All.MaxNumSinks;
      id[nt++] = INT;

      strcpy(tag[nt], "RefinementMass");
      addr[nt] = &All.RefinementMass;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "max_dens");
      addr[nt] = &All.max_dens;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ray_crit_dens");
      addr[nt] = &All.ray_crit_dens;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ray_r_max_sink");
      addr[nt] = &All.ray_r_max_sink;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ray_flag_sun");
      addr[nt] = &All.ray_flag_sun;
      id[nt++] = INT;

      if((fd = fopen(fname, "r")))
	{
	  sprintf(buf, "%s%s", fname, "-usedvalues");
	  if(!(fdout = fopen(buf, "w")))
	    {
	      printf("error opening file '%s' \n", buf);
	      errorFlag = 1;
	    }
	  else
	    {
	      while(!feof(fd))
		{
		  *buf = 0;
		  fgets(buf, 200, fd);
		  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
		    continue;

		  if(buf1[0] == '%')
		    continue;

		  for(i = 0, j = -1; i < nt; i++)
		    if(strcmp(buf1, tag[i]) == 0)
		      {
			j = i;
			tag[i][0] = 0;
			break;
		      }

		  if(j >= 0)
		    {
		      switch (id[j])
			{
			case DOUBLE:
			  *((double *) addr[j]) = atof(buf2);
			  fprintf(fdout, "%-35s%g\n", buf1, *((double *) addr[j]));
			  break;
			case STRING:
			  strcpy(addr[j], buf2);
			  fprintf(fdout, "%-35s%s\n", buf1, buf2);
			  break;
			case INT:
			  *((int *) addr[j]) = atoi(buf2);
			  fprintf(fdout, "%-35s%d\n", buf1, *((int *) addr[j]));
			  break;
			}
		    }
		  else
		    {
		      fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
			      fname, buf1);
		      errorFlag = 1;
		    }
		}
	      fclose(fd);
	      fclose(fdout);

	      i = strlen(All.OutputDir);
	      if(i > 0)
		if(All.OutputDir[i - 1] != '/')
		  strcat(All.OutputDir, "/");

	      sprintf(buf1, "%s%s", fname, "-usedvalues");
	      sprintf(buf2, "%s%s", All.OutputDir, "parameters-usedvalues");
	      sprintf(buf3, "cp %s %s", buf1, buf2);
	      system(buf3);
	    }
	}
      else
	{
	  printf("\nParameter file %s not found.\n\n", fname);
	  errorFlag = 2;
	}

      if(errorFlag != 2)
	for(i = 0; i < nt; i++)
	  {
	    if(*tag[i])
	      {
		printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
		errorFlag = 1;
	      }
	  }

      if(All.OutputListOn && errorFlag == 0)
	errorFlag += read_outputlist(All.OutputListFilename);
      else
	All.OutputListLength = 0;

#ifdef XRAY_BACKGROUND
#ifdef XRAY_VARIABLE_HEATING
      if(errorFlag == 0)
	errorFlag += read_xrbIntensity(All.xrbFile);
#endif /* XRAY_VARIABLE_HEATING */
#endif /* XRAY_BACKGROUND */

#ifdef COSMIC_RAY_BACKGROUND
      initialize_cosmic_ray_background();
#ifdef CR_VARIABLE_HEATING
      if(errorFlag == 0)
	errorFlag += read_crbIntensity(All.crbFile);
#endif /* CR_VARIABLE_HEATING */
#endif /* COSMIC_RAY_BACKGROUND */
    }


  if(errorFlag)
    {
      exit(0);
    }



  if(All.NumFilesWrittenInParallel < 1)
    {
      if(ThisTask == 0)
	printf("NumFilesWrittenInParallel MUST be at least 1\n");
      exit(0);
    }

  if(All.NumFilesWrittenInParallel > NTask)
    {
      if(ThisTask == 0)
	printf("NumFilesWrittenInParallel MUST be smaller than number of processors\n");
      exit(0);
    }

#ifdef TURBULENCE
  if(All.TurbulenceOn==0)
    {
      if(ThisTask==0)
        {
          fprintf(stdout,"Code was compiled with turbulence switched on.\n");
          fprintf(stdout,"You must set `TurbulenceOn=1', or recompile the code.\n");
        }
          exit(0);
    }
#else
  if(All.TurbulenceOn==1)
    {
      if(ThisTask==0)
        {
          fprintf(stdout,"Code was compiled with turbulence switched off.\n");
          fprintf(stdout,"You must set `TurbulenceOn=0', or recompile the code.\n");
        }
      exit(0);
    }
#endif 


#ifdef PERIODIC
  if(All.PeriodicBoundariesOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with periodic boundary conditions switched on.\n");
	  printf("You must set `PeriodicBoundariesOn=1', or recompile the code.\n");
	}
      exit(0);
    }
#else
  if(All.PeriodicBoundariesOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with periodic boundary conditions switched off.\n");
	  printf("You must set `PeriodicBoundariesOn=0', or recompile the code.\n");
	}
      exit(0);
    }
#endif


  if(All.TypeOfTimestepCriterion >= 1)
    {
      if(ThisTask == 0)
	{
	  printf("The specified timestep criterion\n");
	  printf("is not valid\n");
	}
      exit(0);
    }

#if defined(LONG_X) ||  defined(LONG_Y) || defined(LONG_Z)
#ifndef NOGRAVITY
  if(ThisTask == 0)
    {
      printf("Code was compiled with LONG_X/Y/Z, but not with NOGRAVITY.\n");
      printf("Stretched periodic boxes are not implemented for gravity yet.\n");
    }
  exit(0);
#endif
#endif

#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS
}


/*! this function reads a table with a list of desired output times. The
 *  table does not have to be ordered in any way, but may not contain more
 *  than MAXLEN_OUTPUTLIST entries.
 */
int read_outputlist(char *fname)
{
  FILE *fd;

  if(!(fd = fopen(fname, "r")))
    {
      printf("can't read output list in file '%s'\n", fname);
      return 1;
    }

  All.OutputListLength = 0;
  do
    {
      if(fscanf(fd, " %lg ", &All.OutputListTimes[All.OutputListLength]) == 1)
	All.OutputListLength++;
      else
	break;
    }
  while(All.OutputListLength < MAXLEN_OUTPUTLIST);

  fclose(fd);

  printf("\nfound %d times in output-list.\n", All.OutputListLength);

  return 0;
}

#ifdef XRAY_BACKGROUND
#ifdef XRAY_VARIABLE_HEATING
/*! this function reads a table containing the average X-ray background intensity
 *  as a function of redshift. Table must be sorted in descending redshift order,
 *  and may not contain more than MAXLEN_HEATLIST entries.
 */
int read_xrbIntensity(char *fname)
{
  FILE *fd;

  if(!(fd = fopen(fname, "r")))
    {
      printf("can't read X-ray background intensity list in file '%s'\n", fname);
      return 1;
    }

  All.xrbLength = 0;
  do
    {
      if(fscanf(fd, "%lg %lg", &All.Jz[All.xrbLength], &All.Jxr[All.xrbLength]) == 2)
	All.xrbLength++;
      else
	break;
    }
  while(All.xrbLength < MAXLEN_HEATLIST);
  fclose(fd);

  printf("\nfound %d redshift points in X-ray background intensity list.\n", All.xrbLength);

  return 0;
}
#endif /* XRAY_VARIABLE_HEATING */
#endif /* XRAY_BACKGROUND */

#ifdef COSMIC_RAY_BACKGROUND
#ifdef CR_VARIABLE_HEATING
/*! this function reads a table containing the average Cosmic Ray background energy density
 *  as a function of redshift. Table must be sorted in descending redshift order,
 *  and may not contain more than MAXLEN_HEATLIST entries.
 */
int read_crbIntensity(char *fname)
{
  FILE *fd;

  if(!(fd = fopen(fname, "r")))
    {
      printf("can't read cosmic ray background intensity list in file '%s'\n", fname);
      return 1;
    }

  All.crbLength = 0;
  do
    {
      if(fscanf(fd, "%lg %lg", &All.U_CRz[All.crbLength], &All.U_CR[All.crbLength]) == 2)
	All.crbLength++;
      else
	break;
    }
  while(All.crbLength < MAXLEN_HEATLIST);
  fclose(fd);

  printf("\nfound %d redshift points in cosmic ray background intensity list.\n", All.crbLength);

  return 0;
}
#endif /* CR_VARIABLE_HEATING */
#endif /* COSMIC_RAY_BACKGROUND */
