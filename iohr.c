#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

#include "allvars.h"
#include "proto.h"



/*! \file io.c
 *  \brief Routines for reading a snapshot file from disk.
 */

static int n_type[numtype];
static long long ntot_type_all[numtype];

/*! This function tells the size of one data entry in each of the blocks
 *  defined for the output file. If one wants to add a new output-block, this
 *  function should be augmented accordingly.
 */
int get_bytes_per_blockelement(enum iofields blocknr)
{
  int bytes_per_blockelement = 0;

  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
      bytes_per_blockelement = 3 * sizeof(double);
      break;

    case IO_ID:
#ifdef LONGIDS
      bytes_per_blockelement = sizeof(long long);
#else
      bytes_per_blockelement = sizeof(int);
#endif
      break;

    case IO_MASS:
    case IO_U:
    case IO_RHO:
    case IO_HSML:
    case IO_POT:
    case IO_DTENTR:
    case IO_TSTP:
#ifdef POLYTROPE
    case IO_PRESSURE:
#endif
#ifdef CHEMCOOL
    case IO_GAMMA:
#endif
#ifdef METALS_TG
    case IO_METALLICITY:
#endif
#ifdef SINKVAL
    case IO_SINK:
#endif
      bytes_per_blockelement = sizeof(double);
      break;
#ifdef CHEMCOOL
    case IO_CHEM:
      bytes_per_blockelement = TRAC_NUM * sizeof(double);
      break;
#endif
#ifdef RAYTRACE
    case IO_COLUMN:
#ifdef CO_SHIELDING
      bytes_per_blockelement = 18 * sizeof(float);
#else
      bytes_per_blockelement = 12 * sizeof(float);
#endif
      break;
#endif
    }

  return bytes_per_blockelement;
}


/*! This function returns the type of the data contained in a given block of
 *  the output file. If one wants to add a new output-block, this function
 *  should be augmented accordingly.
 */
int get_datatype_in_block(enum iofields blocknr)
{
  int typekey;

  switch (blocknr)
    {
    case IO_ID:
#ifdef LONGIDS
      typekey = 2;		/* native long long */
#else
      typekey = 0;		/* native int */
#endif
      break;

    default:
      typekey = 1;		/* native float */
      break;
    }

  return typekey;
}


/*! This function informs about the number of elements stored per particle for
 *  the given block of the output file. If one wants to add a new
 *  output-block, this function should be augmented accordingly.
 */
int get_values_per_blockelement(enum iofields blocknr)
{
  int values = 0;

  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
      values = 3;
      break;

    case IO_ID:
    case IO_MASS:
    case IO_U:
    case IO_RHO:
    case IO_HSML:
    case IO_POT:
    case IO_DTENTR:
    case IO_TSTP:
#ifdef POLYTROPE
    case IO_PRESSURE:
#endif
#ifdef CHEMCOOL
    case IO_GAMMA:
#endif
#ifdef METALS_TG
    case IO_METALLICITY:
#endif
#ifdef SINKVAL
    case IO_SINK:
#endif
      values = 1;
      break;
#ifdef CHEMCOOL
    case IO_CHEM:
      values = TRAC_NUM;
      break;
#endif

#ifdef RAYTRACE
    case IO_COLUMN:
#ifdef CO_SHIELDING
      values = 18;
#else
      values = 12;
#endif
      break;
#endif
    }

  return values;
}


/*! This function determines how many particles there are in a given block,
 *  based on the information in the header-structure.  It also flags particle
 *  types that are present in the block in the typelist array. If one wants to
 *  add a new output-block, this function should be augmented accordingly.
 */
int get_particles_in_block(enum iofields blocknr, int *typelist)
{
  int i, nall, ntot_withmasses, ngas, nstars;

  nall = 0;
  ntot_withmasses = 0;

  for(i = 0; i < numtype; i++)
    {
      typelist[i] = 0;

      if(header.npart[i] > 0)
	{
	  nall += header.npart[i];
	  typelist[i] = 1;
	}

      if(All.MassTable[i] == 0)
	ntot_withmasses += header.npart[i];
    }

  ngas = header.npart[0];
  nstars = header.npart[4];


  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
    case IO_TSTP:
    case IO_ID:
    case IO_POT:

      return nall;
      break;

    case IO_MASS:
      for(i = 0; i < numtype; i++)
	{
	  typelist[i] = 0;
	  if(All.MassTable[i] == 0 && header.npart[i] > 0)
	    typelist[i] = 1;
	}
      return ntot_withmasses;
      break;

    case IO_U:
    case IO_RHO:
    case IO_HSML:
    case IO_DTENTR:
#ifdef POLYTROPE
    case IO_PRESSURE:
#endif
#ifdef CHEMCOOL
    case IO_CHEM:
    case IO_GAMMA:
#endif
#ifdef METALS_TG
    case IO_METALLICITY:
#endif
#ifdef SINKVAL
     case IO_SINK:
#endif
#ifdef RAYTRACE
    case IO_COLUMN:
#endif
      for(i = 1; i < numtype; i++)
	typelist[i] = 0;
      return ngas;
      break;
    }

  exit(212);
  return 0;
}



/*! This function tells whether or not a given block in the output file is
 *  present, depending on the type of simulation run and the compile-time
 *  options. If one wants to add a new output-block, this function should be
 *  augmented accordingly.
 */
int blockpresent(enum iofields blocknr)
{

#ifndef OUTPUTPOTENTIAL
  if(blocknr == IO_POT)
    return 0;
#endif

#ifndef OUTPUTACCELERATION
  if(blocknr == IO_ACCEL)
    return 0;
#endif

#ifndef OUTPUTCHANGEOFENTROPY
  if(blocknr == IO_DTENTR)
    return 0;
#endif

#ifndef OUTPUTTIMESTEP
  if(blocknr == IO_TSTP)
    return 0;
#endif

#ifdef POLYTROPE
#ifndef OUTPUTPRESSURE
  if(blocknr == IO_PRESSURE)
    return 0;
#endif
#endif

#ifdef RAYTRACE
#ifndef OUTPUTCOLUMN
  if (blocknr == IO_COLUMN)
    return 0;
#endif
#endif

  return 1;			/* default: present */
}




/*! This function associates a short 4-character block name with each block
 *  number.  This is stored in front of each block for snapshot
 *  FileFormat=2. If one wants to add a new output-block, this function should
 *  be augmented accordingly.
 */
void fill_Tab_IO_Labels(void)
{
  enum iofields i;

  for(i = 0; i < IO_NBLOCKS; i++)
    switch (i)
      {
      case IO_POS:
	strncpy(Tab_IO_Labels[IO_POS], "POS ", 4);
	break;
      case IO_VEL:
	strncpy(Tab_IO_Labels[IO_VEL], "VEL ", 4);
	break;
      case IO_ID:
	strncpy(Tab_IO_Labels[IO_ID], "ID  ", 4);
	break;
      case IO_MASS:
	strncpy(Tab_IO_Labels[IO_MASS], "MASS", 4);
	break;
      case IO_U:
	strncpy(Tab_IO_Labels[IO_U], "U   ", 4);
	break;
      case IO_RHO:
	strncpy(Tab_IO_Labels[IO_RHO], "RHO ", 4);
	break;
      case IO_HSML:
	strncpy(Tab_IO_Labels[IO_HSML], "HSML", 4);
	break;
      case IO_POT:
	strncpy(Tab_IO_Labels[IO_POT], "POT ", 4);
	break;
      case IO_ACCEL:
	strncpy(Tab_IO_Labels[IO_ACCEL], "ACCE", 4);
	break;
      case IO_DTENTR:
	strncpy(Tab_IO_Labels[IO_DTENTR], "ENDT", 4);
	break;
      case IO_TSTP:
	strncpy(Tab_IO_Labels[IO_TSTP], "TSTP", 4);
	break;
#ifdef POLYTROPE
      case IO_PRESSURE:
	strncpy(Tab_IO_Labels[IO_PRESSURE], "PRES", 4);
	break;
#endif
#ifdef CHEMCOOL
      case IO_CHEM:
	strncpy(Tab_IO_Labels[IO_CHEM], "CHEM", 4);
	break;
      case IO_GAMMA:
	strncpy(Tab_IO_Labels[IO_GAMMA], "GAMM", 4);
	break;
#endif
#ifdef RAYTRACE
      case IO_COLUMN:
	strncpy(Tab_IO_Labels[IO_COLUMN], "COLN", 4);
	break;
#endif
#ifdef METALS_TG
      case IO_METALLICITY:
	strncpy(Tab_IO_Labels[IO_METALLICITY], "METL", 4);
	break;
#endif
#ifdef SINKVAL
      case IO_SINK:
        strncpy(Tab_IO_Labels[IO_SINK], "SINK", 4);
        break;
#endif
      }
}

/*! This function returns a descriptive character string that describes the
 *  name of the block when the HDF5 file format is used.  If one wants to add
 *  a new output-block, this function should be augmented accordingly.
 *
 *  N.B. Names can be no more than 500 characters long! 
 */
void get_dataset_name(enum iofields blocknr, char *buf)
{

  strcpy(buf, "default");

  switch (blocknr)
    {
    case IO_POS:
      strcpy(buf, "Coordinates");
      break;
    case IO_VEL:
      strcpy(buf, "Velocities");
      break;
    case IO_ID:
      strcpy(buf, "ParticleIDs");
      break;
    case IO_MASS:
      strcpy(buf, "Masses");
      break;
    case IO_U:
      strcpy(buf, "InternalEnergy");
      break;
    case IO_RHO:
      strcpy(buf, "Density");
      break;
    case IO_HSML:
      strcpy(buf, "SmoothingLength");
      break;
    case IO_POT:
      strcpy(buf, "Potential");
      break;
    case IO_ACCEL:
      strcpy(buf, "Acceleration");
      break;
    case IO_DTENTR:
      strcpy(buf, "RateOfChangeOfEntropy");
      break;
    case IO_TSTP:
      strcpy(buf, "TimeStep");
      break;
#ifdef POLYTROPE
    case IO_PRESSURE:
      strcpy(buf, "Pressure");
      break;
#endif
#ifdef CHEMCOOL
    case IO_CHEM:
      strcpy(buf, "ChemicalAbundances");
      break;
    case IO_GAMMA:
      strcpy(buf, "Adiabatic index");
      break;
#endif
#ifdef RAYTRACE
    case IO_COLUMN:
      strcpy(buf, "Column densities");
      break;
#endif
#ifdef METALS_TG
    case IO_METALLICITY:
      strcpy(buf, "Metallicity");
      break;
#endif
#ifdef SINKVAL
    case IO_SINK:
      strcpy(buf, "SinkValue");
      break;
#endif
    }
}

/*! This catches I/O errors occuring for fread(). In this case we
 *  better stop.
 */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fread) on task=%d has occured: %s\n", ThisTask, strerror(errno));
      fflush(stdout);
      exit(778);
    }
  return nread;
}
