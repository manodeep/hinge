#pragma once

#include <stdio.h>
#include <stdlib.h>

#include "defs.h"
#include "set_cosmology.h"

struct params_data
{ 
  int MIN_SNAPSHOT_NUM;
  int MAX_SNAPSHOT_NUM;

  char SNAPSHOT_DIR[MAXLEN];
  char SNAPSHOT_BASE[MAXLEN];

  char GROUP_DIR[MAXLEN];
  char GROUP_BASE[MAXLEN];

  char OUTPUT_DIR[MAXLEN];

  double LINKLENGTH;// linking length for the FOF algorithm
  int INDIVIDUAL_MERGERS;//whether or not merger histories for individual haloids are wanted. 

  /* Populated from the snapshots (not read in from parameter file) */
  double BOXSIZE;
  double MASSARR[6];//Gadget massarr

  //Makefile options
  int fof_only;
  int get_groupvel;
  int get_meanvel;
  int bigsim;
  int longids;
  int use_most_bound_for_centre;

  /* Simulation Parameters*/
  float Softening;
  float PhysicalLinkLength;
  float BoxSize;
  float *RedShift;
  double *Age;
  double *MassArr;

  /*Cosmology struct*/
  struct cosmology_data *COSMO;
};

//global variable
extern struct params_data PARAMS;

/* public functions in read_param.c*/
void set_simulation_params(struct params_data *params);
void read_params(const char *fname,struct params_data *params);
void output_params(const char *fname,struct params_data *params);
void fill_config_params(struct params_data *params);
void sanity_check_params(struct params_data *params);

