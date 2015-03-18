#pragma once

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

  int MAX_DECR_GROUPS;//fillprogenitor will load (MAX_DECR_GROUPS+2) groups files. Earliest "fixed" match will be MAX_DECR_GROUPS+1 snapshots away
  int64 MAX_RANK_LOC;
  double MIN_FCOMMON_THRESH;
  int LOAD_FOUND_PROGENITORS;
  int LOAD_PARTIAL_FOUND_PROGENITORS;

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
