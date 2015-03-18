#pragma once

#include <stdio.h>
#include <stdlib.h>

#include "defs.h"

struct params_data
{
  int MIN_SNAPSHOT_NUM;
  int MAX_SNAPSHOT_NUM;
/*   int SNAPSHOT_NUMBER; */
  
  char SNAPSHOT_DIR[MAXLEN];
  char SNAPSHOT_BASE[MAXLEN];
  
  char GROUP_DIR[MAXLEN];
  char GROUP_BASE[MAXLEN];

  char OUTPUT_DIR[MAXLEN];
  
  int MAX_INCR;
  int64 MAX_RANK_LOC;

  double MIN_FCOMMON_FINDPROGENITOR_THRESH;//findprogenitor
  int64  MIN_NUMPART_IN_FINDPROGENITOR_HALO;//findprogenitor

  double MIN_FCOMMON_SWITCHFOF_THRESH;//switchfof
  int64 MIN_NUMPART_IN_SWITCHFOF_HALO;//switchfof

  /* Populated from the GADGET snapshots (not read in from parameter file) */
  double BOXSIZE;
  double MASSARR[6];//Gadget massarr

  //Makefile options
  int fof_only;
  int get_groupvel;
  int bigsim;
  int longids;
  int make_lean;
};

//global variable
extern struct params_data PARAMS;


void read_params(const char *fname,struct params_data *params);
void output_params(const char *fname,struct params_data *params);
void fill_config_params(struct params_data *params);
void sanity_check_params(struct params_data *params);
