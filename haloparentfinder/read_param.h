#pragma once

#include <stdio.h>

#include "defs.h"
// global variable
// extern struct params_data PARAMS;

#define DOUBLE 1
#define STRING 2
#define INT 3
#define INT64 4
#define MAXTAGS 20

extern void read_params(const char *fname, struct params_data *params, void special_params(struct params_data *params, const int maxtags, void **addr, int *id, char **tag, int *nt));
extern void output_params(const char *fname, struct params_data *params);
extern void fill_config_params(struct params_data *params);
extern void sanity_check_params(struct params_data *params);
extern void haloparentfinder_fill_params(struct params_data *params, const int maxtags, void **addr, int *id, char **tag, int *nt_out);
