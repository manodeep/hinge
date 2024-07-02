#pragma once

#include "hinge.h"

extern void haloparentfinder_fill_params(struct params_data *params, const int maxtags, void **addr, int *id,
                                         char (*tag)[MAXLEN], int *nt_out);
extern void haloparentfinder_write_params(FILE *fp, struct params_data *params);
