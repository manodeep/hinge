#pragma once

#include <stdio.h>

#include "defs.h"

extern void haloparentfinder_fill_params(struct params_data *params, const int maxtags, void **addr, int *id,
                                         char (*tag)[MAXLEN], int *nt_out);
