#pragma once

#include <stdio.h>
#include <stdlib.h>

#include "hinge.h"
#include "set_cosmology.h"

/* public functions in read_param.c*/
extern void mergertree_fill_params(struct params_data *params, const int maxtags, void **addr, int *id,
                                         char (*tag)[MAXLEN], int *nt_out);
