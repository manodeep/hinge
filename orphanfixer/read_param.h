#pragma once

#include "hinge.h"

extern void orphanfixer_fill_params(struct params_data *params, const int maxtags, void **addr, int *id,
                                    char (*tag)[MAXLEN], int *nt_out);
extern void orphanfixer_write_params(FILE *fp, struct params_data *params);
