#pragma once

#include <stdio.h>
#include <stdlib.h>

#include "defs.h"
#include "maketree.h"

/* functions in impulse.c*/
double get_impulse(struct node_data *g, struct node_data *f, float rsep,
                   float vsep);
double get_external_delpot_prim(struct node_data *g, struct node_data *f,
                                float rsep, float vsep);
double get_external_delpot_sec(struct node_data *g, struct node_data *f,
                               float rsep, float vsep);

double get_internal_delpot_prim(struct node_data *g, struct node_data *f,
                                float rsep, float vsep);
double get_internal_delpot_sec(struct node_data *g, struct node_data *f,
                               float rsep, float vsep);
