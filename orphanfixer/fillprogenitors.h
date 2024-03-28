#pragma once

#include "defs.h"
#include "io.h"
#include <float.h>

#include "maketree.h" //for assign_haloid function definition and for node_data struct definition

/* functions in fillprogenitor */
short load_found_progenitors(struct node_data *tree[], int64 *Ngroups,
                             const char *fname, int64 *startgroupnum);
void fillprogenitors(struct node_data *tree[], int64 *Ngroups);
