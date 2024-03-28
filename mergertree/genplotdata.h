#pragma once

#include <stdio.h>
#include <stdlib.h>

#include "defs.h"
#include "maketree.h"

/* functions in genplotdata*/
void print_header_for_fma(FILE *fp);
void print_header_for_massloss(FILE *fp);
void output_plot_data(struct node_data *tree[], int64 *Ngroups);
void print_header_for_subsub_flyby(FILE *fp);
void print_header_for_fof_fof_flyby(FILE *fp);
void print_header_for_halfmass(FILE *fp);
int check_for_flyby(struct node_data *first, struct node_data *second,
                    float *rsep, float *vsep);
int find_fof_flyby_future(struct node_data *const n1,
                          struct node_data *const n2);
