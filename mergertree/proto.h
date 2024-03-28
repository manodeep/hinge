#ifndef __PROTO_H
#define __PROTO_H

#include <stdio.h>
#include <stdlib.h>

#include "io.h"
#include "loadsnapshot.h"
#include "maketree.h"

/* functions in main */
void print_makefile_options(void);
void readsubhalo_hierarchy_levels(const char *fname, struct group_data *group0);
void output_parents(struct node_data *tree[], int64 *Ngroups);
void print_flyby_future_header(FILE *fp);
void output_flyby_futures(struct node_data *tree[], int64 *Ngroups);
void output_all_haloids(struct node_data *tree[], int64 *Ngroups);
void output_children(struct node_data *tree[], int64 *Ngroups);
void output_mergers(struct node_data *tree[], int64 *Ngroups);
void output_mergers_for_haloid(struct node_data *tree[], int64 *Ngroups,
                               const int64 haloid, const float mineta);
void assign_node(struct group_data *group0, int64 Ngroups0,
                 struct parent_data *parent, struct node_data *node,
                 int isnapshot);
void writeids(struct node_data *tree[], int64 *Ngroups);
void assign_vxcm(struct group_data *group, int64 Ngroups,
                 struct particle_data *P, float redshift);
void cumulative_merger_history(struct node_data *tree[], int64 *Ngroups);
void output_interesting_halos(struct node_data *tree[], int64 *Ngroups);
void output_suspicious_halos(struct node_data *tree[], int64 *Ngroups);
void output_milkyway_halos(struct node_data *tree[], int64 *Ngroups);
void output_gill_data(struct node_data *tree[], int64 *Ngroups);
struct node_data *partial_walk_tree(struct node_data *this,
                                    struct node_data *start);
struct node_data *walk_tree(struct node_data *start);

#endif
