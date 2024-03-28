#ifndef __PROTO_H
#define __PROTO_H

#include "io.h"
#include "loadparents.h"
#include "loadsnapshot.h"
#include "maketree.h"

/* functions in main */
void output_mergers(struct node_data *tree[], int64 *Ngroups, const char *outpath);
void output_mergers_for_haloid(struct node_data *tree[], int64 *Ngroups, const char *outpath, const int64 haloid,
                               const float mineta);
void assign_node(struct group_data *group0, int64 Ngroups0, struct parent_data *parent, struct node_data *node,
                 int isnapshot);
void writeids(struct node_data *tree[], int64 *Ngroups, const char *outpath);
void assign_vxcm(struct group_data *group, int64 Ngroups, struct particle_data *P, float redshift);
void readsubhalo_hierarchy_levels(const char *fname, struct group_data *group0);
void output_data_for_SUSSING_TREES(struct node_data *tree[], int64 *Ngroups);
void print_makefile_options(void);

#endif
