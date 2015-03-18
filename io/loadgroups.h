#pragma once

#include "io.h"


/* functions in loadgroups */
void reorder_groups_on_array(const int64 Ngroups,struct group_data *group);
void remove_duplicate_particles(const int64 Ngroups, struct group_data *groups);
int64 returnNhalo(const char* buf);
int64 returnNhalo_SUSSING(const char* buf,const int fof_only);
void loadgroups(int num,struct group_data* group);
void assign_parentlevel(struct group_data *this_group, const int64 Ngroups, struct group_data *groups);
void remove_particles_from_all_groups(struct group_data *groups, const int Ngroups);
