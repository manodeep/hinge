#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

#include "io.h"

    extern int64 returnNhalo_hinge_binary(const struct params_data *params, const int snapnum, const int fof_only);
    extern void loadgroups_hinge_binary(const struct params_data *params, const int snapnum, struct group_data *group);

#ifdef __cplusplus
}
#endif