#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

#include "defs.h"
#include "io.h"

    extern int64 returnNhalo_hinge_ascii(const struct params_data *params, const int snapnum, const int fof_only);
    extern void loadgroups_hinge_ascii(const int snapnum, const struct params_data *params, struct group_data *group);

#ifdef __cplusplus
}
#endif