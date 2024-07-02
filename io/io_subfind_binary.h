#ifdef __cplusplus
extern "C"
{
#endif

#include "hinge.h"
#include "io.h"

    extern int64 returnNhalo_subfind_binary(const struct params_data *params, const int snapnum, const int fof_only);
    extern void loadgroups_subfind_binary(const struct params_data *params, int snapnum, struct group_data *group);

#ifdef __cplusplus
}
#endif