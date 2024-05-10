#ifdef __cplusplus
extern "C"
{
#endif

#include "defs.h"
#include "io.h"

    extern int64 returnNhalo_subfind_binary(const struct params_data *params, const int snapnum, const int fof_only);
    extern void loadgroups_subfind_binary(int num, const struct params_data *params, struct group_data *group);

#ifdef __cplusplus
}
#endif