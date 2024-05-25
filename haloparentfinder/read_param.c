#include "read_param.h"
#include "defs.h"
#include "utils.h"
#include "utils_read_params.h"

void haloparentfinder_fill_params(struct params_data *params, const int maxtags, void **addr, int *id, char (*tag)[MAXLEN],
                                  int *nt_out)
{
    int nt = *nt_out;
    my_snprintf(tag[nt], MAXLEN, "MAX_INCR");
    addr[nt] = &(params->MAX_INCR);
    id[nt++] = INT;
    assert(nt < maxtags);

    my_snprintf(tag[nt], MAXLEN, "MAX_RANK_LOC");
    addr[nt] = &(params->MAX_RANK_LOC);
    id[nt++] = INT64;
    assert(nt < maxtags);

    my_snprintf(tag[nt], MAXLEN, "MIN_FCOMMON_FINDPROGENITOR_THRESH");
    addr[nt] = &(params->MIN_FCOMMON_FINDPROGENITOR_THRESH);
    id[nt++] = DOUBLE;
    assert(nt < maxtags);

    my_snprintf(tag[nt], MAXLEN, "MIN_NUMPART_IN_FINDPROGENITOR_HALO");
    addr[nt] = &(params->MIN_NUMPART_IN_FINDPROGENITOR_HALO);
    id[nt++] = INT64;
    assert(nt < maxtags);

    my_snprintf(tag[nt], MAXLEN, "MIN_FCOMMON_SWITCHFOF_THRESH");
    addr[nt] = &(params->MIN_FCOMMON_SWITCHFOF_THRESH);
    id[nt++] = DOUBLE;
    assert(nt < maxtags);

    my_snprintf(tag[nt], MAXLEN, "MIN_NUMPART_IN_SWITCHFOF_HALO");
    addr[nt] = &(params->MIN_NUMPART_IN_SWITCHFOF_HALO);
    id[nt++] = INT64;
    assert(nt < maxtags);

    *nt_out = nt;
}
