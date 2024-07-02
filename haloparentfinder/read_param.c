#include "read_param.h"
#include "utils.h"
#include "utils_read_params.h"

void haloparentfinder_fill_params(struct params_data *params, const int maxtags, void **addr, int *id,
                                  char (*tag)[MAXLEN], int *nt_out)
{
    int nt = *nt_out;
    my_snprintf(tag[nt], MAXLEN, "MAX_INCR");
    addr[nt] = &(params->MAX_INCR);
    id[nt++] = INT;
    XASSERT(nt < maxtags, "Error: parameter index =%d must be less than max number of parameters=%d\n", nt, maxtags);

    my_snprintf(tag[nt], MAXLEN, "MAX_RANK_LOC");
    addr[nt] = &(params->MAX_RANK_LOC);
    id[nt++] = INT64;
    XASSERT(nt < maxtags, "Error: parameter index =%d must be less than max number of parameters=%d\n", nt, maxtags);

    my_snprintf(tag[nt], MAXLEN, "MIN_FCOMMON_FINDPROGENITOR_THRESH");
    addr[nt] = &(params->MIN_FCOMMON_FINDPROGENITOR_THRESH);
    id[nt++] = DOUBLE;
    XASSERT(nt < maxtags, "Error: parameter index =%d must be less than max number of parameters=%d\n", nt, maxtags);

    my_snprintf(tag[nt], MAXLEN, "MIN_NUMPART_IN_FINDPROGENITOR_HALO");
    addr[nt] = &(params->MIN_NUMPART_IN_FINDPROGENITOR_HALO);
    id[nt++] = INT64;
    XASSERT(nt < maxtags, "Error: parameter index =%d must be less than max number of parameters=%d\n", nt, maxtags);

    my_snprintf(tag[nt], MAXLEN, "MIN_FCOMMON_SWITCHFOF_THRESH");
    addr[nt] = &(params->MIN_FCOMMON_SWITCHFOF_THRESH);
    id[nt++] = DOUBLE;
    XASSERT(nt < maxtags, "Error: parameter index =%d must be less than max number of parameters=%d\n", nt, maxtags);

    my_snprintf(tag[nt], MAXLEN, "MIN_NUMPART_IN_SWITCHFOF_HALO");
    addr[nt] = &(params->MIN_NUMPART_IN_SWITCHFOF_HALO);
    id[nt++] = INT64;
    XASSERT(nt < maxtags, "Error: parameter index =%d must be less than max number of parameters=%d\n", nt, maxtags);

    *nt_out = nt;
}


void haloparentfinder_write_params(FILE *fp, struct params_data *params)
{
    fprintf(fp, "MAX_INCR                          %d\n", params->MAX_INCR);
    fprintf(fp, "MAX_RANK_LOC                      %" STR_FMT "\n", params->MAX_RANK_LOC);

    fprintf(fp, "MIN_FCOMMON_FINDPROGENITOR_THRESH                       %lf\n",  params->MIN_FCOMMON_FINDPROGENITOR_THRESH);
    fprintf(fp, "MIN_NUMPART_IN_FINDPROGENITOR_HALO                      %" STR_FMT "\n", params->MIN_NUMPART_IN_FINDPROGENITOR_HALO);

    fprintf(fp, "MIN_FCOMMON_SWITCHFOF_THRESH                           %lf\n", params->MIN_FCOMMON_SWITCHFOF_THRESH);
    fprintf(fp, "MIN_NUMPART_IN_SWITCHFOF_HALO                          %" STR_FMT "\n", params->MIN_NUMPART_IN_SWITCHFOF_HALO);

    return;
}
