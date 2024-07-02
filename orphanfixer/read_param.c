#include "read_param.h"
#include "hinge.h"
#include "macros.h"
#include "set_cosmology.h"
#include "utils.h"
#include "utils_read_params.h"

void orphanfixer_fill_params(struct params_data *params, const int maxtags, void **addr, int *id, char (*tag)[MAXLEN],
                             int *nt_out)
{
    int nt = *nt_out;
    my_snprintf(tag[nt], MAXLEN, "MAX_DECR_GROUPS");
    addr[nt] = &(params->MAX_DECR_GROUPS);
    id[nt++] = INT;
    XASSERT(nt < maxtags, "Error: parameter index =%d must be less than max number of parameters=%d\n", nt, maxtags);

    my_snprintf(tag[nt], MAXLEN, "MAX_RANK_LOC");
    addr[nt] = &(params->MAX_RANK_LOC);
    id[nt++] = INT64;
    XASSERT(nt < maxtags, "Error: parameter index =%d must be less than max number of parameters=%d\n", nt, maxtags);

    my_snprintf(tag[nt], MAXLEN, "MIN_FCOMMON_THRESH");
    addr[nt] = &(params->MIN_FCOMMON_THRESH);
    id[nt++] = DOUBLE;
    XASSERT(nt < maxtags, "Error: parameter index =%d must be less than max number of parameters=%d\n", nt, maxtags);

    my_snprintf(tag[nt], MAXLEN, "LOAD_FOUND_PROGENITORS");
    addr[nt] = &(params->LOAD_FOUND_PROGENITORS);
    id[nt++] = INT;
    XASSERT(nt < maxtags, "Error: parameter index =%d must be less than max number of parameters=%d\n", nt, maxtags);

    my_snprintf(tag[nt], MAXLEN, "LOAD_PARTIAL_FOUND_PROGENITORS");
    addr[nt] = &(params->LOAD_PARTIAL_FOUND_PROGENITORS);
    id[nt++] = INT;
    XASSERT(nt < maxtags, "Error: parameter index =%d must be less than max number of parameters=%d\n", nt, maxtags);

    *nt_out = nt;
}


void orphanfixer_write_params(FILE *fp, struct params_data *params)
{
    fprintf(fp, "MAX_DECR_GROUPS                   %d\n", params->MAX_DECR_GROUPS);
    fprintf(fp, "MAX_RANK_LOC                      %" STR_FMT "\n", params->MAX_RANK_LOC);

    fprintf(fp, "MIN_FCOMMON_THRESH                %lf\n", params->MIN_FCOMMON_THRESH);
    fprintf(fp, "LOAD_FOUND_PROGENITORS            %d\n", params->LOAD_FOUND_PROGENITORS);
    fprintf(fp, "LOAD_PARTIAL_FOUND_PROGENITORS    %d\n", params->LOAD_PARTIAL_FOUND_PROGENITORS);

    return;
}