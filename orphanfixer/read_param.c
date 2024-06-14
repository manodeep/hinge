#include "hinge.h"
#include "read_param.h"
#include "utils_read_params.h"
#include "utils.h"
#include "macros.h"
#include "set_cosmology.h"

void set_simulation_params(struct params_data *params)
{
    params->RedShift = REDSHIFT;
    params->Age = my_malloc(sizeof(*(params->Age)), PARAMS.MAX_SNAPSHOT_NUM + 1);
    for (int i = PARAMS.MIN_SNAPSHOT_NUM; i <= PARAMS.MAX_SNAPSHOT_NUM; i++)
    {
        params->Age[i] = get_age(REDSHIFT[i]); /* Age of the Universe in GYR
                                                  corresponding to the redshifts */
        fprintf(stderr, "redshift = %10.4g    age = %12.4g Gyr\n", REDSHIFT[i], params->Age[i]);
    }
}


void orphanfixer_fill_params(struct params_data *params, const int maxtags, void **addr, int *id,
                                  char (*tag)[MAXLEN], int *nt_out)
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
