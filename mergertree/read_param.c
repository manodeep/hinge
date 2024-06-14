#include "read_param.h"
#include "hinge.h"
#include "macros.h"
#include "utils.h"
#include "utils_read_params.h"

void mergertree_fill_params(struct params_data *params, const int maxtags, void **addr, int *id, char (*tag)[MAXLEN],
                            int *nt_out)
{
    int nt = *nt_out;

    my_snprintf(tag[nt], MAXLEN, "LINKLENGTH");
    addr[nt] = &(params->LINKLENGTH);
    id[nt++] = DOUBLE;
    XASSERT(nt < maxtags, "Error: parameter index =%d must be less than max number of parameters=%d\n", nt, maxtags);

    my_snprintf(tag[nt], MAXLEN, "INDIVIDUAL_MERGERS");
    addr[nt] = &(params->INDIVIDUAL_MERGERS);
    id[nt++] = INT;
    XASSERT(nt < maxtags, "Error: parameter index =%d must be less than max number of parameters=%d\n", nt, maxtags);

    *nt_out = nt;
}
