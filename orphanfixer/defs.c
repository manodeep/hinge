#include "defs.h"
#include "read_param.h"

struct params_data PARAMS;
float *REDSHIFT;
/* int Nbins = 100; /\*Bins in log r. Gets passed to numerical fitting routine
 * (where it might get modified for very small halos)*\/ */
int64 NUMPART;
int64 MaxHaloId = 0; /* defined in assign_haloid function in maketree.c */
