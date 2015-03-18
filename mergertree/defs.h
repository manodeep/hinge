#ifndef __DEFS_H
#define __DEFS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <stdint.h>
#include <stdarg.h>
#include <inttypes.h>

#define   SQR(X)       ((X)*(X))
#define   CUBE(X)       ((X)*(X)*(X))
#define   PI           (4.0*atan(1.0))
#define   DOUBLE_EPS    (1e-10)

#ifdef BIGSIM
#define int64 int64_t
#define STR_FMT  PRId64
#else
#define int64 int
#define STR_FMT  "d"
#endif


#ifdef LONGIDS
#define id64 int64_t
#define STR_ID_FMT   PRId64
#else
#define id64  int
#define STR_ID_FMT "d"
#endif

//max. length for filenames (and small read buffers)
#define MAXLEN             (1000)
//max. length for read buffers
#define MAXLINESIZE  (10000)


extern double ActualMassUnits;
extern float *REDSHIFT;
extern int Nbins;
extern int64 NUMPART;
extern int64 MaxHaloId; /* defined in assign_haloid function in maketree.c */

#endif
