#ifndef __DEFS_H
#define __DEFS_H

#include <assert.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if __GNUC__
#if __x86_64__ || __ppc64__
#define ENVIRONMENT64
#else
#define ENVIRONMENT32
#endif
#endif

// if 32 bit, enable large file macros
#ifdef ENVIRONMENT32
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
#endif

#define SQR(X) ((X) * (X))
#define CUBE(X) ((X) * (X) * (X))
#define PI (4.0 * atan(1.0))
#define DOUBLE_EPS (1e-10)

#ifdef BIGSIM
#define int64 int64_t
#define STR_FMT PRId64
#else
#define int64 int
#define STR_FMT "d"
#endif

#ifdef LONGIDS
#define id64 int64_t
#define STR_ID_FMT PRId64
#else
#define id64 int
#define STR_ID_FMT "d"
#endif

// filenames will be assumed to be at most MAXLEN characters long
#define MAXLEN (1000)
#define MAXLINESIZE (10000)

extern float *REDSHIFT;
extern int Nbins;
extern int64 NUMPART;
extern int64 MaxHaloId; /* defined in assign_haloid function in maketree.c */

#endif
