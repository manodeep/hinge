#ifndef __DEFS_H
#define __DEFS_H

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

#ifdef BIGSIM
#define int64 int64_t
#define STR_FMT PRId64
#else
#define int64 int
#define STR_FMT "d"
#endif

// Do not defined IDs to be unsigned. Code assigns -1 to ids to mark particles
// for deletion later.
#ifdef LONGIDS
#define id64 int64_t
#define STR_ID PRId64
#else
#define id64 int
#define STR_ID "d"
#endif

/* officially work with groups with larger N
   smaller groups will be used in case a match
   is not found in the next snapshot.

*/
#define GROUPMINLEN (100)
#define ABS_PERIODIC(x, y) (fabs(periodic(x - y)))

#define MAXLEN (1000)
#define MAXLINESIZE (10000)

/* #ifdef SUSSING_TREES */
extern float *REDSHIFT;
/* #endif */

#endif
