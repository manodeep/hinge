#pragma once

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
#define RD_FMT SCNd64
#else
#define int64 int
#define STR_FMT "d"
#define RD_FMT "d"
#endif

// Do not defined IDs to be unsigned. Code assigns -1 to ids to mark particles
// for deletion later.
#ifdef LONGIDS
#define id64 int64_t
#define STR_ID_FMT PRId64
#define RD_ID_FMT SCNd64
#else
#define id64 int
#define STR_ID_FMT "d"
#define RD_ID_FMT "d"
#endif

/* officially work with groups with larger N
   smaller groups will be used in case a match
   is not found in the next snapshot.

*/
#define GROUPMINLEN (100)
#define ABS_PERIODIC(x, y) (fabs(periodic(x - y)))

#define MAXLEN (1000)
#define MAXLINESIZE (10000)

enum valid_group_formats
{
    /* The number of input group catalogs supported
       This consists of two parts, the first part
       dictates the catalog kind (i.e., what the bytes mean), while
       the second part dictates the actual format on disk (i.e.,
       how to read/cast the bytes from disk) */
    subfind_binary = 0,
    hinge_ascii = 1,
    hinge_binary = 2,
    num_group_formats
};

struct params_data
{
    int MIN_SNAPSHOT_NUM;
    int MAX_SNAPSHOT_NUM;
    /*   int SNAPSHOT_NUMBER; */

    char SNAPSHOT_DIR[MAXLEN];
    char SNAPSHOT_BASE[MAXLEN];

    char GROUP_DIR[MAXLEN];
    char GROUP_BASE[MAXLEN];
    enum valid_group_formats GROUP_FORMAT;

    char OUTPUT_DIR[MAXLEN];

    int MAX_INCR;
    int64 MAX_RANK_LOC;

    double MIN_FCOMMON_FINDPROGENITOR_THRESH; // findprogenitor
    int64 MIN_NUMPART_IN_FINDPROGENITOR_HALO; // findprogenitor

    double MIN_FCOMMON_SWITCHFOF_THRESH; // switchfof
    int64 MIN_NUMPART_IN_SWITCHFOF_HALO; // switchfof

    /* Populated from the GADGET snapshots (not read in from parameter file) */
    /* Can be set in the parameter file now: MS 13th June, 2024 */
    double BOXSIZE;
    double MASSARR[6]; // Gadget massarr

    /* Options for Orphanfixer */
    int MAX_DECR_GROUPS;
    int LOAD_FOUND_PROGENITORS;
    int LOAD_PARTIAL_FOUND_PROGENITORS;
    double MIN_FCOMMON_THRESH;

    // Makefile options
    int fof_only;
    int get_groupvel;
    int bigsim;
    int longids;
    int make_lean;
};

/* #ifdef SUSSING_TREES */
extern float *REDSHIFT;
extern struct params_data PARAMS;
/* #endif */
