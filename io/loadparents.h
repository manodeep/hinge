#pragma once

#include <stdio.h>
#include <stdlib.h>

#include "hinge.h"

struct parent_data
{
    /*   int Ngroups; */
    short snapshot;
    int64 groupid;
    short isFof;
    int64 nsub;
    int64 npartinhalo;
    int64 parentid;
    short parentsnapshot;
    int64 npartinparent;
    short parentlevel;
    int64 containerid;
    double rank;
    int64 ncommon;
};

/* functions in loadparents*/
struct parent_data *loadparents(const char *fname, struct parent_data *parent, int64 Ngroups);
