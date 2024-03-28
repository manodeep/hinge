#pragma once

#include "io.h"

/* public functions in findprogenitor.c */
void find_progenitor(struct group_data *nextgroup, int64 NextNsub,
                     struct group_data *prevgroup, int64 PrevNsub,
                     const char *outpath);
