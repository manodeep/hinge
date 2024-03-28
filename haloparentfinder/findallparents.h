#pragma once

#include "io.h"

/* public functions in findallparent */
int64 findallparents(struct group_data *prevgroup, int64 PrevNsub,
                     struct group_data *nextgroup, int64 NextNsub,
                     const int snapshot, const char *outpath);
int64 findfofparents(struct group_data *prevgroup, int64 PrevNsub,
                     struct group_data *nextgroup, int64 NextNsub,
                     const char *outpath);
