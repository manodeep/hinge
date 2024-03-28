#pragma once

#include "io.h"

void check_fof_matches(struct group_data *prevgroup, int64 PrevNsub, struct group_data *nextgroup, int64 NextNsub,
                       const int snapshot, const char *outpath);
