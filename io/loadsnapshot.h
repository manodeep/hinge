#pragma once

#include "hinge.h"
#include "io.h"

struct particle_data
{
    float Pos[3];
    float Vel[3];
    float Mass;
    int32_t Type;
};

struct particle_data *loadsnapshot(const char *fname, struct io_header *header);
void reordering(struct particle_data *P, id64 *Id, int64 N);
