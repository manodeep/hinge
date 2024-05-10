#include <stdio.h>
#include <stdlib.h>

#include "io_hinge_ascii.h"
#include "utils.h"

int64 returnNhalo_hinge_ascii(const struct params_data *params, const int snapnum, const int fof_only)
{
    char catalogue_fname[MAXLEN];
    if (fof_only != 0)
    {
        fprintf(stderr, "%s>: fof_only is not supported for 'HINGE-ASCII' format\n", __FUNCTION__);
        return -1;
    }
    my_snprintf(catalogue_fname, MAXLEN, "%s/%s_%03d_halocat.txt", params->GROUP_DIR, params->GROUP_BASE, snapnum);
    return getnumlines(catalogue_fname, '#');
}

void loadgroups_hinge_ascii(const int snapnum, const struct params_data *params, struct group_data *group)
{
    char catalogue_fname[MAXLEN];
    char particles_fname[MAXLEN];
    FILE *fcat = NULL, *fpart = NULL;

    my_snprintf(particles_fname, MAXLEN, "%s/%s_%03d_particles.txt", params->GROUP_DIR, params->GROUP_BASE, snapnum);
    my_snprintf(catalogue_fname, MAXLEN, "%s/%s_%03d_halocat.txt", params->GROUP_DIR, params->GROUP_BASE, snapnum);

    fcat = my_fopen(catalogue_fname, "rt");
    fpart = my_fopen(particles_fname, "rt");

    while (1)
    {

        break;
    }

    fclose(fcat);
    fclose(fpart);
}