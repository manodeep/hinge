#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

#include "io_hinge_ascii.h"
#include "utils.h"

#ifndef LONGIDS
#error "LONGIDS must be defined"
#endif

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
#ifndef MAXBUFSIZE
#define MAXBUFSIZE 10000
#endif
    char buf1[MAXBUFSIZE];
    FILE *fcat = NULL, *fpart = NULL;

    my_snprintf(particles_fname, MAXLEN, "%s/%s_%03d_particles.txt", params->GROUP_DIR, params->GROUP_BASE, snapnum);
    my_snprintf(catalogue_fname, MAXLEN, "%s/%s_%03d_halocat.txt", params->GROUP_DIR, params->GROUP_BASE, snapnum);

    fcat = my_fopen(catalogue_fname, "rt");
    fpart = my_fopen(particles_fname, "rt");

    int64 ihalo = 0;
    off_t offset = 0;
    while (fgets(buf1, MAXBUFSIZE, fcat) != NULL)
    {
        if (buf1[0] == '#')
        {
            offset = ftell(fcat);
            continue;
        }
    }
    fseek(fcat, offset, SEEK_SET);

    offset = 0;
    while (fgets(buf1, MAXBUFSIZE, fpart) != NULL)
    {
        if (buf1[0] == '#')
        {
            offset = ftell(fpart);
            continue;
        }
    }
    fseek(fpart, offset, SEEK_SET);

    while (fgets(buf1, MAXBUFSIZE, fcat) != NULL)
    {
        // ## haloid  hosthaloid  nsub  mvir npart xc yc zc vxc vyc vzc
        int64 haloid, hosthaloid, nsub, npart;
        double mvir, xc, yc, zc, vxc, vyc, vzc;
        if (sscanf(buf1, "%" RD_FMT " %" RD_FMT " %" RD_FMT " %lf %" RD_FMT " %lf %lf %lf %lf %lf %lf", &haloid,
                   &hosthaloid, &nsub, &mvir, &npart, &xc, &yc, &zc, &vxc, &vyc, &vzc) != 11)
        {
            fprintf(stderr, "%s>: Error reading catalogue file %s\n", __FUNCTION__, catalogue_fname);
            exit(EXIT_FAILURE);
        }

        // group[i].N = SubLen[i];
        // group[i].nodeloc = i;
        // group[i].snapshot = num;
        // group[i].redshift = REDSHIFT[num];

        group[ihalo].Mtot = mvir;
        group[ihalo].groupnum = ihalo;
        group[ihalo].nodeloc = ihalo;
        group[ihalo].snapshot = snapnum;
        group[ihalo].redshift = REDSHIFT[snapnum];
        group[ihalo].N = npart;

        group[ihalo].haloID = haloid;
        int64 fof_hostnum = -1, fof_hostid = -1;
        if (haloid == hosthaloid)
        {
            fof_hostnum = ihalo;
            fof_hostid = haloid;
        }
        // group[ihalo].hosthaloid = hosthaloid;
        group[ihalo].isFof = (haloid == hosthaloid) ? 1 : 0;
        group[ihalo].FOFHalo = fof_hostnum;
        group[ihalo].ContainerIndex = fof_hostnum;

        group[ihalo].ParentLevel = (group[ihalo].isFof == 1) ? 1:-1; // subhalos don't have a parentlevel defined yet

        group[ihalo].xcen = xc;
        group[ihalo].ycen = yc;
        group[ihalo].zcen = zc;
        group[ihalo].vxcen = vxc;
        group[ihalo].vycen = vyc;
        group[ihalo].vzcen = vzc;

        group[ihalo].x = my_malloc(sizeof(group->x[0]), npart);
        group[ihalo].y = my_malloc(sizeof(group->y[0]), npart);
        group[ihalo].z = my_malloc(sizeof(group->z[0]), npart);
        group[ihalo].id = my_malloc(sizeof(group->id[0]), npart);

        assert(fof_hostnum != -1 && fof_hostid != -1 && "Both fofid and fofnum must be set before reading particles");

        group[ihalo].N_per_wedge = 0;
        /* initialise the parent finding variables*/
        group[ihalo].ParentId = -1;
        group[ihalo].NParents = 0;
        group[ihalo].Switched = 0;
        group[ihalo].ParentSnapshot = -1;
        group[ihalo].Ncommon = 0;
        group[ihalo].Rank = 0.0;
        group[ihalo].NpartinParent = 0;

        for (int64 i = 0; i < npart; i++)
        {
            int64_t fofid, part_haloid;
            int type;
            if (fgets(buf1, MAXBUFSIZE, fpart) == NULL)
            {
                fprintf(stderr, "%s>: Error reading particles file %s\n", __FUNCTION__, particles_fname);
                exit(EXIT_FAILURE);
            }
            if (sscanf(buf1, "%" RD_FMT " %" RD_FMT " %" RD_FMT " %d %f %f %f", &fofid, &part_haloid,
                       &group[ihalo].id[i], &type, &group[ihalo].x[i], &group[ihalo].y[i], &group[ihalo].z[i]) != 7)
            {
                fprintf(stderr, "%s>: Error reading particle ids from file %s\n", __FUNCTION__, particles_fname);
                exit(EXIT_FAILURE);
            }
            assert(fofid == fof_hostid && "All particles must belong to the same FOF halo");
            assert(part_haloid == haloid && "All particles must belong to the same halo");
        }

        ihalo++;
    }

    fclose(fcat);
    fclose(fpart);
}