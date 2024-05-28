#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#define __USE_LARGEFILE64
#include <sys/types.h>
#include <unistd.h>

#include "io_hinge_ascii.h"
#include "macros.h"
#include "progressbar.h"
#include "utils.h"

static int64 get_num_fofs(const char *catalogue_fname, char comment);

int64 get_num_fofs(const char *catalogue_fname, char comment)
{
    FILE *fp = my_fopen(catalogue_fname, "rt");
    char buf[BUFSIZ];
    int64 num_fofs = 0;
    while (fgets(buf, BUFSIZ, fp) != NULL)
    {
        if (buf[0] == comment)
            continue;
        int64_t haloid, hostid;
        int nread = sscanf(buf, "%" SCNd64 "%" SCNd64, &haloid, &hostid);
        if (nread != 2)
        {
            fprintf(stderr, "Error: Could not read two ids (haloid, hosthaloid) from line = '%s'\n", buf);
            exit(EXIT_FAILURE);
        }
        if (haloid == hostid)
            num_fofs++;
    }

    fclose(fp);
    return num_fofs;
}

int64 returnNhalo_hinge_ascii(const struct params_data *params, const int snapnum, const int fof_only)
{
#if !defined(LONGIDS) || !defined(BIGSIM)
    fprintf(stderr, "Both 'LONGIDS and 'BIGSIM' need to be defined loading the HINGE-ascii format\n");
    fprintf(stderr, "Please enable those two options in the `common.mk` file and then recompile + re-run\n");
    exit(EXIT_FAILURE);
#endif

    char catalogue_fname[MAXLEN];
    my_snprintf(catalogue_fname, MAXLEN, "%s/%s_halos_z%0.3f.txt", params->GROUP_DIR, params->GROUP_BASE,
                REDSHIFT[snapnum]);

    if (fof_only == 0)
    {
        return getnumlines(catalogue_fname, '#');
    }
    else
    {
        return get_num_fofs(catalogue_fname, '#');
    }
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

    my_snprintf(particles_fname, MAXLEN, "%s/%s_particles_z%0.3f.txt", params->GROUP_DIR, params->GROUP_BASE,
                REDSHIFT[snapnum]);
    my_snprintf(catalogue_fname, MAXLEN, "%s/%s_halos_z%0.3f.txt", params->GROUP_DIR, params->GROUP_BASE,
                REDSHIFT[snapnum]);

    const int64 numgroups = returnNhalo_hinge_ascii(params, snapnum, 0);
    fcat = my_fopen(catalogue_fname, "rt");
    fpart = my_fopen(particles_fname, "rt");

    off_t offset = 0;
    while (fgets(buf1, MAXBUFSIZE, fcat) != NULL)
    {
        if (buf1[0] == '#')
        {
            offset = ftell(fcat);
            continue;
        }
        else
            break;
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
        else
            break;
    }
    fseek(fpart, offset, SEEK_SET);

    int interrupted = 0;
    init_my_progressbar(numgroups, &interrupted);
    int64_t fof_hostnum = -1, fof_hostid = -1;
    int64 ihalo = 0;
    while (fgets(buf1, MAXBUFSIZE, fcat) != NULL)
    {
        // ## haloid  hosthaloid  nsub  mvir npart xc yc zc vxc vyc vzc
        int64_t haloid, hosthaloid, nsub, npart;
        double mvir, xc, yc, zc, vxc, vyc, vzc;
        const int nread = sscanf(buf1, "%" RD_FMT " %" RD_FMT " %" RD_FMT " %lf %" RD_FMT " %lf %lf %lf %lf %lf %lf",
                           &haloid, &hosthaloid, &nsub, &mvir, &npart, &xc, &yc, &zc, &vxc, &vyc, &vzc);
        XASSERT(nread == 11, "Error: Could only read %d elements instead of the expected 11 (file = %s)\n", nread, catalogue_fname);
        my_progressbar(ihalo, &interrupted);
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
        if (haloid == hosthaloid)
        {
            fof_hostnum = ihalo;
            fof_hostid = haloid;
        }
        // group[ihalo].hosthaloid = hosthaloid;
        group[ihalo].isFof = (haloid == hosthaloid) ? 1 : 0;
        group[ihalo].FOFHalo = fof_hostnum;
        group[ihalo].ContainerIndex = fof_hostnum;

        group[ihalo].ParentLevel = (group[ihalo].isFof == 1) ? 1 : -1; // subhalos don't have a parentlevel defined yet

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

        XASSERT(fof_hostnum != -1 && fof_hostid != -1,
                "Error: ihalo %" PRId64 " fof_hostnum = %" PRId64 " and fof_hostid = %" PRId64
                " must both be set (haloid = %" PRId64 ") \n",
                ihalo, fof_hostnum, fof_hostid, haloid);
       
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
            if (sscanf(buf1, "%" SCNd64 " %" SCNd64 " %" SCNd64 " %d %f %f %f", &fofid, &part_haloid,
                       &group[ihalo].id[i], &type, &group[ihalo].x[i], &group[ihalo].y[i], &group[ihalo].z[i]) != 7)
            {
                fprintf(stderr, "%s>: Error reading particle ids from file %s\n", __FUNCTION__, particles_fname);
                exit(EXIT_FAILURE);
            }
            XASSERT(part_haloid == haloid, "All particles must belong to the same halo. ihalo = %"PRId64" part_haloid (particles file) = %"PRId64" haloid (halo catalog) = %"PRId64"\n", ihalo, part_haloid, haloid);
            XASSERT(fofid == fof_hostid, "All particles must belong to the same FOF halo. ihalo = %"PRId64" (haloid = %"PRId64") fofid (particles file) = %"PRId64" fof_hostid (halo catalog) = %"PRId64"\n", ihalo, haloid, fofid, fof_hostid);
        }

        ihalo++;
    }
    finish_myprogressbar(&interrupted);
    fclose(fcat);
    fclose(fpart);
}
