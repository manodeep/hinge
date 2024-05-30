#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#define __USE_LARGEFILE64
#include <sys/types.h>
#include <unistd.h>

#include "io_hinge_binary.h"
#include "macros.h"
#include "progressbar.h"
#include "utils.h"

// static int64_t get_num_fofs(const char *catalogue_fname, char comment);

// int64_t get_num_fofs(const char *catalogue_fname, char comment)
// {
//     FILE *fp = my_fopen(catalogue_fname, "r");
//     fseek(fp, 8, SEEK_SET);//skip 'totnhalos'
//     int64_t num_fofs;
//     my_fread(&num_fofs, sizeof(num_fofs), 1, fp);
//     fclose(fp);
//     return num_fofs;
// }

int64 returnNhalo_hinge_binary(const struct params_data *params, const int snapnum, const int fof_only)
{
#if !defined(LONGIDS) || !defined(BIGSIM)
    fprintf(stderr, "Both 'LONGIDS and 'BIGSIM' need to be defined loading the HINGE-binary format\n");
    fprintf(stderr, "Please enable those two options in the `common.mk` file and then recompile + re-run\n");
    exit(EXIT_FAILURE);
#endif

    char catalogue_fname[MAXLEN];
    my_snprintf(catalogue_fname, MAXLEN, "%s/%s_halos_z%0.3f.bin", params->GROUP_DIR, params->GROUP_BASE,
                REDSHIFT[snapnum]);
    FILE *fp = my_fopen(catalogue_fname, "r");
    int64_t nhalos;
    if(fof_only == 0) {
        my_fread(&nhalos, sizeof(nhalos), 1, fp);
    } else {
        fseek(fp, sizeof(nhalos), SEEK_SET);//skip 'totnhalos'
        my_fread(&nhalos, sizeof(nhalos), 1, fp);//read 'totnfofs'
    }
    fclose(fp);
    return nhalos;
}

void loadgroups_hinge_binary(const int snapnum, const struct params_data *params, struct group_data *group)
{
    char catalogue_fname[MAXLEN];
    my_snprintf(catalogue_fname, MAXLEN, "%s/%s_halos_z%0.3f.bin",
                params->GROUP_DIR, params->GROUP_BASE, REDSHIFT[snapnum]);

    const int64_t numgroups = returnNhalo_hinge_binary(params, snapnum, 0);
    FILE *fcat = my_fopen(catalogue_fname, "r");
    int64_t totnhalos, totnfofs, totnpart;
    my_fread(&totnhalos, sizeof(totnhalos), 1, fcat);
    my_fread(&totnfofs, sizeof(totnfofs), 1, fcat);
    my_fread(&totnpart, sizeof(totnpart), 1, fcat);
    void *buf;
    int64_t numbuf = totnpart > totnhalos ? totnpart:totnhalos;
    buf = malloc(numbuf*sizeof(double));
#define READ_AND_ASSIGN_TO_GROUPS(field, src_type, dst_type) {      \
        my_fread(buf, sizeof(src_type), totnhalos, fcat);           \
        for(int64_t jj=0;jj<totnhalos;jj++) {                       \
            group[jj].field = (dst_type) ((src_type *) buf)[jj];    \
        }                                                           \
    }
    fprintf(stderr,"Reading halo catalog ...\n");
    READ_AND_ASSIGN_TO_GROUPS(haloID, int64_t, long);
    READ_AND_ASSIGN_TO_GROUPS(fofID, int64_t, long);
    READ_AND_ASSIGN_TO_GROUPS(Nsub, int64_t, int64);
    READ_AND_ASSIGN_TO_GROUPS(Mtot, double, float);
    READ_AND_ASSIGN_TO_GROUPS(N, int64_t, int64_t);
    READ_AND_ASSIGN_TO_GROUPS(xcen, double, float);
    READ_AND_ASSIGN_TO_GROUPS(ycen, double, float);
    READ_AND_ASSIGN_TO_GROUPS(zcen, double, float);
    READ_AND_ASSIGN_TO_GROUPS(vxcen, double, float);
    READ_AND_ASSIGN_TO_GROUPS(vycen, double, float);
    READ_AND_ASSIGN_TO_GROUPS(vzcen, double, float);
    fprintf(stderr,"Reading halo catalog ...done\n");

    int interrupted = 0;
    init_my_progressbar(numgroups, &interrupted);
#undef READ_AND_ASSIGN_TO_GROUPS
    for(int64_t ihalo=0;ihalo<numgroups;ihalo++) {
        my_progressbar(ihalo, &interrupted);
        group[ihalo].groupnum = ihalo;
        group[ihalo].nodeloc = ihalo;
        group[ihalo].snapshot = snapnum;
        group[ihalo].redshift = REDSHIFT[snapnum];

        int64_t haloid = group[ihalo].haloID;
        int64_t hosthaloid = group[ihalo].fofID;
        int64 fof_hostnum = -1; //, fof_hostid = -1;
        if (haloid == hosthaloid)
        {
            fof_hostnum = ihalo;
            // fof_hostid = haloid;
        }
        group[ihalo].isFof = (haloid == hosthaloid) ? 1 : 0;
        group[ihalo].FOFHalo = fof_hostnum;
        group[ihalo].ContainerIndex = fof_hostnum;

        group[ihalo].ParentLevel = (group[ihalo].isFof == 1) ? 1 : -1; // subhalos don't have a parentlevel defined yet

        group[ihalo].x = my_malloc(sizeof(group->x[0]), group[ihalo].N);
        group[ihalo].y = my_malloc(sizeof(group->y[0]), group[ihalo].N);
        group[ihalo].z = my_malloc(sizeof(group->z[0]), group[ihalo].N);
        group[ihalo].id = my_malloc(sizeof(group->id[0]), group[ihalo].N);

        // my_fread(group[ihalo].id, sizeof(id64), group[ihalo].N, fcat);
#define READ_AND_ASSIGN_TO_PARTICLES(field, ii, src_type, dst_type, npart) {    \
            my_fread(buf, sizeof(src_type), npart, fcat);                       \
            for (int64_t jj=0;jj<npart;jj++) {                                  \
                group[ii].field[jj] = (dst_type) ((src_type *) buf)[jj];        \
            }                                                                   \
        }
        READ_AND_ASSIGN_TO_PARTICLES(id, ihalo, int64_t, id64, group[ihalo].N)
        fseek(fcat, sizeof(int64_t)*group[ihalo].N, SEEK_CUR);//skip over the particle types
        READ_AND_ASSIGN_TO_PARTICLES(x, ihalo, double, float, group[ihalo].N);
        READ_AND_ASSIGN_TO_PARTICLES(y, ihalo, double, float, group[ihalo].N);
        READ_AND_ASSIGN_TO_PARTICLES(z, ihalo, double, float, group[ihalo].N);

#undef READ_AND_ASSIGN_TO_PARTICLES

        group[ihalo].N_per_wedge = 0;
        /* initialise the parent finding variables*/
        group[ihalo].ParentId = -1;
        group[ihalo].NParents = 0;
        group[ihalo].Switched = 0;
        group[ihalo].ParentSnapshot = -1;
        group[ihalo].Ncommon = 0;
        group[ihalo].Rank = 0.0;
        group[ihalo].NpartinParent = 0;

    }
    fclose(fcat);
    finish_myprogressbar(&interrupted);

}
