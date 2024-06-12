#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
// #define __USE_LARGEFILE64
#include <sys/types.h>
#include <unistd.h>

#include "io_hinge_ascii.h"
#include "io_hinge_binary.h"
#include "io_hinge_utils.h"
#include "macros.h"
#include "progressbar.h"
#include "utils.h"

int64 returnNhalo_hinge_binary(const struct params_data *params, const int snapnum, const int fof_only)
{
    return returnNhalo_hinge_ascii(params, snapnum, fof_only);
}

void loadgroups_hinge_binary(const int snapnum, const struct params_data *params, struct group_data *group)
{
    XASSERT(group != NULL, "group is NULL\n");

    char catalogue_fname[MAXLEN];
    my_snprintf(catalogue_fname, MAXLEN, "%s/%s_halos_z%0.3f.txt", params->GROUP_DIR, params->GROUP_BASE,
                REDSHIFT[snapnum]);
#ifdef FOF_ONLY
    const int fof_only = 1;
#else
    const int fof_only = 0;
#endif

    struct hinge_catalog *halocat = read_hinge_ascii_halo_catalog(catalogue_fname, fof_only);
    if (halocat == NULL)
    {
        fprintf(stderr, "Could not read the halo catalog from file %s\n", catalogue_fname);
        exit(EXIT_FAILURE);
    }
    const int64_t nhalos = halocat->nhalos;
    const int64_t totnpart = halocat->totnpart;

    /* read individual files for each column */

    /* Can't really automate the process - so need to read in individually and assign */

    // const char field_names[][MAXLEN] = {"partid", "xpos", "ypos", "zpos", "haloid", "fofid"};
#define CHECK_NPART_AND_READ_FIELD(field_name, totnpart, buf)                                                          \
    {                                                                                                                  \
        fprintf(stderr, "Reading field = '%s' ...\n", field_name);                                                     \
        int64_t npart_field;                                                                                           \
        char field_fname[MAXLEN];                                                                                      \
        my_snprintf(field_fname, MAXLEN, "%s/%s_particles_z%0.3f_%s.bin", params->GROUP_DIR, params->GROUP_BASE,       \
                    REDSHIFT[snapnum], field_name);                                                                    \
        FILE *fp = my_fopen(field_fname, "r");                                                                         \
        my_fread(&npart_field, sizeof(npart_field), 1, fp);                                                            \
        if (npart_field != totnpart)                                                                                   \
        {                                                                                                              \
            fprintf(stderr,                                                                                            \
                    "Number of particles in file %s = %" PRId64 " does not match total number of particles = %" PRId64 \
                    "\n",                                                                                              \
                    field_fname, npart_field, totnpart);                                                               \
            fclose(fp);                                                                                                \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
        const size_t sizeof_every_field = 8;                                                                           \
        buf = my_malloc(sizeof_every_field, totnpart);                                                                 \
        my_fread(buf, sizeof_every_field, totnpart,                                                                    \
                 fp); /* all the fields are of type int64_t or double (i.e., 8 bytes) */                               \
        fclose(fp);                                                                                                    \
        fprintf(stderr, "Reading field = '%s' ...done\n", field_name);                                                 \
    }

#define ASSIGN_FIELD_TO_GROUPS(field_name, field_type, nhalos, buf, dst_field)                                         \
    {                                                                                                                  \
        int64_t offset = 0;                                                                                            \
        int interrupted = 0;                                                                                           \
        fprintf(stderr, "Assigning field '%s' to groups ...\n", field_name);                                           \
        init_my_progressbar(nhalos, &interrupted);                                                                     \
        for (int64_t i = 0; i < nhalos; i++)                                                                           \
        {                                                                                                              \
            my_progressbar(i, &interrupted);                                                                           \
            const int64_t npart_field = group[i].N;                                                                    \
            group[i].dst_field = my_malloc(sizeof(field_type), npart_field);                                           \
            for (int64_t j = 0; j < npart_field; j++)                                                                  \
            {                                                                                                          \
                group[i].dst_field[j] = ((field_type *)buf)[offset + j];                                               \
            }                                                                                                          \
            offset += npart_field;                                                                                     \
        }                                                                                                              \
        free(buf);                                                                                                     \
        finish_myprogressbar(&interrupted);                                                                            \
        fprintf(stderr, "Assigning field '%s' to groups ...done\n", field_name);                                       \
    }
    void *buf;
    CHECK_NPART_AND_READ_FIELD("partid", totnpart, buf);
    ASSIGN_FIELD_TO_GROUPS("partid", int64_t, nhalos, buf, id);

    //    CHECK_NPART_AND_READ_FIELD(fp, "part_type", totnpart, buf);
    CHECK_NPART_AND_READ_FIELD("xpos", totnpart, buf);
    ASSIGN_FIELD_TO_GROUPS("xpos", double, nhalos, buf, x);

    CHECK_NPART_AND_READ_FIELD("ypos", totnpart, buf);
    ASSIGN_FIELD_TO_GROUPS("ypos", double, nhalos, buf, y);

    CHECK_NPART_AND_READ_FIELD("zpos", totnpart, buf);
    ASSIGN_FIELD_TO_GROUPS("zpos", double, nhalos, buf, z);

    int64_t *haloids, *fofids;
    CHECK_NPART_AND_READ_FIELD("haloid", totnpart, haloids);
    CHECK_NPART_AND_READ_FIELD("fofid", totnpart, fofids);

    int64_t offset = 0;
    int interrupted = 0;
    init_my_progressbar(nhalos, &interrupted);
    for (int64_t ihalo = 0; ihalo < nhalos; ihalo++)
    {
        my_progressbar(ihalo, &interrupted);
        const int64_t *fids = &fofids[offset];
        const int64_t *hids = &haloids[offset];
        group[ihalo].groupnum = ihalo;
        group[ihalo].nodeloc = ihalo;
        group[ihalo].snapshot = snapnum;
        group[ihalo].redshift = REDSHIFT[snapnum];

        int64_t haloid = group[ihalo].haloID;
        int64_t hosthaloid = group[ihalo].fofID;
        int64_t fof_hostnum = -1;
        if (haloid == hosthaloid)
        {
            fof_hostnum = ihalo;
            // fof_hostid = haloid;
        }
        group[ihalo].isFof = (haloid == hosthaloid) ? 1 : 0;
        group[ihalo].FOFHalo = fof_hostnum;
        group[ihalo].ContainerIndex = fof_hostnum;
        for (int64_t j = 0; j < group[ihalo].N; j++)
        {
            XASSERT(hids[j] == haloid, "Haloid mismatch: %" PRId64 " != %" PRId64, hids[j], haloid);
            XASSERT(fids[j] == hosthaloid, "Fofid mismatch: %" PRId64 " != %" PRId64, fids[j], hosthaloid);
        }

        group[ihalo].ParentLevel = (group[ihalo].isFof == 1) ? 1 : -1; // subhalos don't have a parentlevel defined yet
        group[ihalo].N_per_wedge = 0;
        /* initialise the parent finding variables*/
        group[ihalo].ParentId = -1;
        group[ihalo].NParents = 0;
        group[ihalo].Switched = 0;
        group[ihalo].ParentSnapshot = -1;
        group[ihalo].Ncommon = 0;
        group[ihalo].Rank = 0.0;
        group[ihalo].NpartinParent = 0;

        offset += group[ihalo].N;
    }
    finish_myprogressbar(&interrupted);

    fprintf(stderr, "Removing duplicates ...\n");
    const int64_t num_removed = remove_duplicates(group, nhalos);
    fprintf(stderr, "Removing duplicates ...done. Removed %" PRId64 " duplicates out of %" PRId64 " total particles\n",
            num_removed, halocat->totnpart);

    free(haloids);
    free(fofids);
}
