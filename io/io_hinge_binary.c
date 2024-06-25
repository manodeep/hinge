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

void loadgroups_hinge_binary(const struct params_data *params, const int snapnum, struct group_data *group)
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

    int interrupted = 0;
    fprintf(stderr, "Assigning group-level properties ...\n");
    init_my_progressbar(nhalos, &interrupted);
    int64_t fof_hostnum = 0;
    for (int64_t ihalo = 0; ihalo < nhalos; ihalo++)
    {
        my_progressbar(ihalo, &interrupted);
        group[ihalo].N = halocat->halos[ihalo].npart;
        group[ihalo].haloID = halocat->halos[ihalo].halo_id;
        group[ihalo].fofID = halocat->halos[ihalo].fof_id;
        group[ihalo].Nsub = halocat->halos[ihalo].nsub;
        group[ihalo].Mtot = halocat->halos[ihalo].Mvir;
        group[ihalo].xcen = halocat->halos[ihalo].Xc;
        group[ihalo].ycen = halocat->halos[ihalo].Yc;
        group[ihalo].zcen = halocat->halos[ihalo].Zc;
        group[ihalo].vxcen = halocat->halos[ihalo].VXc;
        group[ihalo].vycen = halocat->halos[ihalo].VYc;
        group[ihalo].vzcen = halocat->halos[ihalo].VZc;

        group[ihalo].groupnum = ihalo;
        group[ihalo].nodeloc = ihalo;
        group[ihalo].snapshot = snapnum;
        group[ihalo].redshift = REDSHIFT[snapnum];

        int64_t haloid = group[ihalo].haloID;
        int64_t hosthaloid = group[ihalo].fofID;
        if (haloid == hosthaloid)
        {
            fof_hostnum = ihalo;
            // fof_hostid = haloid;
        }
        group[ihalo].isFof = (haloid == hosthaloid) ? 1 : 0;
        XASSERT(fof_hostnum >= 0 && fof_hostnum < nhalos,
                "Invalid fofnum = %" PRId64 " for halo = %" PRId64 " hosthalo = %" STR_FMT " haloid = %" STR_FMT "\n",
                fof_hostnum, ihalo, hosthaloid, haloid);
        group[ihalo].FOFHalo = fof_hostnum;
        group[ihalo].ContainerIndex = fof_hostnum;
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

        group[ihalo].parentgroupforparticle = my_malloc(sizeof(*group[ihalo].parentgroupforparticle), group[ihalo].N);
        group[ihalo].parentsnapshotforparticle =
            my_malloc(sizeof(*group[ihalo].parentsnapshotforparticle), group[ihalo].N);
    }
    finish_myprogressbar(&interrupted);
    fprintf(stderr, "Assigning group-level properties ...done\n");

    fprintf(stderr, "Reading and assigning field: 'partid', 'xpos', 'ypos', 'zpos' ...\n");

#if 0
    /* read individual files for each column */
    /* Can't really automate the process - so need to read in individually and assign */
    // const char field_names[][MAXLEN] = {"partid", "xpos", "ypos", "zpos", "haloid", "fofid"};
#define CHECK_NPART_AND_READ_FIELD(field_name, totnpart, buf)                                                          \
    {                                                                                                                  \
        fprintf(stderr, " '%s', ", field_name);                                                                        \
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
        time_t _t0 = time(NULL);                                                                                       \
        fprintf(stderr, "Reading " field_name "...\n");                                                                \
        buf = my_malloc(sizeof_every_field, totnpart);                                                                 \
        my_fread(buf, sizeof_every_field, totnpart,                                                                    \
                 fp); /* all the fields are of type int64_t or double (i.e., 8 bytes) */                               \
        fclose(fp);                                                                                                    \
        time_t _t1 = time(NULL);                                                                                       \
        fprintf(stderr, "Reading " field_name "...done. Time taken = %0.3f s\n", difftime(_t1, _t0));                  \
    }

#define ASSIGN_FIELD_TO_GROUPS(field_name, field_type, nhalos, buf, dst_field)                                         \
    {                                                                                                                  \
        int64_t offset = 0;                                                                                            \
        for (int64_t i = 0; i < nhalos; i++)                                                                           \
        {                                                                                                              \
            const int64_t npart_field = group[i].N;                                                                    \
            group[i].dst_field = my_malloc(sizeof(field_type), npart_field);                                           \
            for (int64_t j = 0; j < npart_field; j++)                                                                  \
            {                                                                                                          \
                group[i].dst_field[j] = ((field_type *)buf)[offset + j];                                               \
            }                                                                                                          \
            offset += npart_field;                                                                                     \
        }                                                                                                              \
        free(buf);                                                                                                     \
    }

    void *buf;
    time_t t0 = time(NULL);
    fprintf(stderr, "Reading and assigning field: ");
    CHECK_NPART_AND_READ_FIELD("partid", totnpart, buf);
    ASSIGN_FIELD_TO_GROUPS("partid", int64_t, nhalos, buf, id);

    //    CHECK_NPART_AND_READ_FIELD(fp, "part_type", totnpart, buf);
    CHECK_NPART_AND_READ_FIELD("xpos", totnpart, buf);
    ASSIGN_FIELD_TO_GROUPS("xpos", double, nhalos, buf, x);

    CHECK_NPART_AND_READ_FIELD("ypos", totnpart, buf);
    ASSIGN_FIELD_TO_GROUPS("ypos", double, nhalos, buf, y);

    CHECK_NPART_AND_READ_FIELD("zpos", totnpart, buf);
    ASSIGN_FIELD_TO_GROUPS("zpos", double, nhalos, buf, z);
#else
#define OPEN_FILE_AND_CHECK_NPART(field_name, totnpart, fp_out)                                                        \
    {                                                                                                                  \
        int64_t npart_field;                                                                                           \
        char field_fname[MAXLEN];                                                                                      \
        my_snprintf(field_fname, MAXLEN, "%s/%s_particles_z%0.3f_%s.bin", params->GROUP_DIR, params->GROUP_BASE,       \
                    REDSHIFT[snapnum], field_name);                                                                    \
        fp_out = my_fopen(field_fname, "r");                                                                           \
        my_fread(&npart_field, sizeof(npart_field), 1, fp_out);                                                        \
        if (npart_field != totnpart)                                                                                   \
        {                                                                                                              \
            fprintf(stderr,                                                                                            \
                    "Number of particles in file %s = %" PRId64 " does not match total number of particles = %" PRId64 \
                    "\n",                                                                                              \
                    field_fname, npart_field, totnpart);                                                               \
            fclose(fp_out);                                                                                            \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    }
#define READ_INTO_BUF_OR_GROUP(thisgroup, dst_field, npart_field, fp_in, buf, buf_type)                                \
    {                                                                                                                  \
        const size_t sizeof_every_field = 8;                                                                           \
        if (sizeof(*thisgroup->dst_field) == sizeof_every_field)                                                       \
        {                                                                                                              \
            my_fread(thisgroup->dst_field, sizeof_every_field, npart_field, fp_in);                                    \
        }                                                                                                              \
        else                                                                                                           \
        {                                                                                                              \
            my_fread(buf, sizeof_every_field, npart_field, fp_in);                                                     \
            for (int64_t j = 0; j < npart_field; j++)                                                                  \
            {                                                                                                          \
                thisgroup->dst_field[j] = ((buf_type *)buf)[j];                                                        \
            }                                                                                                          \
        }                                                                                                              \
    }
    time_t t0 = time(NULL);
    void *buf = my_malloc(sizeof(double), totnpart);
    FILE *fp_xpos, *fp_ypos, *fp_zpos, *fp_partid;
    OPEN_FILE_AND_CHECK_NPART("xpos", totnpart, fp_xpos);
    OPEN_FILE_AND_CHECK_NPART("ypos", totnpart, fp_ypos);
    OPEN_FILE_AND_CHECK_NPART("zpos", totnpart, fp_zpos);
    OPEN_FILE_AND_CHECK_NPART("partid", totnpart, fp_partid);
    init_my_progressbar(nhalos, &interrupted);
    for (int64_t i = 0; i < nhalos; i++)
    {
        my_progressbar(i, &interrupted);
        const int64_t npart_field = group[i].N;
        struct group_data *thisgroup = &group[i];
        thisgroup->id = my_malloc(sizeof(*thisgroup->id), npart_field);
        thisgroup->x = my_malloc(sizeof(*thisgroup->x), npart_field);
        thisgroup->y = my_malloc(sizeof(*thisgroup->y), npart_field);
        thisgroup->z = my_malloc(sizeof(*thisgroup->z), npart_field);
        READ_INTO_BUF_OR_GROUP(thisgroup, x, npart_field, fp_xpos, buf, double);
        READ_INTO_BUF_OR_GROUP(thisgroup, y, npart_field, fp_ypos, buf, double);
        READ_INTO_BUF_OR_GROUP(thisgroup, z, npart_field, fp_zpos, buf, double);
        READ_INTO_BUF_OR_GROUP(thisgroup, id, npart_field, fp_partid, buf, int64_t);
    }
    finish_myprogressbar(&interrupted);
    free(buf);
    fclose(fp_xpos);
    fclose(fp_ypos);
    fclose(fp_zpos);
    fclose(fp_partid);
#endif

    time_t t1 = time(NULL);
    fprintf(stderr, "\rReading and assigning field: 'partid', 'xpos', 'ypos', 'zpos', 'haloid', 'fofid' ... done.\n");
    print_time(t0, t1, "Reading and assigning fields");

#ifdef CHECK_HALO_FOF_IDS
    int64_t *haloids, *fofids;
    CHECK_NPART_AND_READ_FIELD("haloid", totnpart, haloids);
    CHECK_NPART_AND_READ_FIELD("fofid", totnpart, fofids);

    int64_t offset = 0;
    fprintf(stderr, "Checking the consistency of the halo and fof ids ...\n");
    for (int64_t ihalo = 0; ihalo < nhalos; ihalo++)
    {
        const int64_t *fids = &fofids[offset];
        const int64_t *hids = &haloids[offset];
        const int64_t haloid = group[ihalo].haloID;
        const int64_t hosthaloid = group[ihalo].fofID;
        for (int64_t j = 0; j < group[ihalo].N; j++)
        {
            XASSERT(hids[j] == haloid, "Haloid mismatch: %" PRId64 " != %" PRId64, hids[j], haloid);
            XASSERT(fids[j] == hosthaloid, "Fofid mismatch: %" PRId64 " != %" PRId64, fids[j], hosthaloid);
        }

        offset += group[ihalo].N;
    }
    fprintf(stderr, "Checking the consistency of the halo and fof ids ...done\n");
    free(haloids);
    free(fofids);
#endif

    // fprintf(stderr, "Removing duplicates ...\n");
    t0 = time(NULL);
    const int64_t num_removed = remove_duplicates(group, nhalos);
    t1 = time(NULL);
    // fprintf(stderr, "Removing duplicates ...done. Removed %" PRId64 " duplicates out of %" PRId64 " total
    // particles\n",
    //         num_removed, halocat->totnpart);
    print_time(t0, t1, "Removing duplicates");

    free_hinge_halocat(halocat);
}
