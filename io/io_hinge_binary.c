#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#define _FILE_OFFSET_BITS 64
#include <sys/types.h>
#include <unistd.h>

#include <fcntl.h> //for open and close

#if __APPLE__
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/uio.h>
#else
#include <sys/sendfile.h>
#endif

#include "io_hinge_ascii.h"
#include "io_hinge_binary.h"
#include "io_hinge_utils.h"
#include "macros.h"
#include "progressbar.h"
#include "utils.h"

void save_unique_particles(const struct params_data *params, const int snapnum, struct group_data *group,
                           const int64 nhalos);
void load_unique_particles(struct params_data *params, const int snapnum, struct group_data *group);

int64 returnNhalo_hinge_binary(const struct params_data *params, const int snapnum, const int fof_only)
{
    return returnNhalo_hinge_ascii(params, snapnum, fof_only);
}

void loadgroups_hinge_binary(struct params_data *params, const int snapnum, struct group_data *group)
{
    XASSERT(group != NULL, "group is NULL\n");

    if (params->LOAD_UNIQUE_PARTICLES)
    {
        fprintf(stderr, "Calling load_unique_particles\n");
        return load_unique_particles(params, snapnum, group);
    }
    else
    {
        fprintf(stderr, "Not calling load_unique_particles. params->load_unique_particles = %d\n",
                params->LOAD_UNIQUE_PARTICLES);
    }

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
            buf_type *_tmp_buf = (buf_type *)buf;                                                                      \
            for (int64_t j = 0; j < npart_field; j++)                                                                  \
            {                                                                                                          \
                thisgroup->dst_field[j] = *_tmp_buf++;                                                                 \
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
    init_my_progressbar(totnpart, &interrupted);
    int64_t npart_read = 0;
    for (int64_t i = 0; i < nhalos; i++)
    {
        my_progressbar(npart_read, &interrupted);
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
        npart_read += npart_field;
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
    // const int64_t num_removed = remove_duplicates(group, nhalos);
    remove_duplicates(group, nhalos);
    t1 = time(NULL);
    // fprintf(stderr, "Removing duplicates ...done. Removed %" PRId64 " duplicates out of %" PRId64 " total
    // particles\n",
    //         num_removed, halocat->totnpart);
    print_time(t0, t1, "Removing duplicates");

    free_hinge_halocat(halocat);

    if (params->SAVE_UNIQUE_PARTICLES)
    {
        fprintf(stderr, "Saving unique particles ...\n");
        t0 = time(NULL);
        save_unique_particles(params, snapnum, group, nhalos);
        t1 = time(NULL);
        fprintf(stderr, "Saving unique particles ...done. Saved %" STR_FMT " halos\n", nhalos);
        print_time(t0, t1, "Saving unique particles");
    }
}

void save_unique_particles(const struct params_data *params, const int snapnum, struct group_data *group,
                           const int64 nhalos)
{
    char unique_fname[MAXLEN];
    fprintf(stderr, "Saving unique particles ... nhalos = %" PRId64 "\n", nhalos);
    my_snprintf(unique_fname, MAXLEN, "%s/%s_unique_particles_partids_z%0.3f.bin", params->OUTPUT_DIR,
                params->GROUP_BASE, REDSHIFT[snapnum]);
    FILE *fp_ids = my_fopen(unique_fname, "w+");
    fwrite(&nhalos, sizeof(nhalos), 1, fp_ids);

    my_snprintf(unique_fname, MAXLEN, "%s/%s_unique_particles_xpos_z%0.3f.bin", params->OUTPUT_DIR, params->GROUP_BASE,
                REDSHIFT[snapnum]);
    FILE *fp_xpos = my_fopen(unique_fname, "w+");
    fwrite(&nhalos, sizeof(nhalos), 1, fp_xpos);

    my_snprintf(unique_fname, MAXLEN, "%s/%s_unique_particles_ypos_z%0.3f.bin", params->OUTPUT_DIR, params->GROUP_BASE,
                REDSHIFT[snapnum]);
    FILE *fp_ypos = my_fopen(unique_fname, "w+");
    fwrite(&nhalos, sizeof(nhalos), 1, fp_ypos);

    my_snprintf(unique_fname, MAXLEN, "%s/%s_unique_particles_zpos_z%0.3f.bin", params->OUTPUT_DIR, params->GROUP_BASE,
                REDSHIFT[snapnum]);
    FILE *fp_zpos = my_fopen(unique_fname, "w+");
    fwrite(&nhalos, sizeof(nhalos), 1, fp_zpos);

    my_snprintf(unique_fname, MAXLEN, "%s/%s_unique_particles_halocat_z%0.3f.bin", params->OUTPUT_DIR,
                params->GROUP_BASE, REDSHIFT[snapnum]);
    const size_t sizeof_group_data = sizeof(struct group_data);
    FILE *fp_cat = my_fopen(unique_fname, "w+");
    fwrite(&nhalos, sizeof(nhalos), 1, fp_cat);
    fwrite(&sizeof_group_data, sizeof(sizeof_group_data), 1, fp_cat);
    int64 totnpart = 0, totnpart_all = 0;
    for (int64 i = 0; i < nhalos; i++)
    {
        //Creating a copy to save because otherwise the current version in memory will get impacted
        //by the change in group.N
        struct group_data thisgroup;
        memcpy(&thisgroup, &group[i], sizeof(struct group_data));
        if (thisgroup.haloID < 0)
        {
            /* While it may be tempting to do 'nhalos--' here, doing so
            would mean that we would to refix *all* the fofnum/hostnum etc, i.e., anything
            that is an index to the group catalog will need to be fixed.
            Certainly doable but not worth the effort for now. MS: 26th June, 2024
            */
            thisgroup.N = 0;
        }
        totnpart_all += thisgroup.N;
        int64 num_dups = 0;
        for (int64 j = 0; j < thisgroup.N; j++)
        {
            num_dups += thisgroup.id[j] < 0 ? 1 : 0;
        }

#define _WRITE_ARRAY_ELEMENTS(start, num_left)                       \
{                                                                           \
    fwrite(&thisgroup.id[start], sizeof(*thisgroup.id), num_left, fp_ids);  \
    fwrite(&thisgroup.x[start], sizeof(*thisgroup.x), num_left, fp_xpos);   \
    fwrite(&thisgroup.y[start], sizeof(*thisgroup.y), num_left, fp_ypos);   \
    fwrite(&thisgroup.z[start], sizeof(*thisgroup.z), num_left, fp_zpos);   \
}
#define WRITE_REMAINING_WHEN_NUM_DUPS_ZERO(num_dups, thisgroup, start, totN, num_left)                                     \
    {                                                                                                            \
        XASSERT(num_dups == 0, "Error: num_dups = %" PRId64 "\n", num_dups);                                     \
        XASSERT(start + num_left == totN, "Error: start = %" PRId64 " num_left = %" PRId64 " totN = %" PRId64 "\n", \
                start, num_left, totN);                                                                          \
        XASSERT(start >=0 && start < totN, "Error: start = %" PRId64 " totN = %" PRId64 "\n", start, totN);      \
        XASSERT(num_left > 0 && num_left <= totN, "Error: num_left = %" PRId64 "\n", num_left);                  \
        _WRITE_ARRAY_ELEMENTS(start, num_left);                                                                  \
    }


        if (num_dups == 0)
        {
            WRITE_REMAINING_WHEN_NUM_DUPS_ZERO(num_dups, thisgroup, 0, thisgroup.N, thisgroup.N);
        }
        else
        {
            int64 num_written = 0, num_dups_remaining = num_dups;
            for (int64 j = 0; j < thisgroup.N; j++)
            {
                if (thisgroup.id[j] < 0)
                {
                    num_dups_remaining--;
                    if (num_dups_remaining == 0 && j < (thisgroup.N - 1))
                    {
                        const int64 num_left = thisgroup.N - 1 - j;
                        if (num_left > 0)
                        {
                            WRITE_REMAINING_WHEN_NUM_DUPS_ZERO(num_dups_remaining, thisgroup, j + 1, thisgroup.N, num_left);
                            num_written += num_left;
                        }
                        break;
                    }
                    continue;
                }
                _WRITE_ARRAY_ELEMENTS(j, 1);
                // fwrite(&thisgroup.id[j], sizeof(thisgroup.id[j]), 1, fp_ids);
                // fwrite(&thisgroup.x[j], sizeof(thisgroup.x[j]), 1, fp_xpos);
                // fwrite(&thisgroup.y[j], sizeof(thisgroup.y[j]), 1, fp_ypos);
                // fwrite(&thisgroup.z[j], sizeof(thisgroup.z[j]), 1, fp_zpos);
                num_written++;
            }
            thisgroup.N -= num_dups;
            XASSERT(num_written == thisgroup.N,
                    "Error: For halo number %" PRId64 " number of particles written = %" PRId64 " != %" PRId64
                    ". Number of duplicates = %" PRId64 "\n",
                    i, num_written, thisgroup.N, num_dups);
            XASSERT(num_dups_remaining == 0, "Number of duplicates left = %" PRId64 "\n", num_dups_remaining);
        }
        fwrite(&thisgroup, sizeof_group_data, 1, fp_cat);
        totnpart += thisgroup.N;
    }
    fprintf(stderr, "In %s> totnpart = %" PRId64 " totnpart_all = %" PRId64 "\n", __FUNCTION__, totnpart, totnpart_all);
    fclose(fp_cat);

    fflush(fp_ids);
    fflush(fp_xpos);
    fflush(fp_ypos);
    fflush(fp_zpos);

    /* concat all the files together */
    my_snprintf(unique_fname, MAXLEN, "%s/%s_unique_particles_allprops_z%0.3f.bin", params->OUTPUT_DIR,
                params->GROUP_BASE, REDSHIFT[snapnum]);
    int fout = open(unique_fname, O_CREAT | O_WRONLY, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if (fout < 0)
    {
        fprintf(stderr, "Could not open file %s for writing\n", unique_fname);
        exit(EXIT_FAILURE);
    }
    write(fout, &nhalos, sizeof(nhalos));
    write(fout, &totnpart, sizeof(totnpart));

#define USE_SENDFILE_TO_WRITE_PROPS(fd_out, fd_in, ptr_offset, len)                                                    \
    {                                                                                                                  \
        ssize_t tot_nbytes_written = 0;                                                                                \
        while (len > 0)                                                                                                \
        {                                                                                                              \
            ssize_t nbytes_written = sendfile(fd_out, fd_in, ptr_offset, len);                                         \
            XASSERT(nbytes_written >= 0, "Error writing to file %s. nbytes_written = %zd\n", unique_fname,             \
                    nbytes_written);                                                                                   \
            len -= nbytes_written;                                                                                     \
            tot_nbytes_written += nbytes_written;                                                                      \
        }                                                                                                              \
        XASSERT(len == 0, "Error: len = %zu\n", len);                                                                  \
    }

    off_t start_offset = sizeof(int64); // to skip over numpart (of type int64) at the start of each file
    size_t len = totnpart * sizeof(group->id[0]);
    USE_SENDFILE_TO_WRITE_PROPS(fout, fileno(fp_ids), &start_offset, len);

    len = totnpart * sizeof(group->x[0]);
    start_offset = sizeof(int64); // since USE_SENDFILE can modify start_offset, start_offset needs to be set everytime
                                  // before calling USE_SENDFILE_TO_WRITE_PROPS
    USE_SENDFILE_TO_WRITE_PROPS(fout, fileno(fp_xpos), &start_offset, len);

    len = totnpart * sizeof(group->y[0]);
    start_offset = sizeof(int64); // since USE_SENDFILE can modify start_offset, start_offset needs to be set everytime
                                  // before calling USE_SENDFILE_TO_WRITE_PROPS
    USE_SENDFILE_TO_WRITE_PROPS(fout, fileno(fp_ypos), &start_offset, len);

    start_offset = sizeof(int64); // since USE_SENDFILE can modify start_offset, start_offset needs to be set everytime
                                  // before calling USE_SENDFILE_TO_WRITE_PROPS
    len = totnpart * sizeof(group->z[0]);
    USE_SENDFILE_TO_WRITE_PROPS(fout, fileno(fp_zpos), &start_offset, len);

    close(fout);
    fclose(fp_ids);
    fclose(fp_xpos);
    fclose(fp_ypos);
    fclose(fp_zpos);

    /* Now delete the individual files */
    my_snprintf(unique_fname, MAXLEN, "%s/%s_unique_particles_partids_z%0.3f.bin", params->OUTPUT_DIR,
                params->GROUP_BASE, REDSHIFT[snapnum]);
    unlink(unique_fname);

    my_snprintf(unique_fname, MAXLEN, "%s/%s_unique_particles_xpos_z%0.3f.bin", params->OUTPUT_DIR, params->GROUP_BASE,
                REDSHIFT[snapnum]);
    unlink(unique_fname);

    my_snprintf(unique_fname, MAXLEN, "%s/%s_unique_particles_ypos_z%0.3f.bin", params->OUTPUT_DIR, params->GROUP_BASE,
                REDSHIFT[snapnum]);
    unlink(unique_fname);

    my_snprintf(unique_fname, MAXLEN, "%s/%s_unique_particles_zpos_z%0.3f.bin", params->OUTPUT_DIR, params->GROUP_BASE,
                REDSHIFT[snapnum]);
    unlink(unique_fname);

    return;
}

void load_unique_particles(struct params_data *params, const int snapnum, struct group_data *group)
{
    char catalog_fname[MAXLEN];
    char unique_fname[MAXLEN];
    my_snprintf(catalog_fname, MAXLEN, "%s/%s_unique_particles_halocat_z%0.3f.bin", params->OUTPUT_DIR,
                params->GROUP_BASE, REDSHIFT[snapnum]);
    my_snprintf(unique_fname, MAXLEN, "%s/%s_unique_particles_allprops_z%0.3f.bin", params->OUTPUT_DIR,
                params->GROUP_BASE, REDSHIFT[snapnum]);
    fprintf(stderr, "Loading unique particles (catalog = '%s', particles = '%s')...\n", catalog_fname, unique_fname);

    time_t t0 = time(NULL);
    FILE *fp = NULL;
    FILE *fp_cat = fopen(catalog_fname, "rb");
    if (fp_cat == NULL)
        goto error;

    fp = fopen(unique_fname, "rb");
    if (fp == NULL)
        goto error;

    int64 nhalos;
    int64 status = fread(&nhalos, sizeof(nhalos), 1, fp_cat);
    if (status != 1)
        goto error;

    size_t sizeof_group_data;
    status = fread(&sizeof_group_data, sizeof(sizeof_group_data), 1, fp_cat);
    if (status != 1)
        goto error;

    if (sizeof_group_data != sizeof(struct group_data))
    {
        fprintf(stderr, "Warning: sizeof_group_data = %zu != %zu\n", sizeof_group_data, sizeof(struct group_data));
        fprintf(stderr, "Calling loadgroups again with `save_unique_particles`\n");
        params->SAVE_UNIQUE_PARTICLES = 1;
        goto error;
    }

    // read the group catalog
    status = fread(group, sizeof_group_data, nhalos, fp_cat);
    if (status != nhalos)
        goto error;

    // check nhalos from the unique particles properties file
    int64 nhalos_check;
    status = fread(&nhalos_check, sizeof(nhalos_check), 1, fp);
    if (status != 1)
        goto error;
    XASSERT(nhalos == nhalos_check,
            "nhalos (from file '%s')= %" PRId64 " should be equal to nhalos %" PRId64 " (from file '%s')\n",
            catalog_fname, nhalos, nhalos_check, unique_fname);

    int64 totnpart;
    status = fread(&totnpart, sizeof(totnpart), 1, fp);
    if (status != 1)
        goto error;

    int interrupted = 0;
    fprintf(stderr, "Reading and assigning field (from unique particles file): 'partid', 'xpos', 'ypos', 'zpos' ...\n");
    init_my_progressbar(totnpart, &interrupted);
    int64 numpart_read = 0;
    for (int64 i = 0; i < nhalos; i++)
    {
        my_progressbar(numpart_read, &interrupted);
        struct group_data *thisgroup = &group[i];
        thisgroup->id = my_malloc(sizeof(*thisgroup->id), thisgroup->N);
        thisgroup->x = my_malloc(sizeof(*thisgroup->x), thisgroup->N);
        thisgroup->y = my_malloc(sizeof(*thisgroup->y), thisgroup->N);
        thisgroup->z = my_malloc(sizeof(*thisgroup->z), thisgroup->N);
        status = fread(thisgroup->id, sizeof(*thisgroup->id), thisgroup->N, fp);
        if (status != thisgroup->N)
            goto error;

        status = fread(thisgroup->x, sizeof(*thisgroup->x), thisgroup->N, fp);
        if (status != thisgroup->N)
            goto error;

        status = fread(thisgroup->y, sizeof(*thisgroup->y), thisgroup->N, fp);
        if (status != thisgroup->N)
            goto error;

        status = fread(thisgroup->z, sizeof(*thisgroup->z), thisgroup->N, fp);
        if (status != thisgroup->N)
            goto error;

        thisgroup->parentgroupforparticle = my_malloc(sizeof(*thisgroup->parentgroupforparticle), thisgroup->N);
        thisgroup->parentsnapshotforparticle = my_malloc(sizeof(*thisgroup->parentsnapshotforparticle), thisgroup->N);

        numpart_read += thisgroup->N;
    }
    finish_myprogressbar(&interrupted);
    fprintf(stderr,
            "Reading and assigning field (from unique particles file): 'partid', 'xpos', 'ypos', 'zpos' ...done\n");
    time_t t1 = time(NULL);
    print_time(t0, t1, "Reading and assigning fields (from unique particle files)");
    fclose(fp_cat);
    fclose(fp);
    fprintf(stderr, "Loading unique particles (catalog = '%s', particles = '%s')...done\n", catalog_fname,
            unique_fname);

    return;

error:
    fprintf(stderr, "In %s> Encountered some error while reading `%s' or `%s'. Calling `loadgroups' again\n",
            __FUNCTION__, catalog_fname, unique_fname);
    if (fp_cat != NULL)
        fclose(fp_cat);
    if (fp != NULL)
        fclose(fp);
    params->LOAD_UNIQUE_PARTICLES =
        0; /* This is critical to "unset". Otherwise, infinite loop will occur loadgroups<->load_unique */
    return loadgroups_hinge_binary(params, snapnum, group);
}
