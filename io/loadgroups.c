#include "loadgroups.h"
#include "utils.h"

// To do the re-ordering
#include "sglib.h"

// #ifdef BGC2
// #include "bgc2.h"
// #endif

// #include "read_param.h"

#include "io_hinge_ascii.h"
#include "io_hinge_binary.h"
#include "io_subfind_binary.h"

#define FOFMINLEN (32)
#define SUBHALOMINLEN (20)

#define SQR_PERIODIC(dx) (periodic(dx) * periodic(dx))

#if 0
void reorder_groups_on_array(const int64 Ngroups, struct group_data *group)
{
    // The last two fields, viz., parentgroupforparticle and
    // parentsnapshotforparticle,  are also particle level information but they do
    // not contain any info at this point. But I will swap them as well, since I
    // might (inadvertently) call this function from elsewhere.

    int64 Nallocated = Ngroups; // some "reasonable" estimate of the max. number
                                // of particles in a halo. Tune to suit your needs
    float *sqr_radius = NULL;

    for (int64 igroup = 0; igroup < Ngroups; igroup++)
    {

        if (group[igroup].N > Nallocated)
        {
            free(sqr_radius);
            sqr_radius = NULL;
            Nallocated = group[igroup].N;
        }

        if (sqr_radius == NULL)
        {
            sqr_radius = my_malloc(sizeof(*sqr_radius), Nallocated);
        }

        for (int64 j = 0; j < group[igroup].N; j++)
        {
            const float dx = periodic(group[igroup].xcen - group[igroup].x[j]);
            const float dy = periodic(group[igroup].ycen - group[igroup].y[j]);
            const float dz = periodic(group[igroup].zcen - group[igroup].z[j]);
            sqr_radius[j] = dx * dx + dy * dy + dz * dz;
        }
        struct group_data *thisgroup = &(group[igroup]);
#ifdef SUSSING_TREES

#define MULTIPLE_ARRAY_EXCHANGER(vartype, a, i, j)                                                                     \
    {                                                                                                                  \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, sqr_radius, i, j);                                                       \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(id64, thisgroup->id, i, j);                                                     \
        /* SGLIB_ARRAY_ELEMENTS_EXCHANGER(int, thisgroup->type,i,j);                                                   \
         */                                                                                                            \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->x, i, j);                                                     \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->y, i, j);                                                     \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->z, i, j);                                                     \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->ParticleEnergy, i, j);                                        \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->vx, i, j);                                                    \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->vy, i, j);                                                    \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->vz, i, j);                                                    \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(int64, thisgroup->parentgroupforparticle, i, j);                                \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(int, thisgroup->parentsnapshotforparticle, i, j);                               \
    }

#else
#define MULTIPLE_ARRAY_EXCHANGER(vartype, a, i, j)                                                                     \
    {                                                                                                                  \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, sqr_radius, i, j);                                                       \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(id64, thisgroup->id, i, j);                                                     \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->x, i, j);                                                     \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->y, i, j);                                                     \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->z, i, j);                                                     \
        /*   SGLIB_ARRAY_ELEMENTS_EXCHANGER(int, thisgroup->type,i,j);	*/                                              \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(int64, thisgroup->parentgroupforparticle, i, j);                                \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(int, thisgroup->parentsnapshotforparticle, i, j);                               \
    }
#endif

        /* SGLIB_ARRAY_QUICK_SORT(float, sqr_radius, thisgroup->N,
         * SGLIB_NUMERIC_COMPARATOR , MULTIPLE_ARRAY_EXCHANGER); */
        /*if the sort is taking too long, replace with HEAP_SORT -> the worst case
         * performance is better*/
        SGLIB_ARRAY_HEAP_SORT(float, sqr_radius, thisgroup->N, SGLIB_NUMERIC_COMPARATOR, MULTIPLE_ARRAY_EXCHANGER);
#undef MULTIPLE_ARRAY_EXCHANGER
    }

    if (sqr_radius != NULL)
    {
        free(sqr_radius);
    }
}

static inline void mark_particle_for_deletion(struct group_data *this_group, const int64 partloc)
{
    assert(this_group != NULL && "The group pointer can not be NULL while removing a particle");
    id64 *id = this_group->id;

    id[partloc] = -1;
}

static inline void swap_groups(struct group_data *groups, const int64 first, const int64 second)
{
    struct group_data tmp = groups[first];
    groups[first] = groups[second];
    groups[second] = tmp;
}

void remove_particles_from_all_groups(struct group_data *groups, const int Ngroups)
{
    assert(groups != NULL && "The group pointer can not be NULL");
    for (int64 igroup = 0; igroup < Ngroups; igroup++)
    {
        struct group_data *this_group = &(groups[igroup]);
        id64 *id = this_group->id;
        // int *type = this_group->type;
        float *x = this_group->x;
        float *y = this_group->y;
        float *z = this_group->z;

#ifdef SUSSING_TREES
        float *ParticleEnergy = this_group->ParticleEnergy;
        float *vx = this_group->vx;
        float *vy = this_group->vy;
        float *vz = this_group->vz;
#endif

        int64 *parentgroupforparticle = this_group->parentgroupforparticle;
        int *parentsnapshotforparticle = this_group->parentsnapshotforparticle;

        int64 num_moved = 0;
        const int64 orig_N = this_group->N;
        int64 j = 0;
        while (j < this_group->N)
        {
            if (id[j] != -1)
            {
                j++;
                continue;
            }
            const int64 source = this_group->N - 1;
            id[j] = id[source];
            // type[j] = type[source];
            x[j] = x[source];
            y[j] = y[source];
            z[j] = z[source];

#ifdef SUSSING_TREES
            ParticleEnergy[j] = ParticleEnergy[source];
            vx[j] = vx[source];
            vy[j] = vy[source];
            vz[j] = vz[source];
#endif

            parentgroupforparticle[j] = parentgroupforparticle[source];
            parentsnapshotforparticle[j] = parentsnapshotforparticle[source];

            num_moved++;
            this_group->N--;
        }
        if (num_moved > 0)
        {
            assert(this_group->N + num_moved == orig_N && "Sum of particles removed and particles left should equal "
                                                          "original number of particles");
            /* fprintf(stderr,"igroup = %"STR_FMT" number of duplicates removed =
             * %"STR_FMT" original Npart = %"STR_FMT" current N =
             * %"STR_FMT"\n",igroup,num_moved,orig_N,this_group->N); */
        }
    }

    // Now check that the removal worked
    for (int64 igroup = 0; igroup < Ngroups; igroup++)
    {
        struct group_data *this_group = &(groups[igroup]);
        id64 *id = this_group->id;
        for (int64 j = 0; j < this_group->N; j++)
        {
            if (!(*id != -1))
            {
                fprintf(stderr, "ERROR: about to crash\n");
                fprintf(stderr, "igroup = %" STR_FMT " j = %" STR_FMT " this_group->N = %" STR_FMT "\n", igroup, j,
                        this_group->N);
            }
            assert(*id != -1 && "id's are correct");
            id++;
        }
    }
}

/* This function should be called after the halos have been
         sorted such that fofs and subs are contiguous.
 */
void remove_duplicate_particles(const int64 Ngroups, struct group_data *groups)
{
    int64 FofIndex = 0; // actual index and not a fof haloid since that can be different
    const int mem_increase_fac = 1.1;
    int64 max_npart = 10000000;
    long *partids = my_malloc(sizeof(*partids), max_npart);
    short *halolevel = my_malloc(sizeof(*halolevel), max_npart);
    int64 *haloindex = my_malloc(sizeof(*haloindex), max_npart);
    int64 *partloc_in_halo = my_malloc(sizeof(*partloc_in_halo), max_npart);
    while (FofIndex < Ngroups)
    {
        int64 location = 0;
        /* First take all of the particle ids, the groups they are located in and
         * the total number of particles in the groups */
        for (int64 igroup = FofIndex; igroup < (FofIndex + groups[FofIndex].Nsub); igroup++)
        {
            if (!(groups[igroup].FOFHalo == FofIndex))
            {
                fprintf(stderr,
                        "ERROR: ABout to crash. igroup = %" STR_FMT " FofIndex = %" STR_FMT
                        " group->FOFHalo = %" STR_FMT "\n",
                        igroup, FofIndex, groups[igroup].FOFHalo);
            }
            assert(groups[igroup].FOFHalo == FofIndex && "The FOFHalo indices should have been setup correctly");

            for (int64 j = 0; j < groups[igroup].N; j++)
            {
                const long this_id = groups[igroup].id[j];
                if (location == max_npart)
                {
                    // reallocate memory
                    max_npart *= mem_increase_fac;
                    while (max_npart <= location)
                        max_npart += 1000;

                    partids = my_realloc(partids, sizeof(*partids), max_npart, "partids");
                    halolevel = my_realloc(halolevel, sizeof(*halolevel), max_npart, "halolevel");
                    haloindex = my_realloc(haloindex, sizeof(*haloindex), max_npart, "haloindex");
                    partloc_in_halo =
                        my_realloc(partloc_in_halo, sizeof(*partloc_in_halo), max_npart, "partloc_in_halo");
                }

                /* if(this_id == 2047766) { */
                /* 	fprintf(stderr,"this_id = %ld igroup = %"STR_FMT" j = %"STR_FMT"
                 * halolevel = %d group->N = %"STR_FMT" containerindex = %"STR_FMT"
                 * groupnum = %"STR_FMT"\n", */
                /* 					this_id,igroup,j,groups[igroup].ParentLevel,groups[igroup].N,groups[igroup].ContainerIndex,groups[igroup].groupnum);
                 */
                /* } */

                partids[location] = this_id;
                // Make sure that Parentlevel is set up before
                assert(groups[igroup].ParentLevel >= 1 && "ParentLevels need to have been initialized already");
                halolevel[location] = groups[igroup].ParentLevel;
                haloindex[location] = igroup;
                partloc_in_halo[location] = j;
                location++;
            }
        }
        /* fprintf(stderr,"haloindex[%"STR_FMT"-1] =
         * %"STR_FMT"\n",location,haloindex[location-1]); */
        const int64 Npart = location;

#define MULTIPLE_ARRAY_EXCHANGER(vartype, a, i, j)                                                                     \
    {                                                                                                                  \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(short, halolevel, i, j);                                                        \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(long, partids, i, j);                                                           \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(int64, haloindex, i, j);                                                        \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(int64, partloc_in_halo, i, j);                                                  \
    }

        SGLIB_ARRAY_QUICK_SORT(id64, partids, Npart, SGLIB_NUMERIC_COMPARATOR, MULTIPLE_ARRAY_EXCHANGER);
#undef MULTIPLE_ARRAY_EXCHANGER

        int64 ipart = 0;
        while (ipart < Npart)
        {
            const long this_partid = partids[ipart];
            int64 nmatches = 0;
            int64 partindex = ipart + 1;
            short maxparentlevel = halolevel[ipart];
            /* int64 haloindex_maxparentlevel = haloindex[ipart]; */
            while (partindex < Npart && partids[partindex] == this_partid)
            {
                if (halolevel[partindex] > maxparentlevel)
                {
                    maxparentlevel = halolevel[partindex];
                    /* haloindex_maxparentlevel = haloindex[partindex]; */
                }
                nmatches++;
                partindex++;
            }

            // There are no duplicates -> update ipart and continue with the loop
            if (nmatches == 0)
            {
                ipart++;
                continue;
            }

            /* So we have duplicate particles -> Need to remove this particles from at
             * least one halo */
            /* fprintf(stderr,"ipart = %"STR_FMT" nmatches = %"STR_FMT" this_id =
             * %"STR_ID_FMT"\n",ipart,nmatches,this_partid); */
            for (int64 i = ipart; i <= (ipart + nmatches); i++)
            {
                /* fprintf(stderr,"partids[%"STR_FMT"] = %"STR_ID_FMT" this_partid =
                 * %"STR_ID_FMT"\n",i,partids[i],this_partid); */
                assert(partids[i] == this_partid && "Particle ids must be identical");
                // Check that the same particle does not belong to two different
                // subhalos at the same parentlevel MS 3rd Mar, 2015 - this condition
                // does not hold. Two different subhalos may contain the same particle.
                /* if(halolevel[i] == maxparentlevel) { */
                /* 	if( !(haloindex_maxparentlevel == haloindex[i])) { */
                /* 		fprintf(stderr,"ERROR: About to crash: partid = %"STR_ID_FMT"
                 * haloindex[%"STR_FMT"] = %"STR_FMT"  haloindex_maxparentlevel =
                 * %"STR_FMT"  halolevel[%"STR_FMT"] = %d maxparentlevel = %d\n", */
                /* 						partids[i],
                 * i,haloindex[i],haloindex_maxparentlevel,i,halolevel[i],maxparentlevel);
                 */
                /* 	} */
                /* 	assert(haloindex_maxparentlevel == haloindex[i] && "The same
                 * particle can not be in two different halos at the same parent
                 * level"); */
                /* } */

                /* Do we need to remove this particle from its group */
                if (halolevel[i] < maxparentlevel)
                {
                    // Remove this particle from its group
                    const int64 this_haloindex = haloindex[i];
                    assert(this_haloindex >= 0 && this_haloindex < (FofIndex + groups[FofIndex].Nsub));
                    struct group_data *this_group = &(groups[this_haloindex]);
                    const int64 partloc = partloc_in_halo[i];
                    assert(partloc >= 0 && partloc < this_group->N);
                    /* fprintf(stderr,"i = %"STR_FMT" before removing haloindex =
                     * %"STR_FMT" group->N =
                     * %"STR_FMT"\n",i,this_haloindex,this_group->N); */
                    mark_particle_for_deletion(this_group, partloc);
                    /* fprintf(stderr,"i = %"STR_FMT " removed particle %"STR_FMT" from
                     * haloindex = %"STR_FMT" group->N =
                     * %"STR_FMT"\n",i,partloc,this_haloindex,this_group->N); */
                }
            }
            ipart = ipart + nmatches + 1;
        }

        assert(groups[FofIndex].Nsub >= 1 && "Number of subhalos must be at least 1");
        FofIndex += groups[FofIndex].Nsub;
    }
    assert(FofIndex == Ngroups && "The number of halos processed should be *exactly* equal to the "
                                  "number of halos present");

    free(partids);
    free(halolevel);
    free(haloindex);
    free(partloc_in_halo);

    // The duplciate particles have already been marked -- now remove them
    remove_particles_from_all_groups(groups, Ngroups);
}


// #ifndef SUSSING_TREES
// #ifndef ASCII_DATA
// #ifndef BGC2
// int64 returnNhalo(const char *buf)
// {
//     FILE *fsub = NULL;
//     int64 Nsub;
//     fsub = my_fopen(buf, "r");
//     my_fread(&Nsub, sizeof(int64), 1, fsub);
//     fclose(fsub);
//     return Nsub;
// }
// #endif
// #endif
// #endif

#endif // if 0

int64 returnNhalo(const struct params_data *params, const int snapnum, const int fof_only)
{
    if (params == NULL)
    {
        fprintf(stderr, "%s>: params is NULL\n", __FUNCTION__);
        return -1;
    }
    if (snapnum < 0 || snapnum > params->MAX_SNAPSHOT_NUM)
    {
        fprintf(stderr, "%s>: snapnum = %d is out of range [0, %d]\n", __FUNCTION__, snapnum, params->MAX_SNAPSHOT_NUM);
        return -1;
    }

    switch (params->GROUP_FORMAT)
    {
    case subfind_binary:
        return returnNhalo_subfind_binary(params, snapnum, fof_only);
    case hinge_ascii:
        return returnNhalo_hinge_ascii(params, snapnum, fof_only);
    case hinge_binary:
        return returnNhalo_hinge_binary(params, snapnum, fof_only);

    default:
        fprintf(stderr, "%s>: Unknown group format = %d\n", __FUNCTION__, (int)params->GROUP_FORMAT);
        return -1;
    }

    return -1;
}

void loadgroups(const struct params_data *params, const int snapnum, struct group_data *group)
{
    if (params == NULL || group == NULL)
    {
        fprintf(stderr, "%s>: params or group is NULL\n", __FUNCTION__);
        exit(EXIT_FAILURE);
    }
    if(params->flag_load_only_partids == 1 && params->GROUP_FORMAT != hinge_binary)
    {
        fprintf(stderr,"Error: flag_load_only_partids is set to 1 but the group format is not hinge_binary\n");
        exit(EXIT_FAILURE);
    }

    switch (params->GROUP_FORMAT)
    {
    case subfind_binary:
        loadgroups_subfind_binary(snapnum, params, group);
        break;
    case hinge_ascii:
        loadgroups_hinge_ascii(snapnum, params, group);
        break;

    case hinge_binary:
        loadgroups_hinge_binary(snapnum, params, group);
        break;

    default:
        fprintf(stderr, "%s>: Unknown group format = %d\n", __FUNCTION__, (int)params->GROUP_FORMAT);
        exit(EXIT_FAILURE);
    }

    return;
}

#undef FOFMINLEN
#undef SUBHALOMINLEN
