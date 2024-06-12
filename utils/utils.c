#include "utils.h"
#include "defs.h"
#include "macros.h"
#include "progressbar.h"
#include "read_param.h"
#include "sglib.h"

static void remove_particle_from_group(const int64 group1, const int64 group2, const int64 part1, const int64 part2,
                                       struct group_data *g, int64 *group_to_remove, int64 *part_to_remove);

// A real wrapper to snprintf that will exit() if the allocated buffer length
// was not sufficient. Usage is the same as snprintf
int my_snprintf(char *buffer, int len, const char *format, ...)
{
    va_list args;
    int nwritten = 0;

    va_start(args, format);
    nwritten = vsnprintf(buffer, len, format, args);
    va_end(args);
    if (nwritten > len || nwritten < 0)
    {
        fprintf(stderr,
                "ERROR: printing to string failed (wrote %d characters while only "
                "%d characters were allocated)\n",
                nwritten, len);
        fprintf(stderr, "Increase maxlen in `defs.h' ..exiting\n");
        exit(EXIT_FAILURE);
    }
    return nwritten;
}

int64 getnumlines(const char *fname, const char comment)
{
    FILE *fp = NULL;
    int64 nlines = 0;
    char str_line[MAXLINESIZE];

    fp = my_fopen(fname, "rt");

    while (1)
    {
        if (fgets(str_line, MAXLINESIZE, fp) != NULL)
        {
            // WARNING: this does not remove white-space. You might
            // want to implement that (was never an issue for me)
            if (str_line[0] != comment)
                nlines++;
        }
        else
            break;
    }
    fclose(fp);
    return nlines;
}

FILE *my_fopen(const char *fname, const char *mode)
{
    FILE *fp = NULL;
    fp = fopen(fname, mode);
    if (fp == NULL)
    {
        fprintf(stderr, "Could not open file %s\n", fname);
        exit(EXIT_FAILURE);
    }

    return fp;
}

FILE *my_fopen_carefully(const char *fname, void (*header)(FILE *))
{
    FILE *fp = NULL;
    fp = fopen(fname, "r");
    if (fp == NULL)
    {
        /*file does not exist -> can open with "w" */
        fp = my_fopen(fname, "w");
        (*header)(fp); /* print the header */
    }
    else
    {
        fclose(fp);
        fp = my_fopen(fname, "a+");
    }

    return fp;
}

void move_existing_file(const char *fname)
{
    FILE *fp = NULL;
    char commandstring[MAXLEN];
    fp = fopen(fname, "r");
    if (fp != NULL)
    {
        fclose(fp);
        fprintf(stderr, "Warning: file `%s' already exists\n", fname);
        my_snprintf(commandstring, MAXLEN, "mv %s  %s.bak", fname, fname);
        fprintf(stderr, "executing command `%s'\n", commandstring);
        system(commandstring); // probably should check for exit status..
    }
}

void my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
    size_t nread;
    nread = fread(ptr, size, nmemb, stream);
    if (nread != nmemb)
    {
        fprintf(stderr, "I/O error (fread) has occured.\n");
        fprintf(stderr, "Instead of reading nmemb=%zu, I got nread = %zu ..exiting\n", nmemb, nread);
        exit(EXIT_FAILURE);
    }
}

void print_time(const time_t t0, const time_t t1, const char *s)
{
    double timediff = difftime(t1, t0);
    double ratios[] = {24 * 3600.0, 3600.0, 60.0, 1};
    char units[4][10] = {"days", "hrs", "mins", "secs"};
    int which = 0;

    double timeleft = timediff;
    double time_to_print;
    fprintf(stderr, "\n Time taken to execute '%s'  = ", s);

    if (timediff < ratios[2])
        fprintf(stderr, "%5d secs", (int)timediff);
    else
        while (which < 4)
        {
            time_to_print = floor(timeleft / ratios[which]);
            if (time_to_print > 1)
            {
                timeleft -= (time_to_print * ratios[which]);
                fprintf(stderr, "%5d %s", (int)time_to_print, units[which]);
            }
            which++;
        }

    fprintf(stderr, "\n");
}

void *my_malloc(size_t size, int64 N)
{
    void *x = NULL;
    x = malloc(N * size);
    /* unsigned long megabytes = (unsigned long)N * size / (1024.0 * 1024.0); */
    if (x == NULL)
    {
        fprintf(stderr, "malloc for %" STR_FMT "  elements with %zu size failed..aborting\n", N, size);
        exit(EXIT_FAILURE);
    }
    /* else if (megabytes > 100) */
    /*     fprintf(stderr, "\n Successfully allocated  %" STR_FMT " elements with total size %lu (MB) \n", N,
     * megabytes); */

    return x;
}

void *my_calloc(size_t size, int64 N)
{
    size_t n = (size_t)N;
    void *x = NULL;
    x = calloc(n, size);
    if (x == NULL)
    {
        fprintf(stderr, "malloc for %" STR_FMT " elements with %zu size failed..aborting\n", N, size);
        exit(EXIT_FAILURE);
    }

    return x;
}

void *my_realloc(void *x, size_t size, int64_t N, const char *varname)
{
    void *tmp = realloc(x, N * size);
    /* size_t gigabytes = N * size / (1024.0 * 1024.0 * 1024.0); */

    if (tmp == NULL)
    {
        fprintf(stderr,
                "ERROR: Could not reallocate for %" PRId64 " elements with %zu size for variable `%s' ..aborting\n", N,
                size, varname);
        my_free((void **)&x);
        exit(EXIT_FAILURE);
    }
    /* else */
    /* { */
    /*     if (gigabytes > 1) */
    /*         fprintf(stderr, */
    /*                 "\n Successfully re-allocated  %" PRId64 " elements with total size %zu (GB) for variable `%s'
     * \n", */
    /*                 N, gigabytes, varname); */
    /* } */
    return tmp;
}

void my_free(void **x)
{
    if (*x != NULL)
        free(*x);

    *x = NULL;
}

struct io_header get_gadget_header(const char *fname)
{
    FILE *fp = NULL;
    char buf[MAXLEN], buf1[MAXLEN];
    int dummy;
    size_t one = 1;
    struct io_header header;
    my_snprintf(buf, MAXLEN, "%s.%d", fname, 0);
    my_snprintf(buf1, MAXLEN, "%s", fname);

    if (sizeof(struct io_header) != 256)
    {
        fprintf(stderr,
                "ERROR: Gadget header is not %zu bytes and not *exactly* 256 "
                "bytes..exiting\n",
                sizeof(struct io_header));
        exit(EXIT_FAILURE);
    }

    fp = fopen(buf, "r");
    if (fp == NULL)
    {
        fp = fopen(buf1, "r");
        if (fp == NULL)
        {
            fprintf(stderr,
                    "ERROR: Could not find snapshot file.\n neither as `%s' nor as "
                    "`%s'\n",
                    buf, buf1);
            fprintf(stderr, "exiting..\n");
            exit(EXIT_FAILURE);
        }
    }

    //// Don't really care which file actually succeeded (as long as one, buf or
    /// buf1, is present)
    my_fread(&dummy, sizeof(dummy), one, fp);
    my_fread(&header, sizeof(header), one, fp);
    my_fread(&dummy, sizeof(dummy), one, fp);
    fclose(fp);
    return header;
}

int get_gadget_nfiles(const char *fname)
{
    FILE *fp = NULL;
    char buf[MAXLEN], buf1[MAXLEN];
    int dummy;
    struct io_header header;
    my_snprintf(buf, MAXLEN, "%s.%d", fname, 0);
    my_snprintf(buf1, MAXLEN, "%s", fname);

    fp = fopen(buf, "r");
    if (fp == NULL)
    {
        fp = fopen(buf1, "r");
        if (fp == NULL)
        {
            fprintf(stderr,
                    "ERROR: Could not find snapshot file.\n neither as `%s' nor as "
                    "`%s'\n",
                    buf, buf1);
            fprintf(stderr, "exiting..\n");
            exit(EXIT_FAILURE);
        }
    }

    //// Don't really care which file actually succeeded (as long as one, buf or
    /// buf1, is present)
    fread(&dummy, sizeof(dummy), 1, fp);
    fread(&header, sizeof(header), 1, fp);
    fread(&dummy, sizeof(dummy), 1, fp);
    fclose(fp);
    return header.num_files;
}

int64 get_Numpart(struct io_header *header)
{
    int64 N = 0;
    if (header->num_files <= 1)
        for (int i = 0; i < 6; i++)
            header->npartTotal[i] = header->npart[i];

    for (int i = 0; i < 6; i++)
    {
        N += header->npartTotal[i];
        N += (((int64_t)header->npartTotalHighWord[i]) << 32);
    }

    return N;
}

void free_group(struct group_data *g, int64 N)
{
    if (g == NULL)
        return;

    for (int64 i = 0; i < N; i++)
    {
#ifndef MAKE_LEAN
        my_free((void **)&(g[i].x));
        my_free((void **)&(g[i].y));
        my_free((void **)&(g[i].z));
#endif
        my_free((void **)&(g[i].id));
        //   my_free((void **) &(g[i].type));
        my_free((void **)&(g[i].parentgroupforparticle));
        my_free((void **)&(g[i].parentsnapshotforparticle));

#ifdef SUSSING_TREES
        my_free((void **)&(g[i].ParticleEnergy));
        my_free((void **)&(g[i].vx));
        my_free((void **)&(g[i].vy));
        my_free((void **)&(g[i].vz));
#endif
    }
    my_free((void **)&(g));
}

void free_group_positions(struct group_data *g, int64 N)
{
    if (g == NULL)
        return;

    for (int64 i = 0; i < N; i++)
    {
        my_free((void **)&(g[i].x));
        my_free((void **)&(g[i].y));
        my_free((void **)&(g[i].z));
    }
}

struct group_data *allocate_group(int64 N)
{
    struct group_data *g = NULL;
    g = (struct group_data *)malloc(N * sizeof(struct group_data));
    if (g == NULL)
    {
        fprintf(stderr, "Error: Could not allocate memory for %" STR_FMT " elements for current snapshot \n", N);
        exit(EXIT_FAILURE);
    }
    else
        fprintf(stderr, "\n Allocated %" STR_FMT " elements for group struct\n", N);

    // initialize the group
    group_init(g, N);

    return g;
}

void group_init(struct group_data *g, int64 N)
{
    for (int64 i = 0; i < N; i++)
    {
        g[i].ParentId = -1;
        g[i].ParentSnapshot = -1;
        g[i].Ncommon = 0;
        g[i].Rank = 0.0;
        g[i].NParents = 0;
        g[i].NpartinParent = 0;
        g[i].N_per_wedge = 0;
    }
}

float periodic(float dx)
{
    if (dx > 0.5 * PARAMS.BOXSIZE)
        dx -= PARAMS.BOXSIZE;

    if (dx < -0.5 * PARAMS.BOXSIZE)
        dx += PARAMS.BOXSIZE;

    return dx;
}

float periodic_wrap(float x)
{
    while (x > PARAMS.BOXSIZE)
        x -= PARAMS.BOXSIZE;

    while (x < 0)
        x += PARAMS.BOXSIZE;

    return x;
}

double compute_rank(int64 location)
{
    double temp;
    temp = pow((double)(location + 1), (double)(-2. / 3.));

    return temp;
}

void init_all_ranks(double *rank, int64 *ncommon, int64 N)
{
    for (int64 i = 0; i < N; i++)
    {
        rank[i] = 0.0;
        ncommon[i] = 0;
    }
}

int64 find_max_rank(double *NextAllRanks, int64 NextNsub)
{
    double max_rank = 0.0;
    int64 max_rankid = -1; /* Note it is a  signed int/long long*/

    for (int64 k = 0; k < NextNsub; k++)
    {
        if (NextAllRanks[k] > max_rank)
        {
            max_rankid = k;
            max_rank = NextAllRanks[k];
        }
    }

    return max_rankid;
}

int64 find_max_ncommon(int64 *Ncommon, int64 N)
{
    int64 maxncommon = 0, max_ranknum = -1;
    for (int64 k = 0; k < N; k++)
        if (Ncommon[k] > maxncommon)
        {
            maxncommon = Ncommon[k];
            max_ranknum = k;
        }

    return max_ranknum;
}

void reset_ncommon(int64 *Ncommon, int64 N)
{
    // or you could use memset
    for (int64 i = 0; i < N; i++)
        Ncommon[i] = 0;
}

int64 get_ncommon(struct group_data *prev, struct group_data *next)
{

    short *PrevAllPartIds = NULL;
    int64 ncommon = 0;
    id64 PrevMaxPartId;

    assert(prev != NULL && "prev group pointer can not be NULL");
    assert(next != NULL && "next group pointer can not be NULL");

    PrevMaxPartId = 0;
    for (int64 i = 0; i < prev->N; i++)
        if (prev->id[i] > PrevMaxPartId)
            PrevMaxPartId = prev->id[i];

    PrevMaxPartId++;

    PrevAllPartIds = my_calloc(sizeof(*PrevAllPartIds), PrevMaxPartId); /* Note use of calloc instead of malloc */
    for (int64 i = 0; i < prev->N; i++)
    {
        const id64 this_id = prev->id[i];
        assert(this_id < PrevMaxPartId && "Particle ID must be smaller than max. particle ID");
        PrevAllPartIds[this_id] = 1;
    }

    for (int64 i = 0; i < next->N; i++)
    {
        const id64 this_id = next->id[i];
        if (this_id >= PrevMaxPartId)
            continue;

        if (PrevAllPartIds[this_id] == 1)
            ncommon++;
    }
    free(PrevAllPartIds);

    return ncommon;
}

void remove_particle_from_group(const int64 group1, const int64 group2, const int64 part1, const int64 part2,
                                struct group_data *g, int64 *group_to_remove, int64 *part_to_remove)
{
    if (g[group1].fofID == g[group2].fofID)
    {
        fprintf(stderr, "ERROR: Found a duplicate particle in two halos that have the same fofID\n");
        fprintf(stderr, "fofid: %lld, Group1: %lld, Group2: %lld, Part1: %lld, Part2: %lld\n",
                (long long)g[group1].fofID, (long long)group1, (long long)group2, (long long)part1, (long long)part2);
        exit(EXIT_FAILURE);
    }

    const double dx1 = periodic(g[group1].x[part1] - g[group1].xcen);
    const double dy1 = periodic(g[group1].y[part1] - g[group1].ycen);
    const double dz1 = periodic(g[group1].z[part1] - g[group1].zcen);
    const double dist_from_cen1 = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;

    const double dx2 = periodic(g[group2].x[part2] - g[group2].xcen);
    const double dy2 = periodic(g[group2].y[part2] - g[group2].ycen);
    const double dz2 = periodic(g[group2].z[part2] - g[group2].zcen);
    const double dist_from_cen2 = dx2 * dx2 + dy2 * dy2 + dz2 * dz2;

    // Keep the particle in the halo with the smaller distance from the center
    *group_to_remove = dist_from_cen1 < dist_from_cen2 ? group2 : group1;
    *part_to_remove = dist_from_cen1 < dist_from_cen2 ? part2 : part1;
    g[*group_to_remove].id[*part_to_remove] = -1;
    return;
}

int64 remove_duplicates(struct group_data *g, int64 N)
{
    int64_t totnpart = 0;
    id64 max_id = 0;
    for (int64 i = 0; i < N; i++)
    {
        totnpart += g[i].N;
        for (int64 j = 0; j < g[i].N; j++)
            max_id = g[i].id[j] > max_id ? g[i].id[j] : max_id;
    }
    int64_t *all_id_offset = my_calloc(sizeof(*all_id_offset), max_id + 1);
    for (int64 i = 0; i <= max_id + 1; i++)
    {
        all_id_offset[i] = -1;
    }
    int64 *groupnum = my_malloc(sizeof(*groupnum), totnpart);
    int64 *partindex = my_malloc(sizeof(*partindex), totnpart);
    int64_t *num_removed_per_group = my_calloc(sizeof(*num_removed_per_group), N);
    int64 nremoved = 0;
    int64_t offset = 0;
    int interrupted = 0;
    fprintf(stderr, "Removing duplicate particles ...\n");
    init_my_progressbar(N, &interrupted);
    for (int64 i = 0; i < N; i++)
    {
        my_progressbar(i, &interrupted);
        for (int64 j = 0; j < g[i].N; j++)
        {
            id64 id = g[i].id[j];
            if (id > max_id)
            {
                fprintf(stderr, "ERROR: Particle ID %lld is greater than the maximum particle ID %lld\n", (long long)id,
                        (long long)max_id);
                exit(EXIT_FAILURE);
            }
            if (all_id_offset[id] == -1)
            {
                groupnum[offset] = i;
                partindex[offset] = j;
                all_id_offset[id] = offset;
            }
            else
            {
                // fprintf(stderr, "Found a duplicate with id = %lld in group %lld\n", (long long)id, (long long)i);
                int64 group_to_remove, part_to_remove;
                int64_t prev_offset = all_id_offset[id];
                remove_particle_from_group(groupnum[prev_offset], i, partindex[prev_offset], j, g, &group_to_remove,
                                           &part_to_remove);
                num_removed_per_group[group_to_remove]++;

                // If we are keeping the i'th groups particle, then we need to update the groupnum and partindex
                if (group_to_remove != i)
                {
                    groupnum[offset] = i;
                    partindex[offset] = j;
                    all_id_offset[id] = offset;
                }
                nremoved++;
            }
            offset++;
        }
    }
    finish_myprogressbar(&interrupted);
    fprintf(stderr, "Removing duplicate particles ...done. Removed %lld particles ...\n", (long long)nremoved);
    free(all_id_offset);
    free(groupnum);
    free(partindex);

    fprintf(stderr, "Now fixing group particle counts ...\n");
    init_my_progressbar(N, &interrupted);
    for (int64 i = 0; i < N; i++)
    {
        my_progressbar(i, &interrupted);
        if (num_removed_per_group[i] == 0)
            continue;

#define MULTIPLE_ARRAY_EXCHANGER(vartype, a, i, j)                                                                     \
    {                                                                                                                  \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(id64, thisgroup->id, i, j);                                                     \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->x, i, j);                                                     \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->y, i, j);                                                     \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->z, i, j);                                                     \
    }
        struct group_data *thisgroup = &g[i];
        SGLIB_ARRAY_QUICK_SORT(id64, g[i].id, g[i].N, SGLIB_NUMERIC_COMPARATOR, MULTIPLE_ARRAY_EXCHANGER);
#undef MULTIPLE_ARRAY_EXCHANGER

        const int64 start = 0, end = num_removed_per_group[i] - 1;
        XASSERT(num_removed_per_group[i] <= g[i].N,
                "Error: Can remove at most %lld particles (i.e., all particles) from groupnum = %lld. Instead found = "
                "%lld\n",
                (long long)g[i].N, (long long)i, (long long)num_removed_per_group[i]);
        XASSERT(thisgroup->id[start] == -1, "Error: First particle ID should be -1. Instead found = %lld\n",
                (long long)thisgroup->id[start]);
        XASSERT(thisgroup->id[end] == -1, "Error: Last particle ID should  be -1. Instead found = %lld\n",
                (long long)thisgroup->id[end]);
        const int64_t nmove = g[i].N - num_removed_per_group[i];
        // Do a memmove to preserve the ordering of the particles
        memmove(&g[i].id[start], &g[i].id[end + 1], nmove * sizeof(g[i].id[0]));
        memmove(&g[i].x[start], &g[i].x[end + 1], nmove * sizeof(g[i].x[0]));
        memmove(&g[i].y[start], &g[i].y[end + 1], nmove * sizeof(g[i].y[0]));
        memmove(&g[i].z[start], &g[i].z[end + 1], nmove * sizeof(g[i].z[0]));

        g[i].N -= num_removed_per_group[i];

        float *sqr_radius = my_malloc(sizeof(*sqr_radius), g[i].N);
        for (int64 j = 0; j < g[i].N; j++)
        {
            const double dx = periodic(g[i].x[j] - g[i].xcen);
            const double dy = periodic(g[i].y[j] - g[i].ycen);
            const double dz = periodic(g[i].z[j] - g[i].zcen);
            sqr_radius[j] = dx * dx + dy * dy + dz * dz;
        }
#define MULTIPLE_ARRAY_EXCHANGER(vartype, a, i, j)                                                                     \
    {                                                                                                                  \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(id64, thisgroup->id, i, j);                                                     \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->x, i, j);                                                     \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->y, i, j);                                                     \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->z, i, j);                                                     \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, sqr_radius, i, j);                                                       \
    }

        SGLIB_ARRAY_HEAP_SORT(float, sqr_radius, g[i].N, SGLIB_NUMERIC_COMPARATOR, MULTIPLE_ARRAY_EXCHANGER);
        free(sqr_radius);

        // for (int64 j = 0; j < g[i].N; j++)
        // {
        //     if(g[i].id[j] != -1) continue;

        //     /* This particle needs to be removed */
        //     if (j == (g[i].N - 1))
        //     {
        //         g[i].N--;
        //         break;
        //     }

        //     const int64_t nmove = g[i].N - (j + 1);
        //     // Do a memmove to preserve the ordering of the particles
        //     memmove(&g[i].id[j], &g[i].id[j + 1], nmove * sizeof(g[i].id[0]));
        //     memmove(&g[i].x[j], &g[i].x[j + 1], nmove * sizeof(g[i].x[0]));
        //     memmove(&g[i].y[j], &g[i].y[j + 1], nmove * sizeof(g[i].y[0]));
        //     memmove(&g[i].z[j], &g[i].z[j + 1], nmove * sizeof(g[i].z[0]));

        //     g[i].N--;
        //     j--;
        // }

        // if (num_removed_per_group[i] > 0)
        // {
        //     fprintf(stderr, "\nRemoved %lld particles from group %lld\n", (long long)num_removed_per_group[i],
        //             (long long)i);
        //     interrupted = 1;
        // }
    }
    finish_myprogressbar(&interrupted);
    fprintf(stderr, "Now fixing group particle counts ...done\n");

    free(num_removed_per_group);
    return nremoved;
}

short float_almost_equal(float A, float B, int maxUlps)
{

    /* MS -- taken from
               http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
    */

    /* Make sure maxUlps is non-negative and small enough that the
             default NAN won't compare as equal to anything.*/

    assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);

    int aInt = *(int *)&A;

    /* Make aInt lexicographically ordered as a twos-complement int*/

    if (aInt < 0)
        aInt = 0x80000000 - aInt;

    /* Make bInt lexicographically ordered as a twos-complement int*/

    int bInt = *(int *)&B;
    if (bInt < 0)
        bInt = 0x80000000 - bInt;

    int intDiff = abs(aInt - bInt);
    if (intDiff <= maxUlps)
        return 1;

    return 0;
}
