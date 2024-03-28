#include "loadsnapshot.h"
#include "sglib.h"
#include "utils.h"

/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */

struct particle_data *loadsnapshot(const char *fname, struct io_header *header)
{
    FILE *fd = NULL;
    char buf[MAXLEN];
    int i, dummy;
    int64 k, ntot_withmasses;
    int64 n, pc, pc_new;
    id64 *Id = NULL;
    size_t nmemb = 0;
    const size_t OneElement = 1;
    int files = get_gadget_nfiles(fname);
    struct particle_data *P = NULL;

#define SKIP my_fread(&dummy, sizeof(dummy), OneElement, fd);

    for (i = 0, pc = 0; i < files; i++, pc = pc_new)
    {
        if (files > 1)
        {
            my_snprintf(buf, MAXLEN, "%s.%d", fname, i);
        }
        else
        {
            my_snprintf(buf, MAXLEN, "%s", fname);
        }

        fd = my_fopen(buf, "r");
        fprintf(stderr, "reading `%s' ...\n", buf);

        nmemb = 1;
        my_fread(&dummy, sizeof(dummy), nmemb, fd);
        my_fread(header, sizeof(struct io_header), nmemb, fd);
        my_fread(&dummy, sizeof(dummy), nmemb, fd);

        for (k = 0, ntot_withmasses = 0; k < 5; k++)
        {
            if (fabs(header->mass[k]) < DOUBLE_EPS)
                ntot_withmasses += header->npart[k];
        }

        if (i == 0)
        {
            P = my_malloc(sizeof(struct particle_data), NUMPART);
            Id = my_malloc(sizeof(id64), NUMPART);
        }

        SKIP;
        for (k = 0, pc_new = pc; k < 6; k++)
        {
            for (n = 0; n < header->npart[k]; n++)
            {
                nmemb = 3;
                my_fread(&P[pc_new].Pos[0], sizeof(float), nmemb, fd);
                pc_new++;
            }
        }
        SKIP;

        SKIP;
        for (k = 0, pc_new = pc; k < 6; k++)
        {
            for (n = 0; n < header->npart[k]; n++)
            {
                nmemb = 3;
                my_fread(&P[pc_new].Vel[0], sizeof(float), nmemb, fd);
                pc_new++;
            }
        }
        SKIP;

        SKIP;
        for (k = 0, pc_new = pc; k < 6; k++)
        {
            for (n = 0; n < header->npart[k]; n++)
            {
                my_fread(&Id[pc_new], sizeof(id64), OneElement, fd);
                pc_new++;
            }
        }
        SKIP;

        if (ntot_withmasses > 0)
            SKIP;
        for (k = 0, pc_new = pc; k < 6; k++)
        {
            for (n = 0; n < header->npart[k]; n++)
            {
                P[pc_new].Type = k;

                if (fabs(header->mass[k]) < DOUBLE_EPS)
                {
                    my_fread(&P[pc_new].Mass, sizeof(float), OneElement, fd);
                }
                else
                {
                    P[pc_new].Mass = header->mass[k];
                }
                pc_new++;
            }
        }
        if (ntot_withmasses > 0)
            SKIP;

        fclose(fd);
    }
    reordering(P, Id, NUMPART); // so I can access by particle ids in assign_vxcm
    free(Id);
    return P;
}

void reordering(struct particle_data *P, id64 *Id, int64 N)
{
    fprintf(stderr, "reordering....");
#define MULTIPLE_ARRAY_EXCHANGER(type, a, i, j)                                                                        \
    {                                                                                                                  \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(struct particle_data, P, i, j);                                                 \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(id64, Id, i, j);                                                                \
    }
    SGLIB_ARRAY_QUICK_SORT(id64, Id, N, SGLIB_SAFE_NUMERIC_COMPARATOR, MULTIPLE_ARRAY_EXCHANGER);

    fprintf(stderr, "done.\n");
}
