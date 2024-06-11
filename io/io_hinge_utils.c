#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "io_hinge_utils.h"
#include "macros.h"
#include "utils.h"

struct hinge_catalog *read_hinge_ascii_halo_catalog(const char *fname, const int fof_only)
{
    FILE *fp = my_fopen(fname, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Could not open file %s\n", fname);
        return NULL;
    }
    char buffer[BUFSIZ];
    if (fgets(buffer, BUFSIZ, fp) == NULL)
    {
        fprintf(stderr, "Could not read file %s\n", fname);
        fclose(fp);
        return NULL;
    }
    //"# ID                  fofid                  nsub                  Mvir                  npart Xc Yc Zc VXc VYc
    //VZc                  Rvir                  partid_start_bytes_offset                  partid_end_bytes_offset"

    char wanted_columns[][MAXLEN] = {"ID", "fofid", "nsub", "Mvir", "npart", "Xc",
                                     "Yc", "Zc",    "VXc",  "VYc",  "VZc",   "Rvir"};
    int num_columns = sizeof(wanted_columns) / sizeof(wanted_columns[0]);
    int wanted_columns_indices[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int num_column_indices = sizeof(wanted_columns_indices) / sizeof(wanted_columns_indices[0]);
    XASSERT(num_columns == num_column_indices, "Number of column names = %d does not match number of indices = %d",
            num_columns, num_column_indices);

    // Now parse the header line
    //  size_t buflen = strnlen(buffer, BUFSIZ);
    char *token = &buffer[0];
    while (*token == '#')
        token++;

    int colnum = 0;
    while ((token = strsep((char **)&buffer, " ")) != NULL)
    {
        fprintf(stderr, "column name: '%s' (column number = %d)\n", token, colnum);
        XASSERT(colnum < num_columns, "Column number %d is greater than the number of columns %d", colnum, num_columns);
        XASSERT(strcasecmp(token, wanted_columns[colnum]) == 0,
                "Column name '%s' does not match wanted column name '%s'", token, wanted_columns[colnum]);
        colnum++;
    }
    int64 nhalos = getnumlines(fname, '#');
    struct hinge_catalog *halocat = my_calloc(sizeof(halocat), 1);
    if (halocat == NULL)
    {
        fprintf(stderr, "Could not allocate memory for the halo catalog\n");
        fclose(fp);
        return NULL;
    }

    halocat->halos = my_malloc(sizeof(halocat->halos), nhalos);
    if (halocat->halos == NULL)
    {
        fprintf(stderr, "Could not allocate memory for the halo catalog\n");
        fclose(fp);
        return NULL;
    }
    halocat->nallocated = nhalos;
    halocat->nhalos = nhalos;
    halocat->nfofs = 0;
    halocat->totnpart = 0;
    int64_t index = 0;
    struct hinge_halo *halos = halocat->halos;
    while (fscanf(fp, "%" SCNd64 " %" SCNd64 " %" SCNd64 " %lf %" SCNd64 " %lf %lf %lf %lf %lf %lf %lf",
                  &(halos->halo_id), &halos->fof_id, &halos->nsub, &halos->Mvir, &halos->npart, &halos->Xc, &halos->Yc,
                  &halos->Zc, &halos->VXc, &halos->VYc, &halos->VZc, &halos->Rvir) == num_columns)
    {
        if (fof_only && halos->fof_id != halos->halo_id)
        {
            continue;
        }

        if (halos->fof_id == halos->halo_id)
        {
            halocat->nfofs++;
        }
        halocat->totnpart += halos->npart;

        halos++;
        index++;
    }
    if (fof_only)
    {
        halocat->nhalos = halocat->nfofs;
        halocat->halos = realloc(halocat->halos, halocat->nfofs * sizeof(struct hinge_halo));
    }
    else
    {
        XASSERT(index == nhalos, "Read %lld halos but expected %lld", (long long)index, (long long)nhalos);
    }

    fclose(fp);
    return halocat;
}

void free_hinge_halocat(struct hinge_catalog *halocat)
{
    if (halocat == NULL)
    {
        return;
    }
    if (halocat->halos != NULL)
    {
        free(halocat->halos);
    }
    free(halocat);
    return;
}