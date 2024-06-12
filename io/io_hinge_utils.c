#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

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
    // VZc                  Rvir                  partid_start_bytes_offset                  partid_end_bytes_offset"

    char wanted_columns[][MAXLEN] = {"ID", "fofid", "nsub", "Mvir", "npart", "Xc",
                                     "Yc", "Zc",    "VXc",  "VYc",  "VZc",   "Rvir"};
    int num_columns = sizeof(wanted_columns) / sizeof(wanted_columns[0]);
    int wanted_columns_indices[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    int num_column_indices = sizeof(wanted_columns_indices) / sizeof(wanted_columns_indices[0]);
    XASSERT(num_columns == num_column_indices, "Number of column names = %d does not match number of indices = %d",
            num_columns, num_column_indices);

    // Now parse the header line
    //  size_t buflen = strnlen(buffer, BUFSIZ);
    char *buf = &buffer[0];

    int colnum = 0;
    char *delim = " #\t\n";
    char *token = NULL;
    while ((token = strsep(&buf, delim)) != NULL)
    {
        if (*token == '\0')
        {
            continue; // Skip empty tokens if multiple delimiters are encountered consecutively
        }
        fprintf(stderr, "column name: '%s' (column number = %d)\n", token, colnum);
        if (colnum >= num_columns)
        {
            fprintf(stderr, "Column number %d is greater than the number of columns %d...breaking\n", colnum,
                    num_columns);
            break;
        }
        else
        {
            fprintf(stderr, "token = '%s'\n", token);
            fprintf(stderr, "colnum = %d wanted_columns[colnum] = '%s'\n", colnum, wanted_columns[colnum]);
            XASSERT(strcasecmp(token, wanted_columns[colnum]) == 0,
                    "Column name '%s' does not match wanted column name '%s'", token, wanted_columns[colnum]);
        }

        colnum++;
    }
    const off_t offset = ftell(fp);
    fprintf(stderr, "offset = %lld\n", offset);
    XASSERT(fgets(buffer, BUFSIZ, fp) != NULL, "Could not read the first data line from file '%s'\n", fname);
    const size_t buflen = strnlen(buffer, BUFSIZ);
    fprintf(stderr, "buffer = '%s'\n", buffer);
    fseek(fp, 0L, SEEK_END);
    const off_t sz = ftell(fp);
    const int64_t nlines = sz / buflen;
    int64_t nallocated = nlines > 1 ? nlines : 1;

    // go back to the beginning of the first data line
    fseek(fp, offset, SEEK_SET);

    struct hinge_catalog *halocat = my_calloc(sizeof(halocat), 1);
    if (halocat == NULL)
    {
        fprintf(stderr, "Could not allocate memory for the halo catalog\n");
        fclose(fp);
        return NULL;
    }

    halocat->halos = my_malloc(sizeof(halocat->halos), nallocated);
    if (halocat->halos == NULL)
    {
        fprintf(stderr, "Could not allocate memory for the halo catalog\n");
        fclose(fp);
        return NULL;
    }
    halocat->nallocated = nallocated;
    halocat->nhalos = 0;
    halocat->nfofs = 0;
    halocat->totnpart = 0;
    int64_t index = 0;
    struct hinge_halo *halos = halocat->halos;
    while (fgets(buffer, BUFSIZ, fp) != NULL)
    {
        if (index == nallocated)
        {
            nallocated *= 1.1;
            halocat->halos = realloc(halocat->halos, nallocated * sizeof(struct hinge_halo));
            if (halocat->halos == NULL)
            {
                fprintf(stderr, "Could not allocate memory for the halo catalog\n");
                fclose(fp);
                return NULL;
            }
            fprintf(stderr, "Reallocating memory for the halo catalog to %" PRId64 "\n", nallocated);
            halos = halocat->halos + index;
        }

        int nread = sscanf(buffer, "%" SCNd64 " %" SCNd64 " %" SCNd64 " %lf %" SCNd64 " %lf %lf %lf %lf %lf %lf %lf",
                           &(halos->halo_id), &halos->fof_id, &halos->nsub, &halos->Mvir, &halos->npart, &halos->Xc,
                           &halos->Yc, &halos->Zc, &halos->VXc, &halos->VYc, &halos->VZc, &halos->Rvir);
        if (nread != num_column_indices)
        {
            fprintf(stderr, "Error: Could not read %d columns from line = '%s'\n", num_column_indices, buffer);
            exit(EXIT_FAILURE);
        }

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
    fclose(fp);

    halocat->nhalos = fof_only ? halocat->nfofs : index;
    halocat->halos = realloc(halocat->halos, halocat->nhalos * sizeof(struct hinge_halo));
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
