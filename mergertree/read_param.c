#include "read_param.h"
#include "defs.h"
#include "utils.h"

void set_simulation_params(struct params_data *params)
{
    params->RedShift = REDSHIFT;
    params->Age = my_malloc(sizeof(*(params->Age)), PARAMS.MAX_SNAPSHOT_NUM + 1);
    for (int i = PARAMS.MIN_SNAPSHOT_NUM; i <= PARAMS.MAX_SNAPSHOT_NUM; i++)
    {
        params->Age[i] = get_age(REDSHIFT[i]); /* Age of the Universe in GYR
                                                  corresponding to the redshifts */
    }
    params->PhysicalLinkLength = params->LINKLENGTH * params->BOXSIZE / pow((double)NUMPART, 1. / 3.);
}

void sanity_check_params(struct params_data *params)
{

    if ((params->MAX_SNAPSHOT_NUM <= params->MIN_SNAPSHOT_NUM) || params->MAX_SNAPSHOT_NUM < 0 ||
        params->MIN_SNAPSHOT_NUM < 0)
    {
        fprintf(stderr,
                "ERROR: [min,max] snapshot numbers = [%d,%d] are not well-behaved "
                "...exiting\n",
                params->MIN_SNAPSHOT_NUM, params->MAX_SNAPSHOT_NUM);
        exit(EXIT_FAILURE);
    }
}

void fill_config_params(struct params_data *params)
{

#ifdef GET_GROUPVEL
    params->get_groupvel = 1;
#else
    params->get_groupvel = 0;
#endif

#ifdef GET_MEANVEL
    params->get_meanvel = 1;
#else
    params->get_meanvel = 1;
#endif

#ifdef BIGSIM
    params->bigsim = 1;
#else
    params->bigsim = 0;
#endif

#ifdef LONGIDS
    params->longids = 1;
#else
    params->longids = 0;
#endif

#ifdef FOF_ONLY
    params->fof_only = 1;
#else
    params->fof_only = 0;
#endif

#ifdef USE_MOST_BOUND_FOR_CENTRE
    params->use_most_bound_for_centre = 1;
#else
    params->use_most_bound_for_centre = 0;
#endif
}

void output_params(const char *fname, struct params_data *params)
{
    FILE *fp = NULL;
    move_existing_file(fname); // mv old param file if it exists
    fp = my_fopen(fname, "w");

    fprintf(fp, "MIN_SNAPSHOT_NUM              %d\n", params->MIN_SNAPSHOT_NUM);
    fprintf(fp, "MAX_SNAPSHOT_NUM              %d\n", params->MAX_SNAPSHOT_NUM);

    fprintf(fp, "SNAPSHOT_DIR                  %s\n", params->SNAPSHOT_DIR);
    fprintf(fp, "SNAPSHOT_BASE                 %s\n", params->SNAPSHOT_BASE);

    fprintf(fp, "GROUP_DIR                  %s\n", params->GROUP_DIR);
    fprintf(fp, "GROUP_BASE                 %s\n", params->GROUP_BASE);

    fprintf(fp, "OUTPUT_DIR                 %s\n", params->OUTPUT_DIR);

    fprintf(fp, "BOXSIZE                    %lf\n", params->BOXSIZE);
    fprintf(fp, "LINKLENGTH                        %lf\n", params->LINKLENGTH);
    fprintf(fp, "INDIVIDUAL_MERGERS                %d\n", params->INDIVIDUAL_MERGERS);

    fprintf(fp, "\n## config options [from Makefile] ##\n");
    fprintf(fp, "## FOF_ONLY                       %d\n", params->fof_only);
    fprintf(fp, "## GET_GROUPVEL                   %d\n", params->get_groupvel);
    fprintf(fp, "## GET_MEANVEL                    %d\n", params->get_meanvel);
    fprintf(fp, "## BIGSIM                         %d\n", params->bigsim);
    fprintf(fp, "## LONGIDS                        %d\n", params->longids);
    fprintf(fp, "## USE_MOST_BOUND_FOR_CENTRE      %d\n", params->use_most_bound_for_centre);

#ifdef SUBFIND
    fprintf(fp, "\n## Parameters read in from Gadget snapshot ##\n");
    for (int i = 0; i < 6; i++)
        fprintf(fp, "## MASSARR[%d]         %e\n", i, params->MASSARR[i]);
#endif

    fclose(fp);
}

void read_params(const char *fname, struct params_data *params)
{

#define DOUBLE 1
#define STRING 2
#define INT 3
#define INT64 4
#define MAXTAGS 20

    FILE *fd = NULL;
    char buf[MAXLINESIZE], buf1[MAXLINESIZE], buf2[MAXLINESIZE], buf3[MAXLINESIZE];
    int i, j, nt = 0;
    int id[MAXTAGS];
    void *addr[MAXTAGS];
    char tag[MAXTAGS][MAXLEN];
    int errorFlag = 0;
    int64 tmp_64_int;

    nt = 0;

    my_snprintf(tag[nt], MAXLEN, "MIN_SNAPSHOT_NUM");
    addr[nt] = &(params->MIN_SNAPSHOT_NUM);
    id[nt++] = INT;

    my_snprintf(tag[nt], MAXLEN, "MAX_SNAPSHOT_NUM");
    addr[nt] = &(params->MAX_SNAPSHOT_NUM);
    id[nt++] = INT;

    my_snprintf(tag[nt], MAXLEN, "SNAPSHOT_DIR");
    addr[nt] = params->SNAPSHOT_DIR;
    id[nt++] = STRING;

    my_snprintf(tag[nt], MAXLEN, "SNAPSHOT_BASE");
    addr[nt] = params->SNAPSHOT_BASE;
    id[nt++] = STRING;

    my_snprintf(tag[nt], MAXLEN, "GROUP_DIR");
    addr[nt] = params->GROUP_DIR;
    id[nt++] = STRING;

    my_snprintf(tag[nt], MAXLEN, "GROUP_BASE");
    addr[nt] = params->GROUP_BASE;
    id[nt++] = STRING;

    my_snprintf(tag[nt], MAXLEN, "OUTPUT_DIR");
    addr[nt] = params->OUTPUT_DIR;
    id[nt++] = STRING;

    my_snprintf(tag[nt], MAXLEN, "BOXSIZE");
    addr[nt] = &(params->BOXSIZE);
    id[nt++] = DOUBLE;

    my_snprintf(tag[nt], MAXLEN, "LINKLENGTH");
    addr[nt] = &(params->LINKLENGTH);
    id[nt++] = DOUBLE;

    my_snprintf(tag[nt], MAXLEN, "INDIVIDUAL_MERGERS");
    addr[nt] = &(params->INDIVIDUAL_MERGERS);
    id[nt++] = INT;

    fd = my_fopen(fname, "r");
    while (!feof(fd))
    {
        *buf = 0;
        fgets(buf, MAXLEN, fd);
        if (sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
            continue;

        if (buf1[0] == '%' || buf1[0] == '#')
            continue;

        for (i = 0, j = -1; i < nt; i++)
        {
            if (strcmp(buf1, tag[i]) == 0)
            {
                j = i;
                tag[i][0] = '0';
                break;
            }
        }

        if (j >= 0)
        {
            switch (id[j])
            {
            case DOUBLE:
                *((double *)addr[j]) = atof(buf2);
                break;
            case STRING:
                strcpy(addr[j], buf2);
                break;
            case INT:
                *((int *)addr[j]) = atoi(buf2);
                break;

            case INT64:
                sscanf(buf2, "%" STR_FMT, &tmp_64_int);
                *((int64 *)addr[j]) = tmp_64_int;
                break;
            }
        }
        else
        {
            fprintf(stderr, "Error in file %s:   Tag '%s' not allowed or multiply defined.\n", fname, buf1);
            errorFlag = 1;
        }
    }
    fclose(fd);

    if (errorFlag != 2)
    {
        for (i = 0; i < nt; i++)
        {
            if (tag[i][0] != '0')
            {
                fprintf(stderr, "Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
                errorFlag = 1;
            }
        }
    }

    if (errorFlag != 0)
    {
        fprintf(stderr, "Error: Parameter file not as expected..bailing\n");
        exit(EXIT_FAILURE);
    }

#undef DOUBLE
#undef STRING
#undef INT
#undef LONG
#undef MAXTAGS
}
