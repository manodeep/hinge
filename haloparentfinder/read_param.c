#include "read_param.h"
#include "defs.h"
#include "utils.h"

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

#ifdef MAKE_LEAN
    params->make_lean = 1;
#else
    params->make_lean = 0;
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
    fprintf(fp, "GROUP_FORMAT               %s\n", GROUP_FORMAT_NAMES[params->GROUP_FORMAT]);

    fprintf(fp, "OUTPUT_DIR                 %s\n", params->OUTPUT_DIR);

    fprintf(fp, "MAX_INCR                          %d\n", params->MAX_INCR);
    fprintf(fp, "MAX_RANK_LOC                      %" STR_FMT "\n", params->MAX_RANK_LOC);

    fprintf(fp, "MIN_FCOMMON_FINDPROGENITOR_THRESH                       %lf\n",
            params->MIN_FCOMMON_FINDPROGENITOR_THRESH);
    fprintf(fp, "MIN_NUMPART_IN_FINDPROGENITOR_HALO                      %" STR_FMT "\n",
            params->MIN_NUMPART_IN_FINDPROGENITOR_HALO);

    fprintf(fp, "MIN_FCOMMON_SWITCHFOF_THRESH                           %lf\n", params->MIN_FCOMMON_SWITCHFOF_THRESH);
    fprintf(fp, "MIN_NUMPART_IN_SWITCHFOF_HALO                          %" STR_FMT "\n",
            params->MIN_NUMPART_IN_SWITCHFOF_HALO);

    fprintf(fp, "\n## config options [from Makefile] ##\n");
    fprintf(fp, "## FOF_ONLY              %d\n", params->fof_only);
    fprintf(fp, "## GET_GROUPVEL          %d\n", params->get_groupvel);
    fprintf(fp, "## BIGSIM                %d\n", params->bigsim);
    fprintf(fp, "## LONGIDS               %d\n", params->longids);
    fprintf(fp, "## MAKE_LEAN             %d\n", params->make_lean);

    fprintf(fp, "\n## Parameters read in from Gadget snapshot ##\n");
    fprintf(fp, "## BOXSIZE               %lf\n", params->BOXSIZE);
    for (int i = 0; i < 6; i++)
        fprintf(fp, "## MASSARR[%d]         %e\n", i, params->MASSARR[i]);
    fclose(fp);
}

/*
        The following function is shamelessly taken from Gadget2 */
void haloparentfinder_fill_params(struct params_data *params, const int maxtags, void **addr, int *id, char **tag,
                                  int *nt_out)
{
    int nt = *nt_out;
    my_snprintf(tag[nt], MAXLEN, "MAX_INCR");
    addr[nt] = &(params->MAX_INCR);
    id[nt++] = INT;
    assert(nt < maxtags);

    my_snprintf(tag[nt], MAXLEN, "MAX_RANK_LOC");
    addr[nt] = &(params->MAX_RANK_LOC);
    id[nt++] = INT64;
    assert(nt < maxtags);

    my_snprintf(tag[nt], MAXLEN, "MIN_FCOMMON_FINDPROGENITOR_THRESH");
    addr[nt] = &(params->MIN_FCOMMON_FINDPROGENITOR_THRESH);
    id[nt++] = DOUBLE;
    assert(nt < maxtags);

    my_snprintf(tag[nt], MAXLEN, "MIN_NUMPART_IN_FINDPROGENITOR_HALO");
    addr[nt] = &(params->MIN_NUMPART_IN_FINDPROGENITOR_HALO);
    id[nt++] = INT64;
    assert(nt < maxtags);

    my_snprintf(tag[nt], MAXLEN, "MIN_FCOMMON_SWITCHFOF_THRESH");
    addr[nt] = &(params->MIN_FCOMMON_SWITCHFOF_THRESH);
    id[nt++] = DOUBLE;
    assert(nt < maxtags);

    my_snprintf(tag[nt], MAXLEN, "MIN_NUMPART_IN_SWITCHFOF_HALO");
    addr[nt] = &(params->MIN_NUMPART_IN_SWITCHFOF_HALO);
    id[nt++] = INT64;
    assert(nt < maxtags);

    *nt_out = nt;
}

void read_params(const char *fname, struct params_data *params,
                 void special_params(struct params_data *params, const int maxtags, void **addr, int *id, char **tag,
                                     int *nt))
{

    FILE *fd = NULL;
    char buf[MAXLINESIZE], buf1[MAXLINESIZE], buf2[MAXLINESIZE], buf3[MAXLINESIZE];
    int i, j, nt = 0;
    int id[MAXTAGS];
    void *addr[MAXTAGS];
    char tag[MAXTAGS][MAXLEN];
    int errorFlag = 0;
    int64 tmp_64_int;
    char group_format_type[MAXLEN];

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

    my_snprintf(tag[nt], MAXLEN, "GROUP_FORMAT");
    addr[nt] = group_format_type;
    id[nt++] = STRING;

    my_snprintf(tag[nt], MAXLEN, "OUTPUT_DIR");
    addr[nt] = params->OUTPUT_DIR;
    id[nt++] = STRING;

    special_params(params, MAXTAGS, addr, id, (char **)&tag, &nt);

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
            fprintf(stderr, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n", fname, buf1);
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

#define CHECK_VALID_ENUM_IN_PARAM_FILE(paramname, num_enum_types, enum_names, enum_values, string_value)               \
    {                                                                                                                  \
        int found = 0;                                                                                                 \
        for (int ii = 0; ii < num_enum_types; ii++)                                                                    \
        {                                                                                                              \
            if (strcasecmp(string_value, enum_names[ii]) == 0)                                                         \
            {                                                                                                          \
                params->paramname = enum_values[ii];                                                                   \
                found = 1;                                                                                             \
                break;                                                                                                 \
            }                                                                                                          \
        }                                                                                                              \
        if (found == 0)                                                                                                \
        {                                                                                                              \
            fprintf(stderr, #paramname " field contains unsupported value of '%s' is not supported\n", string_value);  \
            fprintf(stderr, " Please choose one of the values -- \n");                                                 \
            for (int ii = 0; ii < num_enum_types; ii++)                                                                \
            {                                                                                                          \
                fprintf(stderr, #paramname " = '%s'\n", enum_names[ii]);                                               \
            }                                                                                                          \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    }

    CHECK_VALID_ENUM_IN_PARAM_FILE(GROUP_FORMAT, nvalid_group_format_names, GROUP_FORMAT_NAMES, GROUP_FORMAT_ENUMS,
                                   group_format_type);

    sanity_check_params(params);
    fill_config_params(params);

#undef DOUBLE
#undef STRING
#undef INT
#undef LONG
#undef MAXTAGS

    return;
}
