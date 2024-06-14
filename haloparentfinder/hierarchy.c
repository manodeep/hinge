#include "hierarchy.h"
#include "progressbar.h"
#include "utils.h"

/*
This function will match the sub-subhalos (and under) to their
container halos.

*/
// functions that are not exposed to the outside
static inline int find_octant(float *pos, float *origin);
void print_subhalolevel_header(FILE *fp);
void find_hierarchy_in_fof(struct group_data *group, const int64 StartFofId);

static inline void rotate_prbar(int PLEN)
{
    static short index = 0;
    char a = '\0';
    switch (index)
    {
    case 0:
        a = '\\';
        index++;
        break;
    case 1:
        a = '|';
        index++;
        break;
    case 2:
        a = '/';
        index++;
        break;
    case 3:
        a = '-';
        index = 0;
        break;
    }
    for (int i = 0; i < PLEN; i++)
        fprintf(stderr, "%c", a);
}

static inline void print_char(const char s, const int LEN)
{
    for (int i = 0; i < LEN; i++)
    {
        fprintf(stderr, "%c", (char)s);
    }
}

static inline void estimate_eta(const time_t t0, const float percent, int *hr, int *min, int *sec)
{
    float left = 100.0 - percent; /* percent is percent_done*/
    double ratios[] = {3600.0, 60.0, 1};
    int temp;
    double timeleft = (double)t0 * left / percent;

    if (timeleft < ratios[2])
        timeleft = (double)t0; /* display elapsed time */

    temp = floor(timeleft / ratios[0]);
    *hr = temp > 0 ? temp : 0;
    temp = floor((timeleft - (*hr) * ratios[0]) / ratios[1]);
    *min = temp > 0 ? temp : 0;
    temp = floor(timeleft - (*hr) * ratios[0] - (*min) * ratios[1]);
    *sec = temp > 0 ? temp : 0;
}

static inline int find_octant(float *pos, float *origin)
{
    int octant = 0;
    /*   x > 0 -> 0,3,4,7;   x < 0 -> 1,2,5,6 ; */
    /*   y > 0 -> 0,1,4,5;   y < 0 -> 2,3,6,7 ; */
    /*   z > 0 -> 0,1,2,3;   z < 0 -> 4,5,6,7 ; */

    for (int i = 0; i < 3; i++)
    {
        if (periodic(pos[i] - origin[i]) < 0.0)
            octant += (int)pow(2.0, i); // really is 2^i but pow accepts floats
    }

    return octant;
}

void find_hierarchy_in_fof(struct group_data *group, const int64 StartFofId)
{
    int64 *smallest_match = NULL;
    int64 Nsub;
    int flag = 1, iter = 0;
    const int MAXITER = 2; /* 1 iteration should be enough since small things go
                              inside bigger things..but just in case*/
    int64 index;
    const int NWEDGES = 8; /*octants. If you want to increase it, make sure that
                              the logic in find_octant works for the new value. */
    int64 Counts[NWEDGES];
    float origin[3], pos[3];
    int octant;
    int flag_no_match = 1;
    /* int PRINTSTEP=0,SMALLPRINTSTEP=0,PRINTLEN = 100; */
    // int SMALLPRINTSTEP = 0, PRINTLEN = 100;
    // float percent_done = 0.0;
    /* char octant_flag[NWEDGES+1],final_flag[NWEDGES+1]; */
    // time_t t_start, time_t t_now, time_taken;
    // int hr, min, sec;
    int N_per_wedge;
    int Min_Nsub_forPBar = 300; /* >= 100 */

    if (group != NULL && group[StartFofId].isFof == 1 && group[StartFofId].Nsub > 2)
    {
        Nsub = group[StartFofId].Nsub;
        smallest_match = my_malloc(sizeof(*smallest_match), Nsub + 1);

        for (int64 k = 0; k <= Nsub; k++)
            smallest_match[k] = -1;

        /* for(int k=0;k<NWEDGES;k++) */
        /* 		final_flag[k]  = '1'; */

        /* final_flag[NWEDGES]  = '\0'; */
        /* octant_flag[NWEDGES] = '\0'; */

        /* Displaying the rotating progress bar */
        int interrupted = 0;
        if (Nsub > Min_Nsub_forPBar)
        {
#if 0
            /* PRINTSTEP = (int) floor(0.1*Nsub); */
            SMALLPRINTSTEP = ceil(0.01 * Nsub) > 1 ? ceil(0.01 * Nsub) : 1;
            fprintf(stderr, "\n\n");
            fprintf(stderr, "\n Working on Fofhalo id %" STR_FMT " with %" STR_FMT " subhalos \n\n", StartFofId, Nsub);
            fprintf(stderr, " %3d%%", 0);
            fprintf(stderr, "\b|");
            print_char(' ', PRINTLEN);
            fprintf(stderr, "\b|");
            fprintf(stderr, "ETA: --:--:--");
#endif
            init_my_progressbar(Nsub, &interrupted);
        }

        // t_start = time(NULL);
        for (int64 igroup = (StartFofId + Nsub - 1); igroup > StartFofId + 1; igroup--)
        {
            N_per_wedge = (int)sqrt((double)group[igroup].N) / NWEDGES;
            N_per_wedge = N_per_wedge > 100 ? 100 : (N_per_wedge < 2 ? 2 : N_per_wedge);

            group[igroup].N_per_wedge = (short)N_per_wedge;

            origin[0] = group[igroup].xcen;
            origin[1] = group[igroup].ycen;
            origin[2] = group[igroup].zcen;
            flag_no_match = 1;
            for (int64 j = igroup - 1; j >= (StartFofId + 1); j--)
            {
                for (int k = 0; k < NWEDGES; k++)
                {
                    /* octant_flag[k] = '0'; */
                    Counts[k] = 0;
                }

                for (int64 k = 0; k < group[j].N; k++)
                {
                    pos[0] = group[j].x[k];
                    pos[1] = group[j].y[k];
                    pos[2] = group[j].z[k];
                    octant = find_octant(pos, origin);
                    if (octant >= 0 && octant < NWEDGES)
                    {
                        Counts[octant]++;
                        /* octant_flag[octant] = '1'; */
                    }
                    else
                    {
                        fprintf(stderr,
                                "Error in calculating octant: I have octant = %d and it "
                                "should have been (0 <= octant < %d) ",
                                octant, NWEDGES);
                        exit(EXIT_FAILURE);
                    }

                    /* This is the more efficient way of finding the container but
                             somehow this does NOT work */

                    /* 				  if(strcmp(octant_flag,final_flag) == 1)
                     */
                    /* 					{ */
                    /* 					  smallest_match[igroup] = j; */
                    /* 					  flag_no_match = 0; */
                    /* 					  break; */
                    /* 					} */
                }

                flag_no_match = 1;
                for (int k = 0; k < NWEDGES; k++)
                {
                    if (Counts[k] < N_per_wedge)
                    {
                        flag_no_match = 1;
                        break;
                    }
                    else
                    {
                        flag_no_match = 0;
                    }
                }

                if (flag_no_match == 0)
                {
                    smallest_match[igroup - StartFofId] = j;
                    flag_no_match = 1;
                    break;
                }
            }

            /* For displaying the rotating bar */
            if (Nsub > Min_Nsub_forPBar)
            {
                my_progressbar(igroup, &interrupted);
#if 0
                if (((igroup - StartFofId) % SMALLPRINTSTEP) == 0)
                {
                    percent_done = 100.0 - (double)(igroup - StartFofId) / Nsub * 100.0;
                    if (percent_done > 99.0)
                        percent_done = 100.0;
                    else if (percent_done < 0.0)
                        percent_done = 0.0;

                    fprintf(stderr, "\r %3d%% ", (int)percent_done);
                    fprintf(stderr, "\b|");
                    rotate_prbar((int)floor(percent_done));
                    print_char(' ', PRINTLEN - (int)percent_done);
                    fprintf(stderr, "\b|");
                    t_now = time(NULL);
                    time_taken = difftime(t_now, t_start);
                    estimate_eta(time_taken, percent_done, &hr, &min, &sec);
                    fprintf(stderr, "ETA:: %02dh:%02dm:%02ds", hr, min, sec);
                }
#endif
            }
        }

        if (Nsub > Min_Nsub_forPBar)
        {
#if 0
            fprintf(stderr, "\r %3d%% ", 100);
            fprintf(stderr, "\b|");
            print_char('|', PRINTLEN);
            fprintf(stderr, "\b|");
            t_now = time(NULL);
            time_taken = difftime(t_now, t_start);
            estimate_eta(time_taken, 100.0, &hr, &min, &sec);
            fprintf(stderr, "Time: %02dh:%02dm:%02ds DONE", hr, min, sec);
#endif
            finish_myprogressbar(&interrupted);
        }

        /* Now figure out the actual value of the levels */
        group[StartFofId + 1].ParentLevel = 2;
        flag = 1;
        iter = 0;
        while (flag == 1 && iter < MAXITER)
        {
            for (int64 igroup = StartFofId + 2; igroup < (StartFofId + Nsub); igroup++)
            {
                /* 			  fprintf(stderr,"igroup = %"STR_FMT"  startfofid =
                 * %"STR_FMT"  Nsub = %"STR_FMT" \n",igroup,StartFofId,Nsub); */
                index = smallest_match[igroup - StartFofId];
                if (index != -1)
                {
                    if (group[index].ParentLevel != -1)
                    { /* Check if the container has ParentLevel already set (this
                         should be the case, almost always) */
                        group[igroup].ParentLevel = group[index].ParentLevel + 1;
                        group[igroup].ContainerIndex = index;
                        if (group[index].ParentLevel > 1)
                            group[index].Nsub++;

                        flag = 0;
                    }
                    else
                    {
                        flag = 1;
                    }
                }
                else
                {
                    group[igroup].ParentLevel = 2;
                    group[igroup].ContainerIndex = StartFofId;
                    flag = 0;
                }
            }
            iter++;
        }
        free(smallest_match);
    }
    else
    {
        fprintf(stderr, "This should not have happened -- trying to match subhalos "
                        "on a NULL or only FOF group\n");
        if (group != NULL)
            fprintf(stderr, "StartFofId = %" STR_FMT " and Nsub = %" STR_FMT " isfof = %3d ..exiting\n", StartFofId,
                    group[StartFofId].Nsub, group[StartFofId].isFof);

        exit(EXIT_FAILURE);
    }
}

void print_subhalolevel_header(FILE *fp)
{

    fprintf(fp, "################################################################"
                "##############\n");
    fprintf(fp, "#   FofNum          GroupNum  ParentLev   ContainerNum     Nsub "
                " N_per_wedge #\n");
    fprintf(fp, "################################################################"
                "##############\n");
}

void find_hierarchy_level(struct group_data *group, const int64 Ngroups, const char *outpath)
{
    int64 FofId = 0;
    FILE *fp = NULL;
    char fname[MAXLEN];
    int64 i;
    char str_line[MAXLINESIZE];
    char comment = '#';

    if (Ngroups > 0)
    {
        my_snprintf(fname, MAXLEN, "%s/subhalolevel_%03d.txt", outpath, group->snapshot);
        fp = fopen(fname, "rt");
        if (fp == NULL)
        {
            fp = my_fopen(fname, "w");
            print_subhalolevel_header(fp);
            while (FofId < Ngroups)
            {
                /*if the Fof has exactly one Nsub, i.e., the Fof halo itself, then
                        we know the ParentLevel and the ContainerIndex (set in
                   loadgroups.c) */
                assert(group[FofId].Nsub >= 1 && "FOF halos must have Nsub >=1 (since they contain themselves)");
                if (group[FofId].Nsub > 1)
                {
                    if (group[FofId].Nsub == 2)
                    {
                        group[FofId + 1].ParentLevel = 2; /* legitimate subhalo */
                        group[FofId + 1].ContainerIndex = FofId;
                    }
                    else
                    {
                        /* Make wedges: but let's start out with octants to begin with*/
                        find_hierarchy_in_fof(group, FofId);
                    }
                } /* else { */
                /* 	//Nsub == 1 */
                /* 	fprintf(stderr,"FofId = %"STR_FMT" Nsub = %"STR_FMT" Parentlevel
                 * = %hd \n",FofId,group[FofId].Nsub,group[FofId].ParentLevel); */
                /* 	fprintf(stderr,"WARNING: Hardwiring values - this should not
                 * happen - check the data load\n"); */
                /* 	group[FofId].ParentLevel = 1; */
                /* 	//Hack -> setting Nsub to 1.  */
                /* 	group[FofId].Nsub = 1; */
                /* 	group[FofId].ContainerIndex = FofId; */
                /* } */

                for (i = FofId; i < (FofId + group[FofId].Nsub); i++)
                    fprintf(fp, "%10" STR_FMT "     %10" STR_FMT "   %6hd    %10" STR_FMT "    %10" STR_FMT "   %5d \n",
                            FofId, i, group[i].ParentLevel, group[i].ContainerIndex, group[i].Nsub,
                            group[i].N_per_wedge);

                FofId += group[FofId].Nsub;
            }

            fclose(fp);
        }
        else
        {
            /* So, the subhalolevel file exists for this snapshots. Let's just
                     read it in.
            */

            fprintf(stderr, "Reading in the parent levels from file `%s' ", fname);
            i = 0;
            while (1)
            {
                if (fgets(str_line, MAXLINESIZE, fp) != NULL)
                {
                    if (str_line[0] != comment)
                    {
                        int64 dummy;
                        /* 					  sscanf(str_line,"%*"STR_FMT" %"STR_FMT" %hd
                         * %"STR_FMT"  %"STR_FMT"     %hd",  */
                        int nread = sscanf(str_line, "%*d   %" STR_FMT " %hd   %" STR_FMT "  %" STR_FMT "     %hd",
                                           &dummy, &group[i].ParentLevel, &group[i].ContainerIndex, &group[i].Nsub,
                                           &group[i].N_per_wedge); /*discarding the first field, Fofid */
                        assert(nread == 5 && "ERROR: Could not read in the 5 values "
                                             "expected from subhalo hierarchy file");
                        if (dummy != i)
                        {
                            fprintf(stderr, "Error: While reading in file `%s' for parentlevel \n ", fname);
                            fprintf(stderr, "expected igroup = %" STR_FMT " instead got %" STR_FMT " ..exiting \n\n", i,
                                    dummy);
                            exit(EXIT_FAILURE);
                        }
                        i++;
                    }
                }
                else
                {
                    break;
                }
            }
            fclose(fp);
            fprintf(stderr, " ..done\n");
        }
    } // Ngroups > 0
}
