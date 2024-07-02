#include "loadparents.h"
#include "hinge.h"
#include "utils.h"

struct parent_data *loadparents(const char *fname, struct parent_data *parent, int64 Ngroups)
{
    int64 nlines = 0;
    char comment = '#';
    FILE *fp = NULL;
    char str_line[MAXLINESIZE];

    fp = my_fopen(fname, "rt");
    /* fprintf(stderr,"in loadparents reading file '%s' ngroups =
     * %"STR_FMT"\n",fname,Ngroups); */
    while (1)
    {
        if (fgets(str_line, MAXLINESIZE, fp) != NULL)
        {
            if (str_line[0] != comment)
            {
                /* 			  fprintf(fd,"%4d  %14"STR_FMT " %14d  %10d
                 * %14"STR_FMT   "  %14"STR_FMT "  %14"STR_FMT " %10d        %14"STR_FMT
                 * " %14.4f  %12"STR_FMT " \n", */
                /* 					  snapshot_number,i,group0[i].ParentLevel,group0[i].ContainerIndex,group0[i].Nsub,
                 * group0[i].N,group0[i].ParentId,group0[i].ParentSnapshot,group0[i].NpartinParent,
                 */
                /* 					  group0[i].Rank,group0[i].Ncommon);
                 */

                sscanf(str_line,
                       "%hd   %" STR_FMT " %hd  %" STR_FMT " %" STR_FMT "  %" STR_FMT "  %" STR_FMT
                       " %hd        %" STR_FMT " %lf   %" STR_FMT " \n",
                       &parent[nlines].snapshot, &parent[nlines].groupid, &parent[nlines].parentlevel,
                       &parent[nlines].containerid, &parent[nlines].nsub, &parent[nlines].npartinhalo,
                       &parent[nlines].parentid, &parent[nlines].parentsnapshot, &parent[nlines].npartinparent,
                       &parent[nlines].rank, &parent[nlines].ncommon);

                /* Updating for the new format for parent_*.txt
                         ParentLevel is replacing isFof; and containerid
                         is being added in.
                */

                parent[nlines].isFof = (parent[nlines].parentlevel == 1) ? 1 : 0;
                nlines++;
            }
        }
        else
        {
            if (nlines != Ngroups)
            {
                fprintf(stderr,
                        "ERROR: Could not read in all the expected number of groups. "
                        "Ngroups %" STR_FMT " , valid data points read = %" STR_FMT "\n",
                        Ngroups, nlines);
                exit(EXIT_FAILURE);
            }
            break;
        }
    }
    fclose(fp);

    return parent;
}
