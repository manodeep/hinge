#ifndef FOF_ONLY

#include "loadfillprogenitors.h"
#include "hinge.h"
#include "read_param.h"
#include "utils.h"

void load_found_progenitors(struct node_data *tree[], int64 *Ngroups, const char *fname)
{
    FILE *fp = NULL;
    char comment = '#';
    char str_line[MAXLINESIZE];
    int64 nlines = 0;
    short snapshot, prevsnap, min_snapshot = PARAMS.MAX_SNAPSHOT_NUM;
    int64 groupnum = -1, groupid, groupN, prevnum;
    struct node_data *thisnode = NULL, *BaseNode = NULL, *prevnode = NULL, *tmp_node = NULL;

    fp = my_fopen(fname, "rt");

    while (1)
    {
        if (fgets(str_line, MAXLINESIZE, fp) != NULL)
        {
            if (str_line[0] != comment)
            {
                sscanf(str_line, "%hd  %" STR_FMT "  %" STR_FMT "  %" STR_FMT "  %hd  %" STR_FMT "  \n", &snapshot,
                       &groupnum, &groupid, &groupN, &prevsnap, &prevnum);

                if (snapshot < PARAMS.MIN_SNAPSHOT_NUM || snapshot > PARAMS.MAX_SNAPSHOT_NUM)
                {
                    fprintf(stderr,
                            "ERROR: I have snapshot = %hd  but max snapshot = %d and min "
                            "snapshot = %d\n",
                            snapshot, PARAMS.MAX_SNAPSHOT_NUM, PARAMS.MIN_SNAPSHOT_NUM);
                    fprintf(stderr, "exiting..\n");
                    exit(EXIT_FAILURE);
                }

                if (prevsnap < PARAMS.MIN_SNAPSHOT_NUM || prevsnap > PARAMS.MAX_SNAPSHOT_NUM)
                {
                    fprintf(stderr,
                            "ERROR: I have prevsnap = %hd  but max snapshot = %d and min "
                            "snapshot = %d\n",
                            prevsnap, PARAMS.MAX_SNAPSHOT_NUM, PARAMS.MIN_SNAPSHOT_NUM);
                    fprintf(stderr, "exiting..\n");
                    exit(EXIT_FAILURE);
                }

                if (groupnum > Ngroups[snapshot] || prevnum > Ngroups[prevsnap])
                {
                    fprintf(stderr, "Inside load_found_progenitors: This should not have "
                                    "happened\n ");
                    fprintf(stderr, "snapshot = %hd   groupnum = %" STR_FMT " prevsnap = %hd prevnum = %" STR_FMT " \n",
                            snapshot, groupnum, prevsnap, prevnum);
                    fprintf(stderr, "Ngroups[%hd] = %" STR_FMT " Ngroups[%hd] = %" STR_FMT " \n", snapshot,
                            Ngroups[snapshot], prevsnap, Ngroups[prevsnap]);
                    fprintf(stderr, "exiting \n\n");
                    exit(EXIT_FAILURE);
                }

                /* necessary for partial found_progenitors.txt */
                if (snapshot < min_snapshot)
                    min_snapshot = snapshot;

                /* now switch prevnode[prevnum] to point to  currnode[groupnum] */
                BaseNode = tree[snapshot];
                thisnode = &BaseNode[groupnum];

                BaseNode = tree[prevsnap];
                prevnode = &BaseNode[prevnum];

                if (prevnode->Parent != NULL)
                {
                    tmp_node = prevnode->Parent->BigChild;
                    while (tmp_node->Sibling != prevnode)
                        tmp_node = tmp_node->Sibling;

                    tmp_node->Sibling = prevnode->Sibling;
                    prevnode->Parent->Nchild--;
                }

                prevnode->Sibling = NULL;
                prevnode->Parent = thisnode;
                prevnode->ParentSnapshot = thisnode->snapshot;
                prevnode->ParentID = thisnode->nodeloc;
                prevnode->ParentZ = thisnode->z;
                thisnode->BigChild = prevnode;
                thisnode->Nchild = 1;
                /* 			  thisnode->FormationRedshift =
                 * prevnode->FormationRedshift; */
                tmp_node = thisnode;
                while (tmp_node != NULL && tmp_node->haloid == thisnode->haloid)
                {
                    tmp_node->FormationRedshift = prevnode->FormationRedshift;
                    tmp_node = tmp_node->Parent;
                }

                /* Need to fix the haloids -- since it could be a partial
                   found_progenitors.txt so the behaviour of a complete fillprogenitor
                   run is mimicked.
                */

                tmp_node = prevnode;
                while (tmp_node != NULL)
                {
                    tmp_node->haloid = thisnode->haloid;
                    tmp_node->DestructionRedshift = thisnode->DestructionRedshift;
                    tmp_node = tmp_node->BigChild;
                }

                nlines++;
            }
        }
        else
        {
            break;
        }
    }

    fclose(fp);
    fprintf(stderr, "In load_found_progenitors:  Read in %" STR_FMT " lines\n", nlines);
}

#endif /* Don't run this if only FOF halos are being loaded in */
