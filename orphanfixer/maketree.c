#include "maketree.h"
#include "defs.h"
#include "loadparents.h"
#include "proto.h"
#include "read_param.h"
#include <assert.h>

void assign_parent(struct node_data *halo, int64 haloid, struct node_data *parenthalo, int64 parentid)
{
    struct node_data *tmp_loc = NULL;

    if (halo == NULL || parenthalo == NULL)
    {
        fprintf(stderr, "ERROR: Either the child halo or the parent halo pointer "
                        "is NULL. Exiting..\n");
        if (halo != NULL)
            fprintf(stderr, "halo groupnum = %" STR_FMT " snapshot = %hd\n", halo->nodeloc, halo->snapshot);
        if (parenthalo != NULL)
            fprintf(stderr, "parenthalo groupnum = %" STR_FMT " snapshot = %hd\n", parenthalo->nodeloc,
                    parenthalo->snapshot);

        exit(EXIT_FAILURE);
    }

    parenthalo[parentid].Nchild++;
    if (parenthalo[parentid].Nchild == 1)
    {
        parenthalo[parentid].BigChild = &halo[haloid];
        halo[haloid].Sibling = NULL;
    }
    else
    {
        tmp_loc = parenthalo[parentid].BigChild;
        while (tmp_loc->Mtot > halo[haloid].Mtot && tmp_loc->Sibling != NULL)
            tmp_loc = tmp_loc->Sibling;

        /* Is the new halo going to be the BigChild ? */
        if (tmp_loc == parenthalo[parentid].BigChild)
        {
            if (tmp_loc->Mtot < halo[haloid].Mtot)
            {
                halo[haloid].Sibling = tmp_loc;
                parenthalo[parentid].BigChild = &halo[haloid];
            }
            else
            {
                halo[haloid].Sibling = tmp_loc->Sibling;
                tmp_loc->Sibling = &halo[haloid];
            }
        }
        else
        {
            if (tmp_loc->Sibling == NULL)
            {
                /* The new halo is at the end*/
                tmp_loc->Sibling = &halo[haloid];
                halo[haloid].Sibling = NULL;
            }
            else
            {
                /* The new halo is somewhere in the middle*/
                halo[haloid].Sibling = tmp_loc->Sibling;
                tmp_loc->Sibling = &halo[haloid];
            }
        }
    }

    halo[haloid].Parent = &parenthalo[parentid];
}

void maketree(struct parent_data *allparents[], int64 *Ngroups, struct node_data *tree[])
{
    struct parent_data *p = NULL;
    struct node_data *ParentNode = NULL, *BaseNode = NULL;
    int parentsnapshot;
    int64 parentid;
    int PRINTSTEP, SMALLPRINTSTEP;
    fprintf(stderr, "\n\nAssigning parents..");

    PRINTSTEP = (int)floor(0.1 * (PARAMS.MAX_SNAPSHOT_NUM - PARAMS.MIN_SNAPSHOT_NUM));
    SMALLPRINTSTEP = ceil(0.01 * (PARAMS.MAX_SNAPSHOT_NUM - PARAMS.MIN_SNAPSHOT_NUM)) > 1
                         ? ceil(0.01 * (PARAMS.MAX_SNAPSHOT_NUM - PARAMS.MIN_SNAPSHOT_NUM))
                         : 1;

    for (int isnapshot = PARAMS.MIN_SNAPSHOT_NUM; isnapshot < PARAMS.MAX_SNAPSHOT_NUM; isnapshot++)
    {

        if (SMALLPRINTSTEP > 0 && PRINTSTEP > 0)
        {
            if (isnapshot == PARAMS.MIN_SNAPSHOT_NUM)
                fprintf(stderr, "\n");

            if ((isnapshot - PARAMS.MIN_SNAPSHOT_NUM) % PRINTSTEP == 0)
                fprintf(stderr, "%d%%", (int)ceil((isnapshot - PARAMS.MIN_SNAPSHOT_NUM) / PRINTSTEP) * 10);
            else if ((isnapshot - PARAMS.MIN_SNAPSHOT_NUM) % SMALLPRINTSTEP == 0)
                fprintf(stderr, ".");
        }

        p = allparents[isnapshot];
        BaseNode = tree[isnapshot];
        for (int64 igroup = 0; igroup < Ngroups[isnapshot]; igroup++)
        {
            parentsnapshot = p[igroup].parentsnapshot;
            parentid = p[igroup].parentid;
            if (parentsnapshot >= PARAMS.MIN_SNAPSHOT_NUM && parentsnapshot <= PARAMS.MAX_SNAPSHOT_NUM &&
                parentid < Ngroups[parentsnapshot])
            {
                ParentNode = tree[parentsnapshot];
                assign_parent(BaseNode, igroup, ParentNode, parentid);
            }
        }
    }
    fprintf(stderr, "..done\n");
}

void assign_haloid(struct node_data *tree[], int64 *Ngroups)
{
    int64 haloid = 0;
    struct node_data *BaseNode = NULL;
    struct node_data *thisnode = NULL;
    struct node_data *tmp_node = NULL;

    fprintf(stderr, "\n\n Assigning halo ids to all halos \n\n");
    haloid = 0;

    /* reset all haloid's (in case the tree has been modified) */
    for (int isnapshot = PARAMS.MIN_SNAPSHOT_NUM; isnapshot <= PARAMS.MAX_SNAPSHOT_NUM; isnapshot++)
    {
        BaseNode = tree[isnapshot];
        if (Ngroups[isnapshot] > 0)
        {
            for (int64 igroup = 0; igroup < Ngroups[isnapshot]; igroup++)
            {
                thisnode = &BaseNode[igroup];
                thisnode->haloid = -1;

                /* Make sure that thisnode->Nchild is correct */
                if (thisnode->BigChild != NULL)
                {
                    thisnode->Nchild = 1;
                    tmp_node = thisnode->BigChild->Sibling;
                    while (tmp_node != NULL)
                    {
                        thisnode->Nchild++;
                        tmp_node = tmp_node->Sibling;
                    }
                }
            }
        }
    }

    for (int isnapshot = PARAMS.MAX_SNAPSHOT_NUM; isnapshot >= PARAMS.MIN_SNAPSHOT_NUM; isnapshot--)
    {
        BaseNode = tree[isnapshot];
        if (Ngroups[isnapshot] > 0)
        {
            for (int64 igroup = 0; igroup < Ngroups[isnapshot]; igroup++)
            {
                thisnode = &BaseNode[igroup];
                if (thisnode->haloid < 0)
                {
                    thisnode->haloid = haloid;
                    while (thisnode->BigChild != NULL)
                    {
                        thisnode = thisnode->BigChild;
                        if (thisnode->haloid >= 0)
                        {
                            fprintf(stderr, "\n This should not have happened.. found a "
                                            "haloid while assigning haloids -- exiting \n");
                            fprintf(stderr,
                                    "snapshot = %d this haloid = %" STR_FMT "  this nodeloc = %" STR_FMT
                                    " this parent haloid  = %" STR_FMT "  this parent nodeloc = %" STR_FMT
                                    " snapshot = %hd\n",
                                    thisnode->snapshot, thisnode->haloid, thisnode->nodeloc, thisnode->Parent->haloid,
                                    thisnode->Parent->nodeloc, thisnode->Parent->snapshot);
                            fprintf(stderr,
                                    "thisnode->parent->nchild = %" STR_FMT
                                    "   thisnode->parent->bigchild->haloid = %" STR_FMT " at snapshot = %d \n",
                                    thisnode->Parent->Nchild, thisnode->Parent->BigChild->haloid,
                                    thisnode->Parent->BigChild->snapshot);

                            exit(EXIT_FAILURE);
                        }
                        else
                        {
                            thisnode->haloid = haloid;
                        }
                    }
                    haloid++;
                }
            }
        }
    }
    MaxHaloId = haloid - 1;
}
