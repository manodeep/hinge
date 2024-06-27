#ifndef FOF_ONLY

#include "fillprogenitors.h"
#include "hinge.h"
#include "io.h"         //for definition of groups struct
#include "loadgroups.h" //for loadgroups function definition
#include "read_param.h"
#include "sglib.h"
#include "utils.h"

/* given a group in "source", it matches all the groups
in dest and returns the best match group number [in terms
of Ncommon/rank [depending on match_flag]. ]
*/

#define MATCH_WITH_RANK 1
#define MATCH_WITH_NCOMMON 0

#ifdef USE_INT64_FOR_DEST_ARRAYS
#define DESTGROUP_TYPE int64
#else
#define DESTGROUP_TYPE uint32_t
#endif


/* This is to match the ids of all particles in halos that
do not have a progenitor and return group numbers for all those
particles that are found in some previous snapshot.
*/

// int64 get_best_groupnum_wids(const id64 *sourceIds, const int64 Nids, struct group_data *dest, const int64
// destNgroups,
//                              const int flag, double *rank, id64 *DestPartIds, int64 *DestGroupIds,
//                              const int64 *DestGroupLoc, const id64 DestMaxPartId, const id64 DestMinPartId);
int64 get_best_groupnum_wids(const id64 *sourceIds, const int64 Nids, struct group_data *dest, const int64 destNgroups,
                             const int flag, double *rank, const id64 *DestPartIds, const int64 DestNumPart,
                             const DESTGROUP_TYPE *DestGroupIds, const DESTGROUP_TYPE *DestGroupLoc, const id64 DestMaxPartId,
                             const id64 DestMinPartId);

int compare_id64(const void *a, const void *b);

short load_found_progenitors(struct node_data *tree[], int64 *Ngroups, const char *fname, int64 *startgroupnum)
{
    char comment = '#';
    char str_line[MAXLINESIZE];
    int64 nlines = 0;
    short snapshot, prevsnap, min_snapshot = PARAMS.MAX_SNAPSHOT_NUM;
    int64 groupnum = -1, groupid, groupN, prevnum;
    struct node_data *thisnode = NULL, *BaseNode = NULL, *prevnode = NULL, *tmp_node = NULL;

    FILE *fp = my_fopen(fname, "rt");
    while (fgets(str_line, MAXLINESIZE, fp) != NULL)
    {
        if (str_line[0] == comment)
            continue;

        sscanf(str_line, "%hd  %" STR_FMT "  %" STR_FMT "  %" STR_FMT "  %hd  %" STR_FMT "  \n", &snapshot, &groupnum,
               &groupid, &groupN, &prevsnap, &prevnum);

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
            fprintf(stderr, "Ngroups[%hd] = %" STR_FMT " Ngroups[%hd] = %" STR_FMT " \n", snapshot, Ngroups[snapshot],
                    prevsnap, Ngroups[prevsnap]);
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

            // BigChild will not be switched -- so prevnode
            // has to be at least the Sibling of BigChild
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
        /* Need to fix the haloids -- since it could be a partial
            found_progenitors.txt so the behaviour of a complete fillprogenitor
            run is mimicked.
        */
        tmp_node = prevnode;
        while (tmp_node != NULL)
        {
            tmp_node->haloid = thisnode->haloid;
            tmp_node = tmp_node->BigChild;
        }

        nlines++;
    }

    fclose(fp);
    if (nlines > 0)
    {
        *startgroupnum = groupnum + 1;
        fprintf(stderr, "In load_found_progenitors:  Read in %" STR_FMT " lines\n", nlines);
    }
    else
    {
        min_snapshot = PARAMS.MIN_SNAPSHOT_NUM;
        *startgroupnum = 0;
        fprintf(stderr,
                "In load_found_progenitors:  Incomplete `%s` found - read in %" STR_FMT
                " lines. Recalculating the orphans ...\n",
                fname, nlines);
        PARAMS.LOAD_FOUND_PROGENITORS = 0;
    }

    return min_snapshot;
}

int compare_id64(const void *a, const void *b)
{
    return (*(const id64 *)a - *(const id64 *)b);
}

int64 get_best_groupnum_wids(const id64 *sourceIds, const int64 Nids, struct group_data *dest, const int64 destNgroups,
                             const int flag, double *rank, const id64 *DestPartIds, const int64 DestNumPart,
                             const DESTGROUP_TYPE *DestGroupIds, const DESTGROUP_TYPE *DestGroupLoc, const id64 DestMaxPartId,
                             const id64 DestMinPartId)
{
    double *DestRanks = NULL;
    int64 *DestNcommon = NULL;
    // int64 index, grp_index;
    int64 max_ranknum = 0;
    double max_rank = 0.0;

    XASSERT(destNgroups > 0, "There must be at least one group in the destination = %" STR_FMT "\n", destNgroups);

    DestRanks = my_calloc(sizeof(*DestRanks), destNgroups);
    DestNcommon = my_calloc(sizeof(*DestNcommon), destNgroups);

    // for (int64 i = 0; i < destNgroups; i++)
    // {
    //     DestNcommon[i] = 0;
    //     DestRanks[i] = 0.0;
    // }

    for (int64 i = 0; i < Nids; i++)
    {
        id64 index = sourceIds[i];
        if (index >= DestMaxPartId || index < 0 || index < DestMinPartId)
            continue;
        // #define compare_id64(a, b) ((a) < (b) ? -1 : (a) > (b))
        id64 *ptr = (id64 *)bsearch(&index, DestPartIds, DestNumPart, sizeof(*DestPartIds), compare_id64);
        // #undef compare_id64
        // if (DestPartIds[index] != -1)
        if (ptr == NULL)
            continue;
        // if(ptr != NULL)
        // {
        /* DestPartIds is indexed by the particle ID (and hence requires a large amount of RAM )
            The other arrays are indexed by the *value* in DestPartids - which is really
            the particle offset where the original particle was found (as in a cumulative count over the
            number of particles located in all previous halos plus the number of particles within the
            originating halo - MS 18th June 2024)
        */
        XASSERT(*ptr == index,
                "Error: The particle id = %" STR_ID_FMT " must be equal to the index = %" STR_ID_FMT
                " in the DestPartIds array\n",
                index, *ptr);
        index = ptr - DestPartIds;
        const int64 grp_index = DestGroupIds[index];
        DestNcommon[grp_index]++;

        if (flag == 1)
        {
            if (PARAMS.MAX_RANK_LOC <= 0 || (PARAMS.MAX_RANK_LOC > 0 && i < PARAMS.MAX_RANK_LOC))
                DestRanks[grp_index] += compute_rank(i); /* matching based on rank */

            if (PARAMS.MAX_RANK_LOC <= 0 || (PARAMS.MAX_RANK_LOC > 0 && DestGroupLoc[index] < PARAMS.MAX_RANK_LOC))
                DestRanks[grp_index] += compute_rank(DestGroupLoc[index]);
        }
        else
        {
            DestRanks[grp_index] += 1.0; /* matching based on Ncommon */
        }

        if (DestRanks[grp_index] > max_rank)
        {
            max_rank = DestRanks[grp_index];
            max_ranknum = grp_index;
        }
        // }
    }

    free(DestRanks);
    free(DestNcommon);

    *rank = max_rank;
    if (max_rank > 0.0)
        return dest[max_ranknum].nodeloc; /* dest might not start at 0 -- use nodeloc to return the
                                             actual group number */

    return -1;
}

void fillprogenitors(struct node_data *tree[], int64 *Ngroups)
{
    // There is a small bug in the code where the memory for [startsnapshot] does
    // not get freed until the end.
    int flag;
    int64 searchnodenum;
    short searchsnapshot, startsnapshot = PARAMS.MAX_SNAPSHOT_NUM;
    short NUM_SNAPSHOTS = PARAMS.MAX_SNAPSHOT_NUM + 1;

    char fname[MAXLEN];
    FILE *fp = NULL;
    struct node_data *BaseNode = NULL, *thisnode = NULL, *NewBaseNode = NULL, *checknode = NULL;
    struct node_data *tmp_node = NULL;
    struct group_data *group0 = NULL, *group1 = NULL;
    struct group_data *allgroups[NUM_SNAPSHOTS];
    int64 startgroup = 0;
    time_t t_sectionstart, t_sectionend;

    id64 *DestPartIds[NUM_SNAPSHOTS];
    DESTGROUP_TYPE *DestGroupIds[NUM_SNAPSHOTS];
    DESTGROUP_TYPE *DestGroupLoc[NUM_SNAPSHOTS];

    id64 DestMaxPartId[NUM_SNAPSHOTS];
    id64 DestMinPartId[NUM_SNAPSHOTS];
    int64 numpart_in_halos[NUM_SNAPSHOTS];

    /* This can potentially be optimised. In case
           of massive memory requirements, reduce the following to unsigned ints

           Yes -1 will wrap to UINT_MAX. So, you will get one less than UINT_MAX
       number of particles.

    */

    my_snprintf(fname, MAXLEN, "%s/found_progenitors.txt", PARAMS.OUTPUT_DIR);

    if (PARAMS.LOAD_PARTIAL_FOUND_PROGENITORS == 1)
    {
        startsnapshot = load_found_progenitors(tree, Ngroups, fname, &startgroup);
        fp = my_fopen(fname, "a");
        flag = 1; /* so that assign_haloid *will* be called at the end of this
                     function */
    }
    else
    {
        fp = my_fopen(fname, "w");
        fprintf(fp, "##############################################################"
                    "##############################################################"
                    "###########################################\n");
        fprintf(fp, "# Snapshot      Groupnum            GroupId              "
                    "GroupN       PrevSnap        PrevNum          PrevN           "
                    " Ncommon               Rank           Max_rank  \n");
        fprintf(fp, "#    i             l                   l                    l "
                    "           i               l                l                "
                    "l                    d                d     \n");
        fprintf(fp, "##############################################################"
                    "##############################################################"
                    "###########################################\n");
    }

    for (short isnapshot = PARAMS.MIN_SNAPSHOT_NUM; isnapshot <= PARAMS.MAX_SNAPSHOT_NUM; isnapshot++)
    {
        allgroups[isnapshot] = NULL;
        DestPartIds[isnapshot] = NULL;
        DestGroupIds[isnapshot] = NULL;
        DestGroupLoc[isnapshot] = NULL;
        DestMaxPartId[isnapshot] = -1;
        DestMinPartId[isnapshot] = NUMPART + 1;
        numpart_in_halos[isnapshot] = 0;
    }

    // short snapshot = 0;
    for (short isnapshot = startsnapshot; isnapshot >= PARAMS.MIN_SNAPSHOT_NUM; isnapshot--)
    {
        fprintf(stderr, "\n\nfillprogenitor: Now working on snapshot # %4d Ngroups = %" STR_FMT "\n\n", isnapshot,
                Ngroups[isnapshot]);
        short snapshot = isnapshot + 1;
        int64_t currRealMem = 0, peakRealMem = 0, currVirtMem = 0, peakVirtMem = 0;
        getMemory(&currRealMem, &peakRealMem, &currVirtMem, &peakVirtMem);
        fprintf(stderr,
                "Memory used: Real = %" PRId64 " (peak = %" PRId64 ") bytes, Virtual = %" PRId64 " (peak = %" PRId64
                ") bytes\n",
                currRealMem, peakRealMem, currVirtMem, peakVirtMem);

        /* 	  if(snapshot < startsnapshot && allgroups[snapshot] !=NULL &&
         * Ngroups[snapshot] > 0) */

        if (snapshot <= startsnapshot && allgroups[snapshot] != NULL && Ngroups[snapshot] > 0)
        {
            fprintf(stderr, "freed group (inside original for loop) # %d\n", snapshot);
            free_group(allgroups[snapshot], Ngroups[snapshot]);
            allgroups[snapshot] = NULL;
        }

        /* 	  if(snapshot < startsnapshot && DestPartIds[snapshot] != NULL) */
        if (snapshot <= startsnapshot && DestPartIds[snapshot] != NULL)
        {
            DestMaxPartId[snapshot] = -1;
            my_free((void **)&(DestPartIds[snapshot]));
            my_free((void **)&(DestGroupIds[snapshot]));
            my_free((void **)&(DestGroupLoc[snapshot]));
        }

        if (DestPartIds[isnapshot] != NULL)
        {
            DestMaxPartId[isnapshot] = -1;
            my_free((void **)&(DestPartIds[isnapshot]));
            my_free((void **)&(DestGroupIds[isnapshot]));
            my_free((void **)&(DestGroupLoc[isnapshot]));
        }

        if (Ngroups[isnapshot] > 0)
        {
            for (short incr = 0; incr <= (PARAMS.MAX_DECR_GROUPS + 2); incr++)
            {
                snapshot = isnapshot - incr;
                if (snapshot >= PARAMS.MIN_SNAPSHOT_NUM && allgroups[snapshot] == NULL && Ngroups[snapshot] > 0)
                {
                    fprintf(stderr,
                            "\nIn fillprogenitor: Loading groups for snapshot  %4d, "
                            "Ngroups = %" STR_FMT "....\n",
                            snapshot, Ngroups[snapshot]);
                    getMemory(&currRealMem, &peakRealMem, &currVirtMem, &peakVirtMem);
                    fprintf(stderr,
                            "Memory used (before allocate_group): Real = %" PRId64 " (peak = %" PRId64
                            ") bytes, Virtual = %" PRId64 " (peak = %" PRId64 ") bytes\n",
                            currRealMem, peakRealMem, currVirtMem, peakVirtMem);

                    group0 = allocate_group(Ngroups[snapshot]);
                    loadgroups(&PARAMS, snapshot, group0);
                    allgroups[snapshot] = group0;

                    fprintf(stderr,
                            "\nIn fillprogenitor: Loading groups for snapshot  %4d, "
                            "Ngroups = %" STR_FMT "....done\n",
                            snapshot, Ngroups[snapshot]);
                    getMemory(&currRealMem, &peakRealMem, &currVirtMem, &peakVirtMem);
                    fprintf(stderr,
                            "Memory used (after loadgroups): Real = %" PRId64 " (peak = %" PRId64
                            ") bytes, Virtual = %" PRId64 " (peak = %" PRId64 ") bytes\n",
                            currRealMem, peakRealMem, currVirtMem, peakVirtMem);

                    // remove_duplicates needs particle positions - that's why
                    // we need to load the particle positions in the first place
                    // Once loadgroups returns, we have no need for particle positions
                    free_group_positions(group0, Ngroups[snapshot]);

                    getMemory(&currRealMem, &peakRealMem, &currVirtMem, &peakVirtMem);
                    fprintf(stderr,
                            "Memory used (after loadgroups + free_group_pos): Real = %" PRId64 " (peak = %" PRId64
                            ") bytes, Virtual = %" PRId64 " (peak = %" PRId64 ") bytes\n",
                            currRealMem, peakRealMem, currVirtMem, peakVirtMem);

                    // Find the max particle id
                    id64 max_part_id = -1, min_part_id = NUMPART + 1;
                    numpart_in_halos[snapshot] = 0;
                    for (int64 i = 0; i < Ngroups[snapshot]; i++)
                    {
                        numpart_in_halos[snapshot] += group0[i].N;
                        XASSERT(group0[i].nodeloc == i,
                                "Error: Expected nodeloc = %" STR_FMT " to be *exactly* equal to i=%" STR_FMT "\n",
                                group0[i].nodeloc, i);
                        for (int64 j = 0; j < group0[i].N; j++)
                        {
                            const id64 this_id = group0[i].id[j];
                            if (this_id < 0)
                                continue;
                            max_part_id = this_id > max_part_id ? this_id : max_part_id;
                            min_part_id = this_id < min_part_id ? this_id : min_part_id;
                        }
                    }
                    max_part_id++;
                    if (DestMaxPartId[snapshot] < max_part_id)
                    {
                        fprintf(stderr,
                                "Replacing DestMaxPartId[%d] = %" STR_ID_FMT " with %" STR_ID_FMT " NUMPART = %" STR_FMT
                                "\n",
                                snapshot, DestMaxPartId[snapshot], max_part_id, NUMPART);
                        DestMaxPartId[snapshot] = max_part_id;
                    }
                    if (DestMinPartId[snapshot] > min_part_id)
                    {
                        fprintf(stderr,
                                "Replacing DestMinPartId[%d] = %" STR_ID_FMT " with %" STR_ID_FMT " NUMPART = %" STR_FMT
                                "\n",
                                snapshot, DestMinPartId[snapshot], min_part_id, NUMPART);
                        DestMinPartId[snapshot] = min_part_id;
                    }
                }

                // fill in the data for the destination groups
                if (snapshot >= PARAMS.MIN_SNAPSHOT_NUM && DestPartIds[snapshot] == NULL && incr >= 2)
                {
                    group0 = allgroups[snapshot];
                    DestPartIds[snapshot] = my_malloc(sizeof(*DestPartIds[snapshot]), numpart_in_halos[snapshot]);
                    fprintf(stderr, "Allocating for %" STR_FMT " particles in snapshot %d\n",
                            numpart_in_halos[snapshot], snapshot);
                    fprintf(stderr, "DestMaxPartId[%d] = %" STR_ID_FMT " DestMinPartId[%d] = %" STR_ID_FMT "\n",
                            snapshot, DestMaxPartId[snapshot], snapshot, DestMinPartId[snapshot]);
                    DestGroupIds[snapshot] = my_malloc(sizeof(*DestGroupIds[snapshot]), numpart_in_halos[snapshot]);
                    DestGroupLoc[snapshot] = my_malloc(sizeof(*DestGroupLoc[snapshot]), numpart_in_halos[snapshot]);

                    int64 offset = 0;
                    for (int64 i = 0; i < Ngroups[snapshot]; i++)
                    {
                        for (int64 j = 0; j < group0[i].N; j++)
                        {
                            const id64 this_id = group0[i].id[j];
                            if (this_id < 0)
                                continue;
                            XASSERT(this_id >= 0 && this_id < DestMaxPartId[snapshot],
                                    "Error: Particle id = %" STR_ID_FMT
                                    "must be less than max. particle id = %" STR_ID_FMT " - "
                                    "strange things must have happened",
                                    this_id, DestMaxPartId[snapshot]);
                            DestPartIds[snapshot][offset] = this_id;
#ifndef USE_INT64_FOR_DEST_ARRAYS
                            XASSERT(i < UINT32_MAX, "Error: Group number = %" STR_FMT " must be less than UINT32_MAX\n", i);
                            XASSERT(j < UINT32_MAX, "Error: Particle index = %" STR_FMT " (within group = %"STR_FMT") must be less than UINT32_MAX\n", j, i);
#endif
                            DestGroupIds[snapshot][offset] = i;
                            DestGroupLoc[snapshot][offset] = j;
                            offset++;
                        }
                    }
                    if (offset != numpart_in_halos[snapshot])
                    {
                        fprintf(stderr,
                                "LOG: offset = %" STR_FMT " is different from numpart_in_halos[%d] = %" STR_FMT "\n",
                                offset, snapshot, numpart_in_halos[snapshot]);
                        fprintf(
                            stderr,
                            "LOG: This is likely due to particles with negative ids. Setting 'numpart' to 'offset'\n");
                    }
                    numpart_in_halos[snapshot] = offset;
#define MULTIPLE_ARRAY_EXCHANGER(type, varname, i, j)                                                                  \
    {                                                                                                                  \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(id64, DestPartIds[snapshot], i, j);                                             \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(DESTGROUP_TYPE, DestGroupIds[snapshot], i, j);                                  \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(DESTGROUP_TYPE, DestGroupLoc[snapshot], i, j);                                  \
    }

                    time_t t0 = time(NULL);
                    SGLIB_ARRAY_QUICK_SORT(id64, DestPartIds[snapshot], numpart_in_halos[snapshot],
                                           SGLIB_NUMERIC_COMPARATOR, MULTIPLE_ARRAY_EXCHANGER);
#undef MULTIPLE_ARRAY_EXCHANGER
                    time_t t1 = time(NULL);
                    print_time(t0, t1, "Quick sort in fillprogenitors");
                }
            }

            group0 = allgroups[isnapshot];
            BaseNode = tree[isnapshot];

            // start from group 0 for all other snapshots
            if (isnapshot != startsnapshot)
                startgroup = 0;

            for (int64 igroup = startgroup; igroup < Ngroups[isnapshot]; igroup++)
            {
                thisnode = &BaseNode[igroup];

                /* Does the halo not have a progenitor (while it's Fof has one) */
                if (thisnode->BigChild == NULL && thisnode->isFof == 0 && thisnode->FofHalo->BigChild != NULL)
                {
                    int64 Nids = group0[thisnode->nodeloc].N;
                    id64 *TrackIds = my_malloc(sizeof(*TrackIds), Nids);
                    for (int64 j = 0; j < Nids; j++)
                    {
                        TrackIds[j] = group0[thisnode->nodeloc].id[j];
                    }

                    double max_rank = 0.0;
                    for (int64 i = 0; i < group0[thisnode->nodeloc].N; i++)
                        max_rank += compute_rank(i);

                    /* the fofmatch code would already have tried to find a good match at
                       the previous snapshot. So, if thisnode does not have a proper
                       progenitor -- then the good progenitor, if it exists, must exist at
                       least two snapshots prior -> should start with skip=2.

                             Never mind. Starting with skip = 2. MS 07/18/2011
                    */

                    for (int skip = 2; skip <= PARAMS.MAX_DECR_GROUPS; skip++)
                    {
                        fprintf(stderr,
                                "Inside fillprogenitors: missing haloid = %" STR_FMT " with groupnum = %" STR_FMT
                                " at snapshot = %d skip = %d \n",
                                thisnode->haloid, thisnode->nodeloc, thisnode->snapshot, skip);
                        searchsnapshot = thisnode->snapshot - skip;
                        if (searchsnapshot < PARAMS.MIN_SNAPSHOT_NUM || Ngroups[searchsnapshot] == 0)
                            continue;

                        group1 = allgroups[searchsnapshot];
                        double rank = 0.0;
                        searchnodenum = get_best_groupnum_wids(
                            TrackIds, Nids, group1, Ngroups[searchsnapshot], MATCH_WITH_RANK, &rank,
                            DestPartIds[searchsnapshot], numpart_in_halos[searchsnapshot], DestGroupIds[searchsnapshot],
                            DestGroupLoc[searchsnapshot], DestMaxPartId[searchsnapshot], DestMinPartId[searchsnapshot]);

                        if (searchnodenum != -1)
                        {
                            NewBaseNode = tree[searchsnapshot];
                            checknode = &NewBaseNode[searchnodenum];
                            assert(searchnodenum < Ngroups[searchsnapshot] && "Possible match must be valid");
                            assert(thisnode->nodeloc < Ngroups[thisnode->snapshot] &&
                                   "Candidate node number must be valid");
                            int64 ncommon = get_ncommon(&group0[thisnode->nodeloc], &group1[searchnodenum]);
                            fprintf(stderr,
                                    "found a possible match:  nodenum = %" STR_FMT
                                    " at snapshot = %d. haloid = %" STR_FMT " with ncommon = %" STR_FMT
                                    " out of Npart = %" STR_FMT "\n",
                                    searchnodenum, searchsnapshot, checknode->haloid, ncommon, group1[searchnodenum].N);
                            /* I am taking out the ParentLevel >= 2. As long as the FOF can be
                               re-assigned here, and it is not the main progenitor -- it's
                               fine. (I think) */
                            if (checknode->Parent == NULL ||
                                (checknode->Parent->BigChild != checknode &&
                                 (double)ncommon / (double)group0[thisnode->nodeloc].N > PARAMS.MIN_FCOMMON_THRESH))
                            {

                                if (checknode->Parent != NULL && checknode->Parent->Nchild > 1)
                                {
                                    tmp_node = checknode->Parent->BigChild;
                                    while (tmp_node->Sibling != checknode)
                                        tmp_node = tmp_node->Sibling;

                                    tmp_node->Sibling = checknode->Sibling;
                                    tmp_node->Parent->Nchild--;
                                }

                                checknode->Parent = thisnode;
                                checknode->ParentID = thisnode->nodeloc;
                                checknode->ParentZ = thisnode->z;
                                checknode->ParentSnapshot = thisnode->snapshot;
                                checknode->Sibling = NULL;
                                thisnode->BigChild = checknode;
                                thisnode->Nchild = 1;

                                tmp_node = checknode;
                                while (tmp_node != NULL)
                                {
                                    tmp_node->haloid = thisnode->haloid;
                                    tmp_node = tmp_node->BigChild;
                                }

                                fprintf(stderr,
                                        "Found a progenitor for node with groupnum %" STR_FMT
                                        " at snapshot %d. The progenitor is at snapshot %d, "
                                        "groupnumber = %" STR_FMT " \n",
                                        thisnode->nodeloc, thisnode->snapshot, checknode->snapshot, checknode->nodeloc);
                                fprintf(fp,
                                        "%6d    %12" STR_FMT "    %14" STR_FMT "     %16" STR_FMT
                                        "    %10d   %12" STR_FMT "   %14" STR_FMT "    %14" STR_FMT
                                        "      %16.4lf  %16.4lf\n",
                                        thisnode->snapshot, thisnode->nodeloc, thisnode->haloid,
                                        group0[thisnode->nodeloc].N, checknode->snapshot, checknode->nodeloc,
                                        group1[checknode->nodeloc].N, ncommon, rank, max_rank);

                                fflush(fp);
                                flag = 1;
                                break;
                            }
                        }
                    }
                    free(TrackIds);
                }
            }
        }
    }

    for (int isnapshot = PARAMS.MIN_SNAPSHOT_NUM; isnapshot <= startsnapshot; isnapshot++)
    {
        if (allgroups[isnapshot] != NULL)
        {
            fprintf(stderr, "freeing group for snapshot # %d\n", isnapshot);
            free_group(allgroups[isnapshot], Ngroups[isnapshot]);
            allgroups[isnapshot] = NULL;
        }

        if (DestPartIds[isnapshot] != NULL)
        {
            fprintf(stderr, "freeing particle ids for snapshot # %d\n", isnapshot);
            DestMaxPartId[isnapshot] = -1;
            my_free((void **)&(DestPartIds[isnapshot]));
            my_free((void **)&(DestGroupIds[isnapshot]));
            my_free((void **)&(DestGroupLoc[isnapshot]));
        }
    }

    fclose(fp);

    if (flag == 1)
    {
        t_sectionstart = time(NULL);
        assign_haloid(tree, Ngroups);
        t_sectionend = time(NULL);
        print_time(t_sectionstart, t_sectionend, "assign_haloid (after filling in progenitors)");
    }
}

#undef MATCH_WITH_RANK
#undef MATCH_WITH_NCOMMON

#endif
