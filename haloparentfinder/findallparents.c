#include "findallparents.h"
#include "progressbar.h"
#include "read_param.h"
#include "utils.h"

// private functions
void print_fofassign(int64 thisnum, struct group_data *prevgroup, struct group_data *nextgroup, int64 NextNsub,
                     double *NextAllRanks, int64 *NextAllCommon, const char *outpath);
void print_reassign(int64 thisnum, struct group_data *prevgroup, struct group_data *nextgroup, int64 NextNsub,
                    double *NextAllRanks, int64 *NextAllCommon, const char *outpath);
void print_reassign_header(FILE *fp);
void print_fofassign_header(FILE *fp);

void print_reassign(int64 thisnum, struct group_data *prevgroup, struct group_data *nextgroup, int64 NextNsub,
                    double *NextAllRanks, int64 *NextAllCommon, const char *outpath)
{
    int snapshot = prevgroup->snapshot;
    FILE *fp = NULL;
    char fname[MAXLEN];

    my_snprintf(fname, MAXLEN, "%s/reassigned_%03d.txt", outpath, snapshot);
    fp = my_fopen_carefully(fname, &print_reassign_header);

    for (int64 i = 0; i < NextNsub; i++)
    {
        if (NextAllRanks[i] > 0.0)
            fprintf(fp,
                    "%7" STR_FMT "%15" STR_FMT "%12" STR_FMT "%11d      %10" STR_FMT "%11" STR_FMT "      %10" STR_FMT
                    "%14.4f   %12" STR_FMT " %14" STR_FMT "    %16" STR_FMT "     %12" STR_FMT "%14" STR_FMT " \n",
                    thisnum, prevgroup[thisnum].N, prevgroup[thisnum].Nsub, nextgroup[i].snapshot, i, nextgroup[i].N,
                    nextgroup[i].Nsub, NextAllRanks[i], NextAllCommon[i], nextgroup[i].ParentId, nextgroup[i].FOFHalo,
                    nextgroup[nextgroup[i].FOFHalo].N, nextgroup[nextgroup[i].FOFHalo].ParentId);
    }

    fclose(fp);
}

void print_reassign_header(FILE *fp)
{
    fprintf(fp, "####################################################################"
                "####################################################################"
                "#########################################################\n");
    fprintf(fp, "# PrevGrpNum      PrevN        PrevNsub   NextSnap     NextNum       "
                "NextN       NextNsub      NextRank        NextCommon       NextParentId "
                "     NextFOFHalo     NextFOFNp       NextFOFParent #\n");
    fprintf(fp, "####################################################################"
                "####################################################################"
                "#########################################################\n");
}

void print_fofassign(int64 thisnum, struct group_data *prevgroup, struct group_data *nextgroup, int64 NextNsub,
                     double *NextAllRanks, int64 *NextAllCommon, const char *outpath)
{
    int snapshot = prevgroup->snapshot;
    FILE *fp = NULL;
    char fname[MAXLEN];

    my_snprintf(fname, MAXLEN, "%s/fofassigned_%03d.txt", outpath, snapshot);
    fp = my_fopen_carefully(fname, &print_fofassign_header);
    for (int64 i = 0; i < NextNsub; i++)
    {
        if (nextgroup[i].isFof == 1 && NextAllRanks[i] > 0.0)
            fprintf(fp,
                    "%7" STR_FMT "      %10" STR_FMT "     %7" STR_FMT "    %6d       %10" STR_FMT "   %10" STR_FMT
                    "    %7" STR_FMT "    %12.4f "
                    "  %12" STR_FMT "  %14" STR_FMT "  %16" STR_FMT "  %12" STR_FMT "  %14" STR_FMT " \n",
                    thisnum, prevgroup[thisnum].N, prevgroup[thisnum].Nsub, nextgroup[i].snapshot, i, nextgroup[i].N,
                    nextgroup[i].Nsub, NextAllRanks[i], NextAllCommon[i], nextgroup[i].ParentId, nextgroup[i].FOFHalo,
                    nextgroup[nextgroup[i].FOFHalo].N, nextgroup[nextgroup[i].FOFHalo].ParentId);
    }

    fclose(fp);
}

void print_fofassign_header(FILE *fp)
{

    fprintf(fp, "####################################################################"
                "####################################################################"
                "#########################################################\n");
    fprintf(fp, "# PrevGrpNum      PrevN        PrevNsub   NextSnap     NextNum       "
                "NextN       NextNsub      NextRank        NextCommon       NextParentId "
                "     NextFOFHalo     NextFOFNp       NextFOFParent #\n");
    fprintf(fp, "####################################################################"
                "####################################################################"
                "#########################################################\n");
}

/*

This matches FOF at next to FOF at previous. Assumes swiss-cheese
model for the subhalos [such that particles present in subs have
to specifically included in the FOF]. FOF's that don't get assigned
in this function will get assigned during the subhalo assignment
based on binding energy rank. This function will *always* get executed,
however, subhalo assignments might not [depending on the value of compile
time flag FOF_ONLY].

*/

int64 findfofparents(struct group_data *prevgroup, int64 PrevNsub, struct group_data *nextgroup, int64 NextNsub,
                     const char *outpath)
{
    int64 i, j;
    /* int PRINTSTEP=0,SMALLPRINTSTEP=0; */
    int64 NFofHalofound;
    int64 tmp_id, tmp_grpid;
    int64 max_rankid, oldchildid;

    short int *NextAllPartIds = NULL;
    int64 *NextAllGroupIds = NULL;
    int64 *NextAllRealGroupIds = NULL;
    int64 *NextAllRealGroupLocs = NULL;
    double *NextAllRanks = NULL;
    int64 *NextAllCommon = NULL;
    int64 NextMaxPartId = 0;
    int64 FOF_Parent;
    char buf[MAXLEN];

    my_snprintf(buf, MAXLEN, "rm -f %s/fofassigned_%03d.txt", outpath, prevgroup->snapshot);
    system(buf);

    for (i = 0; i < NextNsub; i++)
    {
        for (j = 0; j < nextgroup[i].N; j++)
        {
            if (nextgroup[i].id[j] > NextMaxPartId)
            {
                NextMaxPartId = nextgroup[i].id[j];
            }
        }
    }

    NextMaxPartId++; /* should be able to index with NextMaxPartId -> n+1 elements
                      */
    NextAllPartIds = my_calloc(sizeof(*NextAllPartIds), NextMaxPartId); /* Note use of calloc instead of malloc */
    NextAllGroupIds = my_malloc(sizeof(*NextAllGroupIds), NextMaxPartId);
    NextAllRealGroupIds = my_malloc(sizeof(*NextAllRealGroupIds), NextMaxPartId);
    NextAllRealGroupLocs = my_malloc(sizeof(*NextAllRealGroupLocs), NextMaxPartId);
    NextAllRanks = my_malloc(sizeof(*NextAllRanks), NextNsub);
    NextAllCommon = my_malloc(sizeof(*NextAllCommon), NextNsub);

    FOF_Parent = 0;
    int64 flag = 0;
    for (i = 0; i < NextNsub; i++)
    {
        NextAllRanks[i] = 0.0;
        NextAllCommon[i] = 0;
        if (nextgroup[i].isFof == 1)
            FOF_Parent = i;

        for (j = 0; j < nextgroup[i].N; j++)
        {
            /*Check if particle ids are repeated !*/
            if (NextAllPartIds[nextgroup[i].id[j]] == 1)
            {
                /* fprintf(stderr,"There are duplicate particle ids - this code will not
                 * work \n"); */
                /* fprintf(stderr,"duplicate id = %"STR_ID" in group = %"STR_FMT"
                 * previously seen in group =%"STR_FMT"
                 * \n",nextgroup[i].id[j],i,NextAllRealGroupIds[nextgroup[i].id[j]]); */
                flag++;
            }
            NextAllPartIds[nextgroup[i].id[j]] = 1;
            NextAllGroupIds[nextgroup[i].id[j]] = FOF_Parent; /* Note that all the subhalos are faked as if they are
                                                                 located in the FOF container*/
            NextAllRealGroupIds[nextgroup[i].id[j]] = i;
            NextAllRealGroupLocs[nextgroup[i].id[j]] = j;
        }
    }
    if (flag > 0)
    {
        fprintf(stderr, "Found %" STR_FMT " duplicate particles. Code might not work properly\n", flag);
    }

    NFofHalofound = 0;
    /* if(PrevNsub > 100) { */
    /*   PRINTSTEP = (int)floor(0.1*PrevNsub); */
    /*   SMALLPRINTSTEP = ceil(0.01*PrevNsub) > 1 ? ceil(0.01*PrevNsub):1; */
    /*   fprintf(stderr,"\n\n"); */
    /* } */

    int interrupted = 0;
    init_my_progressbar(PrevNsub, &interrupted);

    FOF_Parent = 0;
    i = 0;
    while (i < PrevNsub)
    {
        if (prevgroup[i].isFof == 1)
            FOF_Parent = i;

        if (prevgroup[FOF_Parent].ParentId < 0)
        {
            j = i;
            /* 		  fprintf(stderr,"Now starting on halo with i= %d,
             * NextNsub = %d  with FOF parent = %d  prev.fofhalo = %d
             * \n",i,PrevNsub,FOF_Parent,prevgroup[j].FOFHalo); */
            while (j < PrevNsub && prevgroup[j].FOFHalo == FOF_Parent)
            {
                /* if(PrevNsub > 100) { */
                /*   if(j%PRINTSTEP == 0) */
                /*     fprintf(stderr,"%d%%",(int)ceil(j/PRINTSTEP)*10); */
                /*   else */
                /*     if(j%SMALLPRINTSTEP==0) */
                /*       fprintf(stderr,"."); */
                /* } */

                my_progressbar(j, &interrupted);

                for (int64 k = 0; k < prevgroup[j].N; k++)
                {
                    tmp_id = prevgroup[j].id[k];
                    if (tmp_id < NextMaxPartId)
                    {
                        if (NextAllPartIds[tmp_id] == 1)
                        {
                            tmp_grpid = NextAllGroupIds[tmp_id];
                            if (nextgroup[tmp_grpid].isFof != 1)
                            {
                                /*only execute this check since we are matching FOF->FOF halos*/
                                fprintf(stderr,
                                        "\n\n\n Groupnum %" STR_FMT
                                        " at snapshot %d is not a FOF parent..exiting \n\n\n",
                                        tmp_grpid, nextgroup[tmp_grpid].snapshot);
                                exit(EXIT_FAILURE);
                            }

                            NextAllRanks[tmp_grpid] += 1.0;
                            NextAllCommon[tmp_grpid] += 1;

                            prevgroup[j].parentgroupforparticle[k] =
                                NextAllRealGroupIds[tmp_id]; /* stores the real halo number
                                                                and not the Fof halo number.
                                                                Group matching is still done
                                                                based on the Fof halo number.
                                                             */
                            prevgroup[j].parentsnapshotforparticle[k] = nextgroup[tmp_grpid].snapshot;
                            nextgroup[NextAllRealGroupIds[tmp_id]]
                                .parentgroupforparticle[NextAllRealGroupLocs[tmp_id]] = j;
                            nextgroup[NextAllRealGroupIds[tmp_id]]
                                .parentsnapshotforparticle[NextAllRealGroupLocs[tmp_id]] = prevgroup[j].snapshot;
                        }
                    }
                }

                j++;
            }

            max_rankid = find_max_rank(NextAllRanks, NextNsub);
            if (max_rankid != -1)
            {
                /*
                         Only assign if this next FOF halo doesnt already have a parent.
                   This does require the groups to be processed in decreasing particle
                   number but that does not NECESSARILY hold true.

                      */

                if (nextgroup[max_rankid].ParentId < 0)
                {
                    prevgroup[i].ParentId = max_rankid;
                    prevgroup[i].ParentSnapshot = nextgroup[max_rankid].snapshot;
                    prevgroup[i].Ncommon = NextAllCommon[max_rankid];
                    prevgroup[i].Rank = NextAllRanks[max_rankid];
                    prevgroup[i].NpartinParent = nextgroup[max_rankid].N;
                    prevgroup[i].NParents = 1; /* not required but for completeness sake*/

                    nextgroup[max_rankid].ParentId = FOF_Parent;
                    nextgroup[max_rankid].Ncommon = NextAllCommon[max_rankid];
                    nextgroup[max_rankid].ParentSnapshot = prevgroup[i].snapshot;
                    nextgroup[max_rankid].NpartinParent = prevgroup[i].N;
                    nextgroup[max_rankid].Rank = NextAllRanks[max_rankid]; /* Assumes that the rank is something
                                                                              that is reversible e.g., binding
                                                                              energy rank would not work */
                    nextgroup[max_rankid].NParents++; /* does the next halo have multiple progenitors ->
                                                         initialisation in loadgroups.c is important */
                    print_fofassign(FOF_Parent, prevgroup, nextgroup, NextNsub, NextAllRanks, NextAllCommon, outpath);
                    NFofHalofound++;
                }
                else
                {
                    /*
                       some other FOF halo already claims that max_rankid is the parent
                       halo for them. Now I have to compare the goodness of the match and
                       determine the true child.

                    */
                    if (nextgroup[max_rankid].Rank < NextAllRanks[max_rankid])
                    {
                        /* This new match is better. Switch children. */
                        oldchildid = nextgroup[max_rankid].ParentId;

                        prevgroup[i].ParentId = max_rankid;
                        prevgroup[i].ParentSnapshot = nextgroup[max_rankid].snapshot;
                        prevgroup[i].Ncommon = NextAllCommon[max_rankid];
                        prevgroup[i].Rank = NextAllRanks[max_rankid];
                        prevgroup[i].NpartinParent = nextgroup[max_rankid].N;

                        /* reset the earlier found halo */
                        prevgroup[oldchildid].ParentId = -1;
                        prevgroup[oldchildid].ParentSnapshot = -1;
                        prevgroup[oldchildid].Ncommon = 0;
                        prevgroup[oldchildid].Rank = 0.0;
                        prevgroup[oldchildid].NpartinParent = 0;

                        nextgroup[max_rankid].ParentId = FOF_Parent;
                        nextgroup[max_rankid].Ncommon = NextAllCommon[max_rankid];
                        nextgroup[max_rankid].ParentSnapshot = prevgroup[i].snapshot;
                        nextgroup[max_rankid].NpartinParent = prevgroup[i].N;
                        nextgroup[max_rankid].Rank = NextAllRanks[max_rankid]; /* Assumes that the rank is Ncommon */

                        print_fofassign(FOF_Parent, prevgroup, nextgroup, NextNsub, NextAllRanks, NextAllCommon,
                                        outpath);
                        /* 					  NFofHalofound++; */ /* Dont
                                                                                     update
                                                                                     NFofHalofound
                                                                                     since
                                                                                     one
                                                                                     halo is
                                                                                     being
                                                                                     rejected
                                                                                   */
                    }
                    else
                    {
                        /*
                           Here I could potentially assign it to some other FOF group. But I
                           will leave this prevgroup[i] halo to be assigned by subfind.
                        */
                    }
                }
            }
            init_all_ranks(NextAllRanks, NextAllCommon, NextNsub);
            i = j;
        }
        else
        {
            if (prevgroup[i].isFof == 1)
                NFofHalofound++;

            while (prevgroup[i].FOFHalo == FOF_Parent)
                i++;
        }
    }

    finish_myprogressbar(&interrupted);

    free(NextAllPartIds);
    free(NextAllGroupIds);
    free(NextAllRealGroupIds);
    free(NextAllRealGroupLocs);
    free(NextAllRanks);
    free(NextAllCommon);

    return NFofHalofound;
}

int64 findallparents(struct group_data *prevgroup, int64 PrevNsub, struct group_data *nextgroup, int64 NextNsub,
                     const int snapshot, const char *outpath)
{

    /* int PRINTSTEP=0,SMALLPRINTSTEP=0; */

    /* MS: 19th July, 2010.
           MAX_RANK_LOC imposes a max. number of core particles
           for which the rank is evaluated -- set in the parameter file.
           if MAX_RANK_LOC is set to a non-positive number -- all particles
           will contribute to the rank.
    */

    /* Also, symmetrizing the binding energy rank. The core of the
           next group will also get factored in.
    */

    int64 NumNotFound = 0;
    int64 Nhalofound = 0;
    int64 tmp_id, tmp_grpid;
    int64 max_rankid;
    double tmp_max_rank;
    int64 tmp_max_rankid;

    short int *NextAllPartIds = NULL;
    int64 *NextAllGroupIds = NULL;
    int64 *NextAllGroupLocs = NULL;
    double *NextAllRanks = NULL;
    int64 *NextAllCommon = NULL;

    int64 NextMaxPartId = 0;
    for (int64 i = 0; i < NextNsub; i++)
        for (int64 j = 0; j < nextgroup[i].N; j++)
            if (nextgroup[i].id[j] > NextMaxPartId)
                NextMaxPartId = nextgroup[i].id[j];

    NextMaxPartId++;

    NextAllPartIds = my_calloc(sizeof(*NextAllPartIds), NextMaxPartId); /* Note use of calloc instead of malloc */
    NextAllGroupIds = my_malloc(sizeof(*NextAllGroupIds), NextMaxPartId);
    NextAllGroupLocs = my_malloc(sizeof(*NextAllGroupLocs), NextMaxPartId);
    NextAllRanks = my_malloc(sizeof(*NextAllRanks), NextNsub);
    NextAllCommon = my_malloc(sizeof(*NextAllCommon), NextNsub);

    for (int64 i = 0; i < NextNsub; i++)
    {
        NextAllRanks[i] = 0.0;
        NextAllCommon[i] = 0;
        for (int64 j = 0; j < nextgroup[i].N; j++)
        {
            NextAllPartIds[nextgroup[i].id[j]] = 1;
            NextAllGroupIds[nextgroup[i].id[j]] = i;
            NextAllGroupLocs[nextgroup[i].id[j]] = j;
        }
    }

    Nhalofound = 0;

    /* if(PrevNsub > 100) { */
    /*   PRINTSTEP = (int)floor(0.1*PrevNsub); */
    /*   SMALLPRINTSTEP = ceil(0.01*PrevNsub) > 1 ? ceil(0.01*PrevNsub):1; */
    /*   fprintf(stderr,"\n\n"); */
    /* } */
    int interrupted = 0;
    init_my_progressbar(PrevNsub, &interrupted);

    for (int64 i = 0; i < PrevNsub; i++)
    {
        /* if(PrevNsub > 100) { */
        /*   if(i%PRINTSTEP == 0) */
        /* 	fprintf(stderr,"%d%%",(int)ceil(i/PRINTSTEP)*10); */
        /*   else */
        /* 	if(i%SMALLPRINTSTEP==0) */
        /* 	  fprintf(stderr,"."); */
        /* } */

        my_progressbar(i, &interrupted);

        if (prevgroup[i].ParentId < 0)
        {
            /* so we dont have a parent for this prev group yet  */

            for (int64 j = 0; j < prevgroup[i].N; j++)
            {
                tmp_id = prevgroup[i].id[j];
                tmp_grpid = -1;
                if (tmp_id < NextMaxPartId)
                { /* NextMaxPartId is actually 1 greater
                     than the actual max particle id. Hence
                     the < instead of an <= */
                    if (NextAllPartIds[tmp_id] == 1)
                    {
                        tmp_grpid = NextAllGroupIds[tmp_id];
                        NextAllCommon[tmp_grpid]++;

                        if ((PARAMS.MAX_RANK_LOC <= 0) ||
                            (PARAMS.MAX_RANK_LOC > 0 && j < PARAMS.MAX_RANK_LOC)) // 0 based indexing -> index
                                                                                  // MAX_RANK_LOC-1 is the
                                                                                  // MAX_RANK_LOC'th element
                            NextAllRanks[tmp_grpid] += compute_rank(j);

                        if ((PARAMS.MAX_RANK_LOC <= 0) ||
                            (PARAMS.MAX_RANK_LOC > 0 && NextAllGroupLocs[tmp_id] < PARAMS.MAX_RANK_LOC))
                            NextAllRanks[tmp_grpid] += compute_rank(NextAllGroupLocs[tmp_id]);

                        prevgroup[i].parentgroupforparticle[j] = tmp_grpid;
                        prevgroup[i].parentsnapshotforparticle[j] = snapshot;

                        nextgroup[tmp_grpid].parentgroupforparticle[NextAllGroupLocs[tmp_id]] = i;
                        nextgroup[tmp_grpid].parentsnapshotforparticle[NextAllGroupLocs[tmp_id]] =
                            prevgroup[i].snapshot;
                    }
                }
            }

            max_rankid = find_max_rank(NextAllRanks, NextNsub);

            /* If a FOF halo is trying to get a subhalo as a parent, then
               check the FOF container of the subhalo actually has a FOF child
               and then assign it to the subhalo.
            */
            if (max_rankid != -1)
            {
                tmp_max_rank = NextAllRanks[max_rankid];
                /* Going to exploit the fact that nextgroup FOF finding would have
                 * resulted in their ParentId being filled up */

                /* 21st April, 2010: Do I need the following section ? */
                if ((prevgroup[i].isFof == 1) && (prevgroup[i].Nsub > 1) && (nextgroup[max_rankid].isFof != 1) &&
                    (nextgroup[(nextgroup[max_rankid].FOFHalo)].ParentId < 0) &&
                    (snapshot - prevgroup[i].snapshot) == 1)
                {
                    print_reassign(i, prevgroup, nextgroup, NextNsub, NextAllRanks, NextAllCommon, outpath);
                    /* 				  tmp_nfound =
                     * find_nonzero_ranks(NextAllRanks,NextNsub,tmp_max_rankids); */
                    NextAllRanks[max_rankid] = 0.0;
                    tmp_max_rankid = find_max_rank(NextAllRanks, NextNsub);
                    if (tmp_max_rankid == -1)
                    {
                        fprintf(stderr, "\n\n WARNING: IT MAY BE THAT A FOF HALO NEEDS TO "
                                        "BE ASSIGNED TO A SUBHALO AFTER ALL\n\n");
                        fprintf(stderr,
                                "original  group loc = %" STR_FMT "  the to be parent subhalo is %" STR_FMT
                                " at snapshot %d\n ",
                                i, max_rankid, snapshot);
                        NumNotFound++;
                    }
                    else
                    {
                        fprintf(stderr, "Force switching subhalo to point to fof at next step:\n");
                        fprintf(stderr, "prevgroupnum = %" STR_FMT "  with nextgroupnum = %" STR_FMT "   \n", i,
                                max_rankid);
                        fprintf(stderr,
                                "instead prevgroup will get assigned to fof halonum = %" STR_FMT " with rank = %g\n",
                                nextgroup[max_rankid].FOFHalo, NextAllRanks[max_rankid]);
                        max_rankid = nextgroup[max_rankid].FOFHalo;
                        tmp_max_rank = NextAllRanks[max_rankid];
                    }
                }

                prevgroup[i].ParentId = max_rankid;
                prevgroup[i].ParentSnapshot = snapshot;
                prevgroup[i].Ncommon = NextAllCommon[max_rankid];
                prevgroup[i].Rank = tmp_max_rank; /*NextAllRank might have gotten reset. */
                prevgroup[i].NpartinParent = nextgroup[max_rankid].N;
                prevgroup[i].NParents = 1; /* not reqd. but for completeness sake.
                                              (don't make the code time reversible) */

                if (nextgroup[max_rankid].ParentId >= 0 && nextgroup[max_rankid].ParentId < PrevNsub)
                {
                    if (prevgroup[i].N > prevgroup[nextgroup[max_rankid].ParentId].N)
                    {
                        nextgroup[max_rankid].ParentId = i;
                        nextgroup[max_rankid].Ncommon = NextAllCommon[max_rankid];
                        nextgroup[max_rankid].ParentSnapshot = prevgroup[i].snapshot;
                    }
                }
                else
                {
                    nextgroup[max_rankid].ParentId = i;
                    nextgroup[max_rankid].Ncommon = NextAllCommon[max_rankid];
                    nextgroup[max_rankid].ParentSnapshot = prevgroup[i].snapshot;
                }

                nextgroup[max_rankid].NParents++; /* does the next halo have multiple progenitors ->
                                                     initialisation in loadgroups.c is important */
            }
            else
            {
                /* so max_rank id was -1 -> this can be a side-effect of using
                   some `N' most bound particles. Should check Ncommon and see if there
                   is a possibility of assigning a parent.
                */
            }

            if (prevgroup[i].ParentId >= 0)
                Nhalofound++;
        }
        else
        {
            /* some previous iteration  of findallparents actually found a parent */
            Nhalofound++;
        }

        init_all_ranks(NextAllRanks, NextAllCommon, NextNsub);
    }

    finish_myprogressbar(&interrupted);

    fprintf(stderr, "\nTotal number not found %" STR_FMT " \n", NumNotFound);

    /*   fprintf(stderr,"done...hence aborting\n"); */
    /*   exit(EXIT_FAILURE); */

    free(NextAllPartIds);
    free(NextAllGroupIds);
    free(NextAllGroupLocs);
    free(NextAllRanks);
    free(NextAllCommon);

    return Nhalofound;
}
