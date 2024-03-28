#include "switchfof.h"
#include "read_param.h" //for global definition of PARAMS
#include "utils.h"      //for many helper routines

// private function
void switch_parents(struct group_data *prev, struct group_data *next, const int64 ncommon, const double rank,
                    const int snapshot, FILE *fp);

void switch_parents(struct group_data *prev, struct group_data *next, const int64 ncommon, const double rank,
                    const int snapshot, FILE *fp)
{
    fprintf(fp,
            "\n\n switching prevgroupnum = %" STR_FMT " to point to nextgroupnum = %" STR_FMT
            " with ncommon  = %" STR_FMT " \n\n",
            prev->groupnum, next->groupnum, ncommon);

    prev->ParentId = next->groupnum;
    prev->ParentSnapshot = snapshot;
    prev->Ncommon = ncommon;
    prev->Rank = rank;
    prev->NpartinParent = next->N;
    prev->Switched = 1;

    next->ParentId = prev->groupnum;
    next->ParentSnapshot = prev->snapshot;
    next->Ncommon = ncommon;
    next->Rank = rank;
    next->NpartinParent = prev->N;
    next->Switched = 1;
}

void check_fof_matches(struct group_data *prevgroup, int64 PrevNsub, struct group_data *nextgroup, int64 NextNsub,
                       const int snapshot, const char *outpath)
{
    // snapshot corresponds to nextgroup
    int64 parentid, prevparentid, nswitched = 0;
    int64 *FeasibleNcommon = NULL;
    int64 flag = 0;
    FILE *fp = NULL;
    char fname[MAXLEN];

    my_snprintf(fname, MAXLEN, "%s/switchfofhalos_%03d.txt", outpath, snapshot);

    FeasibleNcommon = my_calloc(sizeof(*FeasibleNcommon), PrevNsub);
    for (int64 igroup = 0; igroup < NextNsub; igroup++)
    {
        parentid = nextgroup[igroup].ParentId;
        if (nextgroup[igroup].N > PARAMS.MIN_NUMPART_IN_SWITCHFOF_HALO && parentid >= 0 && parentid < PrevNsub)
        {
            /* so the halo is big enough and does have at least one valid progenitor
                           Note that even though many prevgroups can point to the same
               nextgroup, nextgroup will only contain the biggest of all such
               prevgroups as it's parentid.
            */
            if ((double)prevgroup[parentid].N / (double)nextgroup[igroup].N < PARAMS.MIN_FCOMMON_SWITCHFOF_THRESH)
            {
                fprintf(stderr,
                        "Conflict: for nextgroup = %" STR_FMT " nextgroup.nparents = %" STR_FMT " parentid = %" STR_FMT
                        "  prevncommon = %" STR_FMT " next.N = %" STR_FMT "  \n",
                        igroup, nextgroup[igroup].NParents, parentid, prevgroup[parentid].Ncommon, nextgroup[igroup].N);

                for (int64 i = 0; i < nextgroup[igroup].N; i++)
                {
                    if (nextgroup[igroup].parentgroupforparticle[i] >= 0 &&
                        nextgroup[igroup].parentsnapshotforparticle[i] == prevgroup->snapshot)
                    {
                        FeasibleNcommon[nextgroup[igroup].parentgroupforparticle[i]]++;
                    }
                }

                flag = 0;
                for (int64 i = 0; i < PrevNsub; i += prevgroup[i].Nsub)
                {
                    if (FeasibleNcommon[i] > 0)
                    {
                        fprintf(stderr, "prevgroup = %" STR_FMT " with ncommon = %" STR_FMT " \n", i,
                                FeasibleNcommon[i]);
                    }

                    if ((double)FeasibleNcommon[i] / (double)nextgroup[igroup].N > PARAMS.MIN_FCOMMON_SWITCHFOF_THRESH)
                    {
                        /* possible match: let's see if the prevgroup[i] can be switched to
                         * point to nextgroup[igroup] */
                        prevparentid = prevgroup[i].ParentId; /* nextgroup[prevparentid] was
                                                                 the original parent */
                        fprintf(stderr,
                                "prevgroup = %" STR_FMT " with ncommon = %" STR_FMT
                                " is a potential match..checking all other subhalos \n",
                                i, FeasibleNcommon[i]);
                        if (prevparentid >= 0)
                        {
                            for (int64 j = 0; j < PrevNsub; j++)
                            {
                                if ((double)prevgroup[j].Ncommon / (double)nextgroup[prevparentid].N >
                                        PARAMS.MIN_FCOMMON_SWITCHFOF_THRESH &&
                                    prevgroup[j].ParentId == prevparentid &&
                                    prevgroup[j].ParentSnapshot == nextgroup[prevparentid].snapshot && j != i)
                                {
                                    fprintf(stderr,
                                            "prevgroupnum = %" STR_FMT " is a valid progenitor for nextgroup=%" STR_FMT
                                            " with ncommon1=%" STR_FMT " so that prevfof=%" STR_FMT
                                            " can point to next=%" STR_FMT " \n",
                                            j, prevparentid, prevgroup[j].Ncommon, i, igroup);

                                    nextgroup[prevparentid].ParentId = j;
                                    nextgroup[prevparentid].Ncommon = prevgroup[j].Ncommon;
                                    nextgroup[prevparentid].Rank = prevgroup[j].Rank;
                                    nextgroup[prevparentid].ParentSnapshot = prevgroup->snapshot;
                                    nextgroup[prevparentid].NpartinParent = prevgroup[j].N;
                                    nextgroup[prevparentid].NParents--;

                                    if (fp == NULL)
                                    {
                                        fp = my_fopen(fname, "w");
                                    }

                                    switch_parents(&prevgroup[i], &nextgroup[igroup], FeasibleNcommon[i],
                                                   (double)FeasibleNcommon[i], snapshot, fp);
                                    nswitched++;
                                    flag = 1;
                                    break;
                                }
                            }
                        }
                        else
                        {
                            if (fp == NULL)
                            {
                                fp = my_fopen(fname, "w");
                            }

                            switch_parents(&prevgroup[i], &nextgroup[igroup], FeasibleNcommon[i],
                                           (double)FeasibleNcommon[i], snapshot, fp);
                            nswitched++;
                            flag = 1;
                            break;
                        }

                        if (flag != 0)
                        {
                            break;
                        }
                    }
                }

                reset_ncommon(FeasibleNcommon, PrevNsub);
            }
        }
    }

    free(FeasibleNcommon);
    if (fp != NULL)
    {
        fclose(fp);
    }

    if (nswitched > 0)
    {
        fprintf(stderr, "Inside check_fof_matches:   nswitched halos = %" STR_FMT " \n", nswitched);
    }
}
