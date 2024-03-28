
#ifndef FOF_ONLY

#include "findprogenitor.h"
#include "defs.h"
#include "read_param.h"
#include "utils.h"

void find_progenitor(struct group_data *nextgroup, int64 NextNsub,
                     struct group_data *prevgroup, int64 PrevNsub,
                     const char *outpath) {
  char outfname[MAXLEN];
  FILE *fp = NULL;
  int64 best_match;
  short found = 0;
  int64 old_parentid;
  int64 Nfixed = 0;

  if (nextgroup->snapshot == PARAMS.MIN_SNAPSHOT_NUM)
    return;

  int64 *PrevNcommon = NULL;
  PrevNcommon = my_calloc(sizeof(*PrevNcommon), PrevNsub);

  for (int64 igroup = 0; igroup < NextNsub; igroup++) {
    if (nextgroup[igroup].ParentId >= 0 || nextgroup[igroup].isFof == 1)
      continue;

    /* nextgroup[igroup] is a subhalo that does not have a progenitor
           It's theoretically impossible for a subhalo to form out of nowhere
           inside a Fof (unless there is splitting of another subhalo).
           So, subhalos *should* really have a progenitor.

           It would be really nice if I could update the .NParents variable
       inside findallparents.c even for nextgroup[igroup]'s that don't get
       assigned.
    */

    if (nextgroup[igroup].ParentId < 0 && nextgroup[igroup].isFof == 0 &&
        nextgroup[igroup].N > PARAMS.MIN_NUMPART_IN_FINDPROGENITOR_HALO) {
      best_match = -1;
      for (int64 j = 0; j < nextgroup[igroup].N; j++)
        if (nextgroup[igroup].parentsnapshotforparticle[j] ==
            prevgroup->snapshot) /* all of prevgroup is at same snapshot, this
                                    compares prevgroup[0] */
          PrevNcommon[nextgroup[igroup].parentgroupforparticle[j]]++;

      found = 0;
      best_match = find_max_ncommon(PrevNcommon, PrevNsub);
      while (found == 0 && best_match != -1) {
        old_parentid = prevgroup[best_match].ParentId;
        if (old_parentid < 0) {
          fprintf(stderr, "\nThis is unusual: Nextgroup is progenitorless but "
                          "the best match is descendant-less \n");
          fprintf(stderr, "A side-effect of using some fixed number of most "
                          "bound particles ..\n");
          fprintf(stderr,
                  "igroup = %" STR_FMT " parentid = %" STR_FMT
                  " best_match = %" STR_FMT " old_parentid = %" STR_FMT " \n",
                  igroup, nextgroup[igroup].ParentId, best_match, old_parentid);

          prevgroup[best_match].ParentId = igroup;
          prevgroup[best_match].ParentSnapshot = nextgroup->snapshot;
          prevgroup[best_match].Ncommon = PrevNcommon[best_match];
          prevgroup[best_match].Rank = (double)PrevNcommon[best_match];
          prevgroup[best_match].NpartinParent = nextgroup[igroup].N;
          prevgroup[best_match].NParents = 1;

          nextgroup[igroup].ParentId = best_match;
          nextgroup[igroup].ParentSnapshot = prevgroup->snapshot;
          nextgroup[igroup].Ncommon = PrevNcommon[best_match];
          nextgroup[igroup].NpartinParent = prevgroup[best_match].N;
          nextgroup[igroup].Rank = (double)PrevNcommon[best_match];
          nextgroup[igroup].NParents = 1;

          break;
        }

        if (nextgroup[old_parentid].NParents > 1 &&
            nextgroup[old_parentid].ParentId != best_match &&
            (double)PrevNcommon[best_match] / (double)nextgroup[igroup].N >
                PARAMS.MIN_FCOMMON_FINDPROGENITOR_THRESH &&
            prevgroup[best_match].Switched == 0) {
          /* Relies on the fact that nextgroup will point to the best possible
             progenitor. If the code is here, nextgroup[old_parentid] claims
             some other prevgroup (not group number best_match) as it's real
             progenitor. In that case, switch prevgroup[best_match] to point to
             nextgroup[igroup]..
          */

          fprintf(stderr,
                  "\n Inside findprogenitor: found a progenitor for groupnum "
                  "%" STR_FMT " : %" STR_FMT " with ncommon = %" STR_FMT "  \n",
                  igroup, best_match, PrevNcommon[best_match]);
          fprintf(stderr, "switching progenitor groupnum = %" STR_FMT " ",
                  best_match);
          fprintf(stderr,
                  " with old_parentid = %" STR_FMT
                  " at snapshot = %d. Parentid for old_parentid = %" STR_FMT
                  "  \n",
                  old_parentid, prevgroup[best_match].ParentSnapshot,
                  nextgroup[old_parentid].ParentId);

          prevgroup[best_match].ParentId = igroup;
          prevgroup[best_match].ParentSnapshot = nextgroup->snapshot;
          prevgroup[best_match].Ncommon = PrevNcommon[best_match];
          prevgroup[best_match].Rank = (double)PrevNcommon[best_match];
          prevgroup[best_match].NpartinParent = nextgroup[igroup].N;
          prevgroup[best_match].NParents = 1;

          nextgroup[igroup].ParentId = best_match;
          nextgroup[igroup].ParentSnapshot = prevgroup->snapshot;
          nextgroup[igroup].Ncommon = PrevNcommon[best_match];
          nextgroup[igroup].NpartinParent = prevgroup[best_match].N;
          nextgroup[igroup].Rank = (double)PrevNcommon[best_match];
          nextgroup[igroup].NParents = 1;

          if (Nfixed == 0) {
            my_snprintf(outfname, MAXLEN, "%s/findprogenitor_%03d.txt", outpath,
                        prevgroup->snapshot);
            fp = my_fopen(outfname, "w");
            fprintf(fp, "######################################################"
                        "######################################\n");
            fprintf(fp, "#     PrevGrpNum      PrevN        NextSnap           "
                        " NextNum       NextN      NextCommon  \n");
            fprintf(fp, "#         l             l             i               "
                        "    l            l             l      \n");
            fprintf(fp, "######################################################"
                        "######################################\n");
          }

          fprintf(fp,
                  " %12" STR_FMT " %12" STR_FMT "  %14d    %14" STR_FMT
                  "  %12" STR_FMT "  %12" STR_FMT "   \n",
                  best_match, prevgroup[best_match].N, nextgroup->snapshot,
                  igroup, nextgroup[igroup].N, PrevNcommon[best_match]);

          /* break out of the while loop */
          found = 1;
          Nfixed++;
        } else {
          PrevNcommon[best_match] = 0;
          best_match = find_max_ncommon(PrevNcommon, PrevNsub);
        }
      }

      reset_ncommon(PrevNcommon, PrevNsub);
    }
  }

  fprintf(stderr, "\n\n Finished findprogenitor: fixed %" STR_FMT " halos \n\n",
          Nfixed);
  free(PrevNcommon);
  if (Nfixed > 0)
    fclose(fp);
}

#endif
