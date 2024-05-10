#include "io_subfind_binary.h"
#include "utils.h"

int64 returnNhalo_subfind_binary(const struct params_data *params, const int snapnum, const int fof_only)
{
    char catalogue_fname[MAXLEN];
    if (fof_only)
    {
        my_snprintf(catalogue_fname, MAXLEN, "%s/%s_%03d.fofcat", params->GROUP_DIR, params->GROUP_BASE, snapnum);
    }
    else
    {
        my_snprintf(catalogue_fname, MAXLEN, "%s/%s_%03d.subcat", params->GROUP_DIR, params->GROUP_BASE, snapnum);
    }

    FILE *fcat = my_fopen(catalogue_fname, "r");
    int64 Ngroups;
    my_fread(&Ngroups, sizeof(int64), 1, fcat);
    fclose(fcat);
    return Ngroups;
}

void loadgroups_subfind_binary(int num, const struct params_data *params, struct group_data *group)
{
    char catalogue_fname[MAXLEN];
    char particles_fname[MAXLEN];
    // char parttypes_fname[MAXLEN];
    char partids_fname[MAXLEN];
    my_snprintf(particles_fname, MAXLEN, "%s/%s_%03d.pos", params->GROUP_DIR, params->GROUP_BASE, num);
    // my_snprintf(parttypes_fname, MAXLEN, "%s/%s_%03d.types", params->GROUP_DIR,
    // params->GROUP_BASE, num);
    my_snprintf(partids_fname, MAXLEN, "%s/%s_%03d.ids", params->GROUP_DIR, params->GROUP_BASE, num);
    my_snprintf(catalogue_fname, MAXLEN, "%s/%s_%03d.fofcat", params->GROUP_DIR, params->GROUP_BASE, num);

    FILE *fcat = my_fopen(catalogue_fname, "r");

    int64 Ngroups;
    my_fread(&Ngroups, sizeof(int64), 1, fcat);
    if (Ngroups == 0)
    {
        fclose(fcat);
        fprintf(stderr, "No groups found in file..freeing the group pointer \n");
        free(group);
        return;
    }

    if (Ngroups < 0)
    {
        fprintf(stderr, "Error: Ngroups = %" STR_FMT " must be >= 0\n", Ngroups);
        exit(EXIT_FAILURE);
    }

    // If we have reached here, then Ngroups > 0. MS 20th Nov, 2023
    FILE *fpos = my_fopen(particles_fname, "r");
    // FILE *ftype = my_fopen(parttypes_fname, "r");
    FILE *fid = my_fopen(partids_fname, "r");

    /* Read in the positions  */
    int64 NumPart;
    my_fread(&NumPart, sizeof(int64), 1, fpos);
    /* 	  fprintf(stderr,"There are " STR_FMT " particles in the groups file
     * `%s'\n",NumPart,particles_fname); */
    int64 xx;
    // my_fread(&xx, sizeof(int64), 1, ftype);
    // assert(xx == NumPart && "Numpart in types file should be equal to NumPart
    // in positions file");
    my_fread(&xx, sizeof(int64), 1, fid);
    assert(xx == NumPart && "Numpart in ids file should be equal to NumPart in positions file");

    float *xyzpos = my_malloc(sizeof(*xyzpos), 3 * NumPart);
    float *xpos = my_malloc(sizeof(*xpos), NumPart);
    float *ypos = my_malloc(sizeof(*ypos), NumPart);
    float *zpos = my_malloc(sizeof(*zpos), NumPart);
    // char *type = my_malloc(sizeof(*type), NumPart);
    id64 *id = my_malloc(sizeof(*id), NumPart);

    /* read positions */
    my_fread(xyzpos, sizeof(*xyzpos), 3 * NumPart, fpos);
    for (int64 i = 0; i < NumPart; i++)
    {
        xpos[i] = xyzpos[3 * i + 0];
        ypos[i] = xyzpos[3 * i + 1];
        zpos[i] = xyzpos[3 * i + 2];
    }
    free(xyzpos);

    /* read types */
    // my_fread(type, sizeof(char), NumPart, ftype);

    /* read ids */
    my_fread(id, sizeof(id64), NumPart, fid);

    fclose(fpos);
    // fclose(ftype);
    fclose(fid);

#ifndef FOF_ONLY
    char subhalos_fname[MAXLEN];
    char subprop_fname[MAXLEN];

    my_snprintf(subhalos_fname, MAXLEN, "%s/%s_%03d.subcat", params->GROUP_DIR, params->GROUP_BASE, num);
    my_snprintf(subprop_fname, MAXLEN, "%s/%s_%03d.subprop", params->GROUP_DIR, params->GROUP_BASE, num);

    FILE *fsubcat = my_fopen(subhalos_fname, "r");
    FILE *fsubprop = my_fopen(subprop_fname, "r");

    int64 *GroupSubs = my_malloc(sizeof(*GroupSubs),
                                 Ngroups); // not recasting to (int64 *) -> should be automatic
    fseek(fcat, sizeof(int64) * 2 * Ngroups,
          SEEK_CUR); // Seek over GroupLen and GroupOffset which are not required
                     // for Subfind + subhalos
    my_fread(GroupSubs, sizeof(int64), Ngroups, fcat);

    /* Read in from subhalo catalogue*/
    int64 Nsub;
    my_fread(&Nsub, sizeof(int64), 1, fsubcat);
    int64 *SubLen = my_malloc(sizeof(*SubLen), Nsub);
    int64 *SubOffset = my_malloc(sizeof(*SubOffset), Nsub);

    my_fread(SubLen, sizeof(int64), Nsub, fsubcat);
    my_fread(SubOffset, sizeof(int64), Nsub, fsubcat);
    fclose(fsubcat);

    /* read in from the subprop file*/
    my_fread(&Nsub, sizeof(int64), 1, fsubprop);
    float *SubMtot = my_malloc(sizeof(*SubMtot), Nsub);
    // float *SubMgas = my_malloc(sizeof(*SubMgas), Nsub);

    float *SubCM[3];
#ifdef GET_GROUPVEL
    float *SubCMV[3];
#endif
    for (int i = 0; i < 3; i++)
    {
        SubCM[i] = my_malloc(sizeof(*SubCM[i]), Nsub);
        my_fread(SubCM[i], sizeof(float), Nsub, fsubprop);

#ifdef GET_GROUPVEL
        SubCMV[i] = my_malloc(sizeof(*(SubCMV[i])), Nsub);
        my_fread(SubCMV[i], sizeof(float), Nsub, fsubprop);
#endif
    }

    my_fread(SubMtot, sizeof(float), Nsub, fsubprop);
    // my_fread(SubMgas, sizeof(float), Nsub, fsubprop);

    fclose(fsubprop);

#else

    // Output comes from FOF-only
    float *GroupCM[3];
#ifdef GET_GROUPVEL
    float *GroupCMV[3];
#endif
    int64 *GroupLen = my_malloc(sizeof(*GroupLen), Ngroups);
    int64 *GroupOffset = my_malloc(sizeof(*GroupOffset), Ngroups);

    my_fread(GroupLen, sizeof(int64), (size_t)Ngroups, fcat);
    my_fread(GroupOffset, sizeof(int64), (size_t)Ngroups, fcat);

    float *GroupMtot = my_malloc(sizeof(*GroupMtot), Ngroups);
    float *GroupMgas = my_malloc(sizeof(*GroupMgas), Ngroups);

    char fofprop_fname[MAXLEN];
    my_snprintf(fofprop_fname, MAXLEN, "%s/%s_%03d.fofprop", params->GROUP_DIR, params->GROUP_BASE, num);
    /* Read in from fof_prop file*/
    FILE *fprop = my_fopen(fofprop_fname, "r");
    int64 yy = 0;
    my_fread(&yy, sizeof(int64), 1, fprop);
    assert(yy == Ngroups && "Ngroups should be equal between fofprop and fofcat");
    for (int i = 0; i < 3; i++)
    {
        GroupCM[i] = my_malloc(sizeof(*GroupCM[i]), Ngroups);
        my_fread(GroupCM[i], sizeof(float), Ngroups, fprop);
#ifdef GET_GROUPVEL
        GroupCMV[i] = my_malloc(sizeof(*GroupCMV[i]), Ngroups);
        my_fread(GroupCMV[i], sizeof(float), Ngroups, fprop); /* Warning: This was written directly as a block of 3. I
                                                                 am reading it back in as 3 sets of 1. */
#endif
    }

    my_fread(GroupMtot, sizeof(float), Ngroups, fprop);
    my_fread(GroupMgas, sizeof(float), Ngroups, fprop);

    fclose(fprop);
#endif // end of #ifndef FOF-ONLY

    fclose(fcat);
    /* So everything has been read and all the files have been closed
         Now, lets insert the subhalos into the proper places in the structure
         If we are reading FOF only halos then pretty much nothing needs to be
       done.

         Keep in mind that the subfind does return the FOF halos, all the
       particles present in the subhalos will add up to be <= the original FOF
       halo, i.e., the position array will probably contain elements that do not
       belong to any subhalo.

    */

#ifndef FOF_ONLY

    /* output came from subfind */

    /* Decide which are the parent halos.
         Notice that the loop below goes to Ngroups and NOT Nsub.
    */

    int64 *tempindex = my_malloc(sizeof(*tempindex), Ngroups);

    for (int64 i = 0; i < Nsub; i++)
    {
        group[i].isFof = 0;
        group[i].Nsub = 0;
    }

    int64 sumgroups = 0;
    for (int64 i = 0; i < Ngroups; i++)
    {
        tempindex[i] = sumgroups;
        sumgroups += GroupSubs[i];
    }

    for (int64 i = 0; i < Ngroups; i++)
    {
        if (tempindex[i] >= Nsub)
        {
            fprintf(stderr, "WARNING: Possible case of FOF group falling below "
                            "resolution limit \n");
            fprintf(stderr,
                    "WARNING: i = %" STR_FMT " tempindex = %" STR_FMT " is greater than Nsub = %" STR_FMT
                    " with Ngroups = %" STR_FMT "  \n",
                    i, tempindex[i], Nsub, Ngroups);
        }
        else
        {
            group[tempindex[i]].isFof = 1;
            group[tempindex[i]].Nsub = GroupSubs[i];
        }
    }
    free(tempindex);

    int64 FOF_Parent = 0;
    for (int64 i = 0; i < Nsub; i++)
    {
        /* only the fof's have a definite
             parentlevel at this point. all
             others (subhalos) need to be assigned
             one at a later stage.
        */

        group[i].ParentLevel = -1;
        if (group[i].isFof == 1)
        {
            group[i].ParentLevel = 1;
            FOF_Parent = i;
        }

        group[i].N = SubLen[i];
        group[i].nodeloc = i;
        group[i].snapshot = num;
        group[i].redshift = REDSHIFT[num];

        group[i].Mtot = SubMtot[i];
        // group[i].Mgas = SubMgas[i];
        group[i].xcen = xpos[SubOffset[i]];
        group[i].ycen = ypos[SubOffset[i]];
        group[i].zcen = zpos[SubOffset[i]];

#ifdef GET_GROUPVEL
        group[i].vxcen = SubCMV[0][i];
        group[i].vycen = SubCMV[1][i];
        group[i].vzcen = SubCMV[2][i];
#endif
        group[i].groupnum = i;
        group[i].N_per_wedge = 0;
        /* initialise the parent finding variables*/
        group[i].ParentId = -1;
        group[i].NParents = 0;
        group[i].Switched = 0;
        group[i].ParentSnapshot = -1;
        group[i].Ncommon = 0;
        group[i].Rank = 0.0;
        group[i].NpartinParent = 0;
        group[i].snapshot = (short)num;
        group[i].FOFHalo = FOF_Parent;
        group[i].ContainerIndex = FOF_Parent;

        group[i].x = my_malloc(sizeof(*group[i].x), SubLen[i]);
        group[i].y = my_malloc(sizeof(*group[i].y), SubLen[i]);
        group[i].z = my_malloc(sizeof(*group[i].z), SubLen[i]);
        group[i].id = my_malloc(sizeof(*group[i].id), SubLen[i]);
        // group[i].type = my_malloc(sizeof(*group[i].type), SubLen[i]);

        group[i].parentgroupforparticle = my_malloc(sizeof(*group[i].parentgroupforparticle), SubLen[i]);
        group[i].parentsnapshotforparticle = my_malloc(sizeof(*group[i].parentsnapshotforparticle), SubLen[i]);

        float rmax = 0.0;
        for (int64 j = 0; j < SubLen[i]; j++)
        {
            group[i].x[j] = xpos[SubOffset[i] + j];
            group[i].y[j] = ypos[SubOffset[i] + j];
            group[i].z[j] = zpos[SubOffset[i] + j];
            group[i].id[j] = id[SubOffset[i] + j];
            // group[i].type[j] = (int)type[SubOffset[i] + j];

            /* 			  fprintf(stderr,"group[%d].id[%d]=%ld\n",i,j,group[i].id[j]);
             */

            float temp_rmax =
                sqrt(SQR_PERIODIC(group[i].x[j] - group[i].xcen) + SQR_PERIODIC(group[i].y[j] - group[i].ycen) +
                     SQR_PERIODIC(group[i].z[j] - group[i].zcen));
            if (temp_rmax > rmax)
                rmax = temp_rmax;

            /* initialise the parent finding variables*/
            group[i].parentgroupforparticle[j] = -1;
            group[i].parentsnapshotforparticle[j] = -1;
        }
        group[i].Rmax = rmax;
    }

#else
    /* output came from fof groups */

    for (int64 i = 0; i < Ngroups; i++)
    {
        group[i].N = GroupLen[i];
        group[i].NpartinParent = 0;
        group[i].nodeloc = i;
        group[i].snapshot = num;
        group[i].redshift = REDSHIFT[num];

        group[i].x = my_malloc(sizeof(*group[i].x), GroupLen[i]);
        group[i].y = my_malloc(sizeof(*group[i].y), GroupLen[i]);
        group[i].z = my_malloc(sizeof(*group[i].z), GroupLen[i]);
        group[i].id = my_malloc(sizeof(*group[i].id), GroupLen[i]);
        // group[i].type = my_malloc(sizeof(*group[i].type), GroupLen[i]);

        group[i].parentgroupforparticle = my_malloc(sizeof(*group[i].parentgroupforparticle), GroupLen[i]);
        group[i].parentsnapshotforparticle = my_malloc(sizeof(*group[i].parentsnapshotforparticle), GroupLen[i]);

        /* Fake value since all the FOF are considered parents. Almost */
        group[i].groupnum = i;
        group[i].NParents = 0;
        group[i].Switched = 0;
        group[i].isFof = 1;
        group[i].Nsub = 0;
        group[i].FOFHalo = i;
        group[i].ParentLevel = 1;
        group[i].ContainerIndex = i;

        group[i].N_per_wedge = 0;
        /* initialise the parent finding variables*/
        group[i].ParentId = -1;
        group[i].ParentSnapshot = -1;
        group[i].Ncommon = 0;
        group[i].Rank = 0.0;
        group[i].snapshot = (short)num;

        /*WARNING: Place-holder for future computation of actual cm in FOF_ONLY
         * case*/
        group[i].xcen = xpos[GroupOffset[i]];
        group[i].ycen = ypos[GroupOffset[i]];
        group[i].zcen = zpos[GroupOffset[i]];

#ifdef GET_GROUPVEL
        group[i].vxcen = GroupCMV[0][i];
        group[i].vycen = GroupCMV[1][i];
        group[i].vzcen = GroupCMV[2][i];
#endif
        // mtot = 0.0;
        // mgas = 0.0;

        /* This rmax should not be taken too seriously for FOF groups. The center
         * needs to be re-calculated */
        float rmax = 0.0;

        for (int64 j = 0; j < GroupLen[i]; j++)
        {
            group[i].x[j] = xpos[GroupOffset[i] + j];
            group[i].y[j] = ypos[GroupOffset[i] + j];
            group[i].z[j] = zpos[GroupOffset[i] + j];
            group[i].id[j] = id[GroupOffset[i] + j];
            // group[i].type[j] = (int)type[GroupOffset[i] + j];

            /* 			  mtot+=params->MASSARR[group[i].type[j]]; */
            /* 			  if (group[i].type[j] == 0) */
            /* 				mgas+=params->MASSARR[0]; */

            float temp_rmax =
                sqrt(SQR_PERIODIC(group[i].x[j] - group[i].xcen) + SQR_PERIODIC(group[i].y[j] - group[i].ycen) +
                     SQR_PERIODIC(group[i].z[j] - group[i].zcen));
            if (temp_rmax > rmax)
                rmax = temp_rmax;

            group[i].parentgroupforparticle[j] = -1;
            group[i].parentsnapshotforparticle[j] = -1;
        }
        /* 		  group[i].Mtot = mtot; */
        /* 		  group[i].Mgas = mgas; */
        group[i].Mtot = GroupMtot[i];
        group[i].Mgas = GroupMgas[i];
        group[i].Rmax = rmax;
    }
#endif

    /*free up memory*/
    free(xpos);
    free(ypos);
    free(zpos);
    // free(type);
    free(id);

#ifndef FOF_ONLY
    free(SubMtot);
    // free(SubMgas);
    free(SubLen);
    free(SubOffset);
    free(GroupSubs);
    for (int i = 0; i < 3; i++)
    {
        free(SubCM[i]);
#ifdef GET_GROUPVEL
        free(SubCMV[i]);
#endif
    }
#else
    free(GroupLen);
    free(GroupOffset);
    free(GroupMtot);
    free(GroupMgas);
    for (int i = 0; i < 3; i++)
    {
        free(GroupCM[i]);
#ifdef GET_GROUPVEL
        free(GroupCMV[i]);
#endif
    }
#endif
}
