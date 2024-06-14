/*

Log started: July 2009: MS.

07/01/2009:   started writing the codes. Complete with a Makefile.

07/19/2009:   loadgroups now written and compiles. Needs to be cross-matched
              with known-good IDL routines to make sure that the data is read
                          in correctly.

08/24/2009:   No more parent checking. Just read in the parents_xxx.txt and
              assign them. Added 3 categories for mergers:
July, 2011: MS

Using the skeleton to make the fixing of the orphans a much easier job

April, 2012:  Complete with a parameter file for the upcoming public release.

*/

#include <stdio.h>
#include <stdlib.h>

#include "hinge.h"
#include "fillprogenitors.h"
#include "io.h"
#include "loadgroups.h"
#include "loadsnapshot.h"
#include "maketree.h"
#include "missinghalos.h"
#include "proto.h"
#include "read_param.h"
#include "set_cosmology.h"
#include "utils.h"
#include "macros.h"
#include "utils_read_params.h"

struct params_data PARAMS; //global variable
float *REDSHIFT;
int64 NUMPART;//Can not figure out how NUMPART is/was calculated but it is used in the code. Needs further investigation.

int main(int argc, char **argv)
{
    char outfname[MAXLEN];

    int64 Ngroups0 = 0;
    int NUM_SNAPSHOTS;

#if (defined(GET_GROUPVEL) == 0 && defined(GET_MEANVEL) == 0)
    struct particle_data *P = NULL;
#endif
    struct group_data *group0 = NULL;
    struct node_data *node = NULL;
    struct parent_data *parent = NULL;
    struct cosmology_data COSMO;

    /* FILE *fp=NULL; */
    /* char str_line[MAXLEN]; */
    /* int line=0; */
    time_t t_codestart, t_codeend, t_sectionstart, t_sectionend, t_bigsectionstart;

    t_codestart = time(NULL);

    /* Check compilation options */
#if ((defined(WMAP1) + defined(WMAP3) + defined(WMAP5)) > 1)
#error Please select only one of the cosmologies from various WMAP options in the Makefile
#endif

#ifdef BIGSIM
    if (sizeof(size_t) != 8)
    {
        fprintf(stderr, "Error: Code needs to be compiled in 64 bit mode \n");
        fprintf(stderr, "Please add -m64 to the options in the Makefile ..(or find "
                        "a 64 bit compiler)\n");
        fprintf(stderr, "Exiting..\n");
        exit(EXIT_FAILURE);
    }
#endif

    // Check command line for the parameter file name.
    if (argc != 2)
    {
        fprintf(stderr, "Usage: %s <parameterfile>\n", argv[0]);
        fprintf(stderr, "\n\n\tCode was compiled with \n\n");
        print_makefile_options();
        exit(EXIT_FAILURE);
    }

    my_snprintf(outfname, MAXLEN, "%s", argv[1]); // parameter file

    // read in the parameter file
    fprintf(stderr, "reading parameter file `%s'...", outfname);
    read_params(outfname, &PARAMS, orphanfixer_fill_params);
    fprintf(stderr, "..done\n");

    // the declarations that depend on MAX_SNAPSHOT_NUM
    // If MAX_SNAPSHOT_NUM is too large, consider mallocing these
    NUM_SNAPSHOTS = PARAMS.MAX_SNAPSHOT_NUM + 1;
    REDSHIFT = my_malloc(sizeof(*REDSHIFT), NUM_SNAPSHOTS);
    int64 Ngroups[NUM_SNAPSHOTS];                            // there are groups corresponding to
                                                             // MAX_SNAPSHOT_NUM
    struct parent_data *allparents[PARAMS.MAX_SNAPSHOT_NUM]; // parents_??? files only go up to
                                                             // MAX_SNAPSHOT_NUM-1
    struct node_data *tree[NUM_SNAPSHOTS];                   // there are groups in MAX_SNAPSHOT_NUM

    // output the parameter file
    my_snprintf(outfname, MAXLEN, "%s/orphanfixer.params", PARAMS.OUTPUT_DIR);
    fprintf(stderr, "output parameter file to `%s'...", outfname);
    output_params(outfname, &PARAMS);
    fprintf(stderr, "..done\n");

    REDSHIFT = my_malloc(sizeof(*REDSHIFT), NUM_SNAPSHOTS);

#ifdef SUSSING_TREES
    my_snprintf(outfname, MAXLEN, "%s/redshifts.list", PARAMS.GROUP_DIR);
#else
    my_snprintf(outfname, MAXLEN, "%s/redshift", PARAMS.GROUP_DIR);
#endif
    int nred = read_redshifts(outfname, REDSHIFT, NUM_SNAPSHOTS);
    XASSERT(nred == NUM_SNAPSHOTS, "Error: Number of redshifts read in = %d is not equal to the number of snapshots = %d. "\
                                    "Please make sure that '%s' file contains redshifts for all snapshots.\n", nred, NUM_SNAPSHOTS, outfname);

    set_cosmology(&COSMO);
    PARAMS.COSMO = &COSMO;
    set_simulation_params(&PARAMS);

#ifdef FOF_ONLY
    const int fof_only = 1;
#else
    const int fof_only = 0;
#endif

    for (int isnapshot = PARAMS.MIN_SNAPSHOT_NUM; isnapshot <= PARAMS.MAX_SNAPSHOT_NUM; isnapshot++)
        Ngroups[isnapshot] = 0;

    /* WARNING: the definition of Ngroups changes with compilation option FOF_ONLY
     */
    for (int isnapshot = PARAMS.MIN_SNAPSHOT_NUM; isnapshot <= PARAMS.MAX_SNAPSHOT_NUM; isnapshot++)
    {
        t_bigsectionstart = time(NULL);
        fprintf(stderr, "\n\n Now working on snapshot# %4d \n", isnapshot);

        Ngroups0 = returnNhalo(&PARAMS, isnapshot, fof_only);
        if(Ngroups0 == 0) continue;

        if (isnapshot < PARAMS.MAX_SNAPSHOT_NUM)
        {
            parent = my_malloc(sizeof(*parent), Ngroups0);
            my_snprintf(outfname, MAXLEN, "%s/parents_%03d.txt", PARAMS.OUTPUT_DIR, isnapshot);

            fprintf(stderr, "Reading in parents from `%s' Ngroups0 = %" STR_FMT "\n", outfname, Ngroups0);
            t_sectionstart = time(NULL);
            parent = loadparents(outfname, parent, Ngroups0);
            t_sectionend = time(NULL);
            fprintf(stderr, " done ...\n\n");
            print_time(t_sectionstart, t_sectionend, "loadparents");
            allparents[isnapshot] = parent;
        }

        group0 = allocate_group(Ngroups0); // allocate and initialize
        fprintf(stderr, "loading group for snapshot # %d with %" STR_FMT " halos ", isnapshot, Ngroups0);

        t_sectionstart = time(NULL);
        loadgroups(&PARAMS, isnapshot, group0);
        fprintf(stderr, " done ...\n\n");
        t_sectionend = time(NULL);
        print_time(t_sectionstart, t_sectionend, "loadgroups");
        Ngroups[isnapshot] = Ngroups0;

        t_sectionstart = time(NULL);
        my_snprintf(outfname, MAXLEN, "%s/subhalolevel_%03d.txt", PARAMS.OUTPUT_DIR, isnapshot);
        readsubhalo_hierarchy_levels(outfname, group0);
        fprintf(stderr, " done ...\n\n");
        t_sectionend = time(NULL);
        print_time(t_sectionstart, t_sectionend, "subhalo hierarchy levels");

        /*Read in the actual snapshot */
#if (defined(GET_GROUPVEL) == 0 && defined(GET_MEANVEL) == 0)
        char snapshotname[MAXLEN];
        my_snprintf(snapshotname, MAXLEN, "%s/%s_%03d", PARAMS.SNAPSHOT_DIR, PARAMS.SNAPSHOT_BASE, isnapshot);
        t_sectionstart = time(NULL);
        P = loadsnapshot(snapshotname, &header);
        t_sectionend = time(NULL);
        print_time(t_sectionstart, t_sectionend, "loadsnapshot");

        t_sectionstart = time(NULL);
        assign_vxcm(group0, Ngroups0, P, REDSHIFT[isnapshot]);
        t_sectionend = time(NULL);
        print_time(t_sectionstart, t_sectionend, "assign meanvel from snapshot");
        my_free((void **)&P);
#endif

        /* assign the groups for their snapshot number to the correct locations*/
        node = my_malloc(sizeof(*node), Ngroups0);
        t_sectionstart = time(NULL);
        assign_node(group0, Ngroups0, parent, node, isnapshot);
        t_sectionend = time(NULL);
        print_time(t_sectionstart, t_sectionend, "assign_node");

        /* Now store the node data in the tree */
        tree[isnapshot] = node;
        free_group(group0, Ngroups0);

        my_snprintf(outfname, MAXLEN, "%s %d ", "Reading groups and assigning to node for i = ",
                    isnapshot); /*use outfname as temporary message variable. Gets reset..*/
        t_sectionend = time(NULL);
        print_time(t_bigsectionstart, t_sectionend, outfname);
    }

    /* All the groups have now been assigned in the node/tree structure. Make the
     * parent/child pointer associations */
    t_sectionstart = time(NULL);
    maketree(allparents, Ngroups, tree);
    t_sectionend = time(NULL);
    print_time(t_sectionstart, t_sectionend, "make_fof_tree");

    /* Give the halos unique ids starting at redshift 0 and tracing them all the
     * way back. */
    t_sectionstart = time(NULL);
    assign_haloid(tree, Ngroups);
    t_sectionend = time(NULL);
    print_time(t_sectionstart, t_sectionend, "assign_haloid");

#ifndef FOF_ONLY
    if (PARAMS.LOAD_FOUND_PROGENITORS == 1)
    {
        my_snprintf(outfname, MAXLEN, "%s/found_progenitors.txt", PARAMS.OUTPUT_DIR);
        t_sectionstart = time(NULL);
        {
            int64 startgroup = 0;
            load_found_progenitors(tree, Ngroups, outfname, &startgroup); /* tmp_snapshot is useless in this section */
        }
        assign_haloid(tree, Ngroups);
        t_sectionend = time(NULL);
        print_time(t_sectionstart, t_sectionend, "load_found_progenitors");
    }
    else
    {
        t_sectionstart = time(NULL);
        fillprogenitors(tree, Ngroups);
        t_sectionend = time(NULL);
        print_time(t_sectionstart, t_sectionend, "fill progenitors");
    }

    /*   t_sectionstart = time(NULL); */
    /*   output_missing_halo_centres(tree,Ngroups); */
    /*   t_sectionend = time(NULL); */
    /*   print_time(t_sectionstart,t_sectionend,"output missing halo centres"); */

#endif

    /* #ifdef SUSSING_TREES   */
    /*   t_sectionstart = time(NULL); */
    /*   output_data_for_SUSSING_TREES(tree,Ngroups); */
    /*   t_sectionend = time(NULL); */
    /*   print_time(t_sectionstart,t_sectionend,"For SUSSING Merger Tree
     * comparison project"); */
    /* #endif   */

    /* free up the memory */
    for (int isnapshot = PARAMS.MIN_SNAPSHOT_NUM; isnapshot <= PARAMS.MAX_SNAPSHOT_NUM; isnapshot++)
    {
        if (Ngroups[isnapshot] > 0)
        {
            if (isnapshot < PARAMS.MAX_SNAPSHOT_NUM)
                free(allparents[isnapshot]);

            free(tree[isnapshot]);
        }
    }

    free(REDSHIFT);
    free(PARAMS.Age);
    fprintf(stderr, "\n\n Done...\n\n");

    t_codeend = time(NULL);
    print_time(t_codestart, t_codeend, "Entire code");

    exit(EXIT_SUCCESS);
}

void assign_vxcm(struct group_data *group, int64 Ngroups, struct particle_data *P, float redshift)
{
    int64 index;
    double sumx, sumy, sumz, summ;
    double sumvx, sumvy, sumvz;
    float mass = 0.0;
    for (int64 i = 0; i < Ngroups; i++)
    {
        sumx = 0.0;
        sumy = 0.0;
        sumz = 0.0;
        sumvx = 0.0;
        sumvy = 0.0;
        sumvz = 0.0;
        summ = 0.0;

#ifdef FOF_ONLY
        group[i].Mtot = 0.0;
#endif
        for (int64 j = 0; j < group[i].N; j++)
        {
            index = group[i].id[j] - 1; // 0 based indexing
            mass = P[index].Mass;
            sumvx += (P[index].Vel[0] * mass);
            sumvy += (P[index].Vel[1] * mass);
            sumvz += (P[index].Vel[2] * mass);
            summ += mass;
#ifndef FOF_ONLY
            sumx += (periodic(P[index].Pos[0] - group[i].xcen) * mass); /* computing the centre of mass -> needs to
                                                                           account for box wrapping. */
            sumy += (periodic(P[index].Pos[1] - group[i].ycen) * mass);
            sumz += (periodic(P[index].Pos[2] - group[i].zcen) * mass);
#else
            sumx += (periodic(P[index].Pos[0] - group[i].xcen) * mass); /* Make sure that group[i].xcen is initiliased
                                                                           by loadgroups under the FOF_ONLY flag. */
            sumy += (periodic(P[index].Pos[1] - group[i].ycen) * mass);
            sumz += (periodic(P[index].Pos[2] - group[i].zcen) * mass);
#endif

#ifdef FOF_ONLY
            group[i].Mtot += mass;
#endif
        }

#ifdef FOF_ONLY
        group[i].xcen = periodic_wrap(sumx / summ);
        group[i].ycen = periodic_wrap(sumy / summ);
        group[i].zcen = periodic_wrap(sumz / summ);
#else
        group[i].Xoffset = group[i].xcen - periodic_wrap(sumx / summ);
        group[i].Yoffset = group[i].ycen - periodic_wrap(sumy / summ);
        group[i].Zoffset = group[i].zcen - periodic_wrap(sumz / summ);
#endif

        group[i].vxcen = sumvx / summ;
        group[i].vycen = sumvy / summ;
        group[i].vzcen = sumvz / summ;

        group[i].vxcen *= sqrt(get_scalefactor(redshift));
        group[i].vycen *= sqrt(get_scalefactor(redshift));
        group[i].vzcen *= sqrt(get_scalefactor(redshift));
    }
}

void assign_node(struct group_data *group0, int64 Ngroups0, struct parent_data *parent, struct node_data *node,
                 int isnapshot)
{

#ifdef GET_MEANVEL
    FILE *fp = NULL;
    int64 tmp_grpnum = 0;
    char str_line[MAXLINESIZE];
    char fname[MAXLEN];
    const char comment = '#';
    my_snprintf(fname, MAXLEN, "%s/allgroupvels_%03d.txt", PARAMS.GROUP_DIR, isnapshot);
    fp = my_fopen(fname, "rt");
    while (1)
    {
        if (fgets(str_line, MAXLINESIZE, fp) != NULL)
            if (str_line[0] != comment)
                break;
    }

#endif

    struct node_data *FOF_Parent = NULL;
    float scale_factor = 1.0 / (1.0 + REDSHIFT[isnapshot]);
    float root_a = sqrt(scale_factor);
    double boundfofmtot = 0.0;

    for (int64 igroup = 0; igroup < Ngroups0; igroup++)
    {
        if (FOF_Parent != NULL && group0[igroup].isFof == 1)
        {
            /* this will assign the sum of all the subhalo masses as the bound fof
               mtot. however, this will miss the last fof halo if it has only itself
               as the bound subhalo -- need to check for that at the end of the loop
             */
            FOF_Parent->BoundFofMtot = boundfofmtot;
        }

        if (group0[igroup].isFof == 1)
        {
            FOF_Parent = &node[igroup];
            boundfofmtot = group0[igroup].Mtot;
        }
        else
        {
            boundfofmtot += group0[igroup].Mtot;
            node[igroup].BoundFofMtot = group0[igroup].Mtot; /* Not really proper -- subhalos have their
                                                                'normal' mass as the bound fof mass */
        }

        node[igroup].isFof = group0[igroup].isFof;
        /* 	  node[igroup].Nsub     = group0[igroup].Nsub; */
        node[igroup].haloid = -1;
#ifdef SUSSING_TREES
        node[igroup].SUSSING_haloID = group0[igroup].haloID;
#endif

        node[igroup].nodeloc = igroup;
        node[igroup].z = REDSHIFT[isnapshot];
        node[igroup].snapshot = isnapshot;
        node[igroup].xcen = group0[igroup].xcen;
        node[igroup].ycen = group0[igroup].ycen;
        node[igroup].zcen = group0[igroup].zcen;

#ifdef GET_GROUPVEL
        node[igroup].vxcen = group0[igroup].vxcen * root_a;
        node[igroup].vycen = group0[igroup].vycen * root_a;
        node[igroup].vzcen = group0[igroup].vzcen * root_a;
#else
#ifndef GET_MEANVEL
        node[igroup].vxcen = group0[igroup].vxcen; /* group vels have been loaded in from snapshots.
                                                      Copy them to the nodes. */
        node[igroup].vycen = group0[igroup].vycen;
        node[igroup].vzcen = group0[igroup].vzcen;
#endif
#endif

        node[igroup].Mtot = group0[igroup].Mtot;
        node[igroup].ContainerId = group0[igroup].ContainerIndex;
        node[igroup].ParentLevel = group0[igroup].ParentLevel;
        node[igroup].Nsub = group0[igroup].Nsub;
        node[igroup].ContainerHalo = &(node[node[igroup].ContainerId]);

        if ((node[igroup].ContainerHalo)->Nsub == 0)
        {
            fprintf(stderr,
                    "This should not have happened. The container claims to have no "
                    "Nsubs. snapshot = %hd igroup = %" STR_FMT " Container id = %" STR_FMT " ..exiting\n",
                    isnapshot, igroup, node[igroup].ContainerId);
            exit(EXIT_FAILURE);
        }

        if (isnapshot < PARAMS.MAX_SNAPSHOT_NUM)
        {
            if (parent[igroup].parentid >= 0 && parent[igroup].parentsnapshot <= PARAMS.MAX_SNAPSHOT_NUM)
            {
                /* 			  fprintf(stderr,"isnapshot = %d igroup = %d parent
                 * snapshot = %d\n",isnapshot,igroup,parent[igroup].parentsnapshot); */
                node[igroup].ParentID = parent[igroup].parentid;
                node[igroup].ParentZ = REDSHIFT[parent[igroup].parentsnapshot];
                node[igroup].ParentSnapshot = parent[igroup].parentsnapshot;
            }
        }
        else
        {
            node[igroup].ParentID = -1;
            node[igroup].ParentZ = -1.0;
            node[igroup].ParentSnapshot = -1;
        }

        node[igroup].FofHalo = FOF_Parent;
        node[igroup].Sibling = NULL;
        node[igroup].Parent = NULL;
        node[igroup].BigChild = NULL;
        node[igroup].Nchild = 0;
#ifdef GET_MEANVEL
        sscanf(str_line, "%" STR_FMT " %f  %f  %f\n", &tmp_grpnum, &(node[igroup].meanvel[0]),
               &(node[igroup].meanvel[1]), &(node[igroup].meanvel[2]));
        if (tmp_grpnum != igroup)
        {
            fprintf(stderr, "Error: While reading in velocites for the groups\n");
            fprintf(stderr, "expected groupnum = %" STR_FMT "  got groupnum = %" STR_FMT " ..exiting\n", igroup,
                    tmp_grpnum);
            exit(EXIT_FAILURE);
        }

        if (fgets(str_line, MAXLINESIZE, fp) == NULL && igroup != (Ngroups0 - 1))
        {
            fprintf(stderr,
                    "Error (eof): expected more lines to read in from group "
                    "velocities file `%s' \n",
                    fname);
            fprintf(stderr, "currently at igroup = %" STR_FMT " ..exiting\n", igroup);
            exit(EXIT_FAILURE);
        }

        /*
                  MS: 5th Aug, 2010. Took out the Gadget sqrt(a) dependence of v
        */

        for (int k = 0; k < 3; k++)
            node[igroup].meanvel[k] *= root_a;

#else
#ifdef GET_GROUPVEL
        node[igroup].meanvel[0] = node[igroup].vxcen;
        node[igroup].meanvel[1] = node[igroup].vycen;
        node[igroup].meanvel[2] = node[igroup].vzcen;
#endif
#endif

        /* Make sure that the last fof halo also has the bound fof halo mass
         * assigned. */
        if (igroup == (Ngroups0 - 1) && node[igroup].isFof == 1)
        {
            node[igroup].BoundFofMtot = boundfofmtot;
        }
    }
#ifdef GET_MEANVEL
    fclose(fp);
#endif
}

void readsubhalo_hierarchy_levels(const char *fname, struct group_data *group0)
{
    FILE *fp = NULL;
    char str_line[MAXLINESIZE];
    const char comment = '#';
    int64 i = 0, dummy;
    short dummy1;
    fp = my_fopen(fname, "r");
    i = 0;
    while (1)
    {
        if (fgets(str_line, MAXLINESIZE, fp) != NULL)
        {
            if (str_line[0] != comment)
            {
                /*                    sscanf(str_line,"%*"STR_FMT" %"STR_FMT" %hd
                 * %"STR_FMT"  %"STR_FMT"     %hd",  */
                sscanf(str_line, "%*d   %" STR_FMT " %hd   %" STR_FMT "  %" STR_FMT "     %hd", &dummy,
                       &(group0[i].ParentLevel), &(group0[i].ContainerIndex), &(group0[i].Nsub),
                       &dummy1); /*discarding the first field, Fofid */
                // The first field in the hierarchy file is of type int64. But since I
                // am suppressing that field -> I can as well use the %d type [and this
                // way the compiler warning goes away.

                if (dummy != i)
                {
                    fprintf(stderr,
                            "Error: While reading in file `%s' for subhalo hierarchy "
                            "level \n ",
                            fname);
                    fprintf(stderr, "expected igroup = %" STR_FMT " instead got %" STR_FMT " ..exiting \n\n", i, dummy);
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
}

#ifdef SUSSING_TREES
void output_data_for_SUSSING_TREES(struct node_data *tree[], int64 *Ngroups)
{
    struct node_data *BaseNode = NULL, *thisnode = NULL;
    FILE *fp = NULL;
    char outfname[MAXLEN];
    int64_t SumTotal_AllHalos = 0;
    int64_t Nbad = 0;

    char date[9];
    time_t t = time(0);
    struct tm *tm;
    tm = gmtime(&t);
    strftime(date, sizeof(date), "%Y%m%d", tm);

    for (int isnapshot = PARAMS.MAX_SNAPSHOT_NUM; isnapshot >= PARAMS.MIN_SNAPSHOT_NUM; isnapshot--)
    {
        SumTotal_AllHalos += Ngroups[isnapshot];
    }

    my_snprintf(outfname, MAXLEN, "%s/HINGE_Sinha_DATASET-II.txt", PARAMS.OUTPUT_DIR);
    fp = my_fopen(outfname, "w");
    fprintf(fp, "1\n");
    fprintf(fp,
            "Halo Interaction Network and Galaxy Evolution - settings used from "
            "previous experience with SUBFIND. v1.0, %s \n",
            date);
    fprintf(fp, "%" PRId64 "\n", SumTotal_AllHalos);

    for (int isnapshot = PARAMS.MAX_SNAPSHOT_NUM; isnapshot >= PARAMS.MIN_SNAPSHOT_NUM; isnapshot--)
    {
        BaseNode = tree[isnapshot];
        for (int64 igroup = 0; igroup < Ngroups[isnapshot]; igroup++)
        {
            thisnode = &BaseNode[igroup];
            fprintf(fp, "%ld  \t %" STR_FMT " \n", thisnode->SUSSING_haloID, thisnode->Nchild);
            if (thisnode->isFof == 0 && thisnode->BigChild == NULL && thisnode->FofHalo->BigChild != NULL)
                Nbad++;

            thisnode = thisnode->BigChild;
            while (thisnode != NULL)
            {
                fprintf(fp, "%ld \n", thisnode->SUSSING_haloID);
                thisnode = thisnode->Sibling;
            }
        }
    }

    fprintf(fp, "END\n");
    fclose(fp);
    fprintf(stderr, "Nbad = %" PRId64 "\n", Nbad);
}
#endif

void print_makefile_options(void)
{

#ifdef FOF_ONLY
    fprintf(stderr, "The code is going to read in FOF groups only\n");
#else
    fprintf(stderr, "The code is going to read in Subfind groups \n");
#endif

#ifdef WMAP1
    fprintf(stderr, "The code was compiled for WMAP 1 parameters \n");
#endif

#ifdef WMAP3
    fprintf(stderr, "The code was compiled for WMAP 3 parameters \n");
#endif

#ifdef WMAP5
    fprintf(stderr, "The code was compiled for WMAP 5 parameters \n");
#endif

#ifdef GET_GROUPVEL
    fprintf(stderr, "The code was compiled to load groups with velocity info\n");
#else
    fprintf(stderr, "The code was *NOT* compiled to read in velocities from the "
                    "group files \n");
#endif

#ifdef GET_MEANVEL
    fprintf(stderr, "The code was compiled to read allgroupvels_???.txt files "
                    "for velocity info\n");
#else
    fprintf(stderr, "The code was *NOT* compiled to read allgroupvels_???.txt "
                    "files for velocity info\n");
#endif

#ifdef LONGIDS
    fprintf(stderr, "Assumes that Gadget particle ids are 8 bytes (but the total "
                    "particle number is inside INT_MAX)\n");
#else
    fprintf(stderr, "Assumes that Gadget particle ids are 4 bytes (and the total "
                    "particle number is inside INT_MAX)\n");
#endif

#ifdef BIGSIM
    fprintf(stderr, "Assumes particle load is larger than INT_MAX (2 billion)\n");
#else
    fprintf(stderr, "Assumes particle load is smaller than INT_MAX (2 billion)\n");
#endif

#ifdef SUSSING_TREES
    fprintf(stderr, "The code will assume data for the SUSSING Mergertree "
                    "Comparison Project\n");
#endif

#ifdef ASCII_DATA
    fprintf(stderr, "The code will read in ASCII input data (only valid with "
                    "-DSUSSING_TREES; ignored otherwise) \n");
#endif
}
