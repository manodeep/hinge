/*
This code will find the parents of all halos at a given redshift.
(The parents are at future redshifts). One snapshot may be skipped
in this quest for an unique parent. The idea is that the parent in
uniquely identified by the value [for a possible j'th parent halo]

          R_j = \sum_k BE_k^(-2/3),   [Ref: Boylan-Kolchin, Millenium
Simulation-II paper]

where k runs over all the common elements between the halo in question
and the would be j'th parent halo. BE is the binding energy rank, we
will proxy it by the knowledge that subfind outputs all the particle
positions already sorted by potential. In this formulation, more highly
bound particles get a higher precedence when determining the unique
parent halo.

While I have coded in the compatibility for only FOF halos, the usefulness
is quite doubtful.


Log started: July 2009: MS.

07/01/2009:   started writing the codes. Complete with a Makefile.

07/19/2009:   loadgroups now written and compiles. Needs to be cross-matched
              with known-good IDL routines to make sure that the data is read
                          in correctly.

03/14/2011:  Wow!! The changes to this code are really few and far between. Now,
I am trying to identify cases where subhalos might be `lost' in the host halo,
i.e, the subhalo appears to dissolve in the FOF [the code should really catch
all escalation of hierarchy levels] -- in such a case the subhalo does not get
assigned a parent. If the subhalo re-appears within the next two snapshots, then
it will be assigned properly -- thereby reducing the load on the mergertree
code.



Important global variable [macro really] definition -- GROUPMINLEN. This sets
the minimum group length for a `legitimate' group. It's probably safe to set
this to around 100 and have the groupfinder identify down to 20-30 particles.
Reduces the clutter at the low mass parent matching.


*/

#include <assert.h>
#include <inttypes.h> //defines PRId64 for printing int64_t
#include <limits.h>
#include <math.h>
#include <stdint.h> //defines int64_t datatype -> *exactly* 8 bytes int
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "findallparents.h"
#include "findprogenitor.h"
#include "hierarchy.h"
#include "hinge.h"
#include "io.h"
#include "loadgroups.h"
#include "read_param.h"
#include "switchfof.h"
#include "utils.h"
#include "utils_read_params.h" //for global variable extern definition + function calls

struct params_data PARAMS; // global variable
float *REDSHIFT;

// functions in main
void print_makefile_options(void);

int main(int argc, char **argv)
{

#ifdef FOF_ONLY
    const int fof_only = 1
#else
    const int fof_only = 0;
#endif

        FILE *fd = NULL;
    int64 Ngroups0 = 0;
    int64 Ngroups1 = 0;
    int64 NFof0 = 0;
    /* int64 NFof1 = 0; */

    int snapshot_number;
    char outfname[MAXLEN];
    int64 notfound = 0;
    int64 Nparentsfound = 0;
    int incr = 1;

    struct group_data *group0 = NULL, *group1 = NULL;
    /* struct io_header header; */
    int NUM_SNAPSHOTS;

    time_t t_codestart, t_codeend, t_sectionstart, t_sectionend, t_bigsectionstart;
    t_codestart = time(NULL);

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

#if ((defined(SUSSING_TREES) + defined(ASCII_DATA) + defined(BGC2) + defined(SUBFIND)) > 1)
#error Only ONE of the MAKEFILE options SUSSING_TREES, ASCII_DATA, BGC2 should be selected
#endif

    // Check command line for the parameter file name.
    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s <parameterfile> <snapshotnumber> \n", argv[0]);
        fprintf(stderr, "Code was compiled with \n");
        print_makefile_options();
        exit(EXIT_FAILURE);
    }

    my_snprintf(outfname, MAXLEN, "%s", argv[1]); // parameter file

    // get the snapshot number
    snapshot_number = atoi(argv[2]);

    // print out the actual commandline arguments used
    fprintf(stderr, "\n\n Running `%s' with the following parameters \n", argv[0]);
    fprintf(stderr, "\t\t\t parameter file   = `%s'\n", outfname);
    fprintf(stderr, "\t\t\t snapshot number  = `%d'\n", snapshot_number);
    fprintf(stderr, "\n\n");

    // read in the parameter file
    fprintf(stderr, "reading parameter file `%s'...", outfname);
    read_params(outfname, &PARAMS, haloparentfinder_fill_params);
    fprintf(stderr, "..done\n");

    // Check if the snapshot is the last one - nothing to do then
    if (snapshot_number >= PARAMS.MAX_SNAPSHOT_NUM)
    {
        fprintf(stderr,
                "snapshot_number = %d is greater than or equal to the last "
                "snapshot (%d). Exiting since there is nothing to do \n",
                snapshot_number, PARAMS.MAX_SNAPSHOT_NUM);
        exit(EXIT_SUCCESS);
    }

    NUM_SNAPSHOTS = PARAMS.MAX_SNAPSHOT_NUM + 1;
    REDSHIFT = my_malloc(sizeof(*REDSHIFT), NUM_SNAPSHOTS);

#ifdef SUSSING_TREES
    my_snprintf(outfname, MAXLEN, "%s/redshifts.list", PARAMS.GROUP_DIR);
#else
    my_snprintf(outfname, MAXLEN, "%s/redshift", PARAMS.GROUP_DIR);
#endif
    int nred = read_redshifts(outfname, REDSHIFT, NUM_SNAPSHOTS);

    // output the parameter file
    my_snprintf(outfname, MAXLEN, "%s/fofmatch.params", PARAMS.OUTPUT_DIR);
    fprintf(stderr, "output parameter file to `%s'...", outfname);
    output_params(outfname, &PARAMS);
    fprintf(stderr, "..done\n");

    NFof0 = returnNhalo(&PARAMS, snapshot_number, 1);
    fprintf(stderr, "NFof0 = %" STR_FMT " outfname = `%s'\n", NFof0, outfname);

    /*Now to actually allocate the memory and load the halos */
    Ngroups0 = returnNhalo(&PARAMS, snapshot_number, fof_only);

    notfound = NFof0;
    if (Ngroups0 > 0)
    {
        group0 = allocate_group(Ngroups0);
        fprintf(stderr, "loading group for snapshot # %d with %" STR_FMT " halos\n", snapshot_number, Ngroups0);
        t_sectionstart = time(NULL);
        loadgroups(&PARAMS, snapshot_number, group0);
        t_sectionend = time(NULL);
        fprintf(stderr, " done ...\n\n");
        print_time(t_sectionstart, t_sectionend, "loadgroups");
        Ngroups1 = returnNhalo(&PARAMS, snapshot_number + incr, fof_only);
        group1 = allocate_group(Ngroups1);
        fprintf(stderr, "loading group for snapshot # %d with %" STR_FMT " halos\n", snapshot_number + incr, Ngroups1);
        t_sectionstart = time(NULL);
        loadgroups(&PARAMS, snapshot_number + incr, group1);
        t_sectionend = time(NULL);
        fprintf(stderr, " done ...\n\n");
        print_time(t_sectionstart, t_sectionend, "loadgroups");

        /* #if !((defined(SUSSING_TREES)) || (defined(AHF_INPUT))) */
        fprintf(stderr,
                "find hierarchy level for the subhalos group for snapshot # %d "
                "with %" STR_FMT " halos\n",
                snapshot_number, Ngroups0);
        t_sectionstart = time(NULL);
        find_hierarchy_level(group0, Ngroups0, PARAMS.OUTPUT_DIR);
        t_sectionend = time(NULL);
        print_time(t_sectionstart, t_sectionend, "hierarchy level at current snapshot ");
        /* #endif		 */

#ifdef MAKE_LEAN
        fprintf(stderr, "freeing memory associated with particle positions \n");
        free_group_positions(group0, Ngroups0);
#endif

        /* #if !defined(SUSSING_TREES) && !defined(AHF_INPUT) */
        fprintf(stderr,
                "find parent level for the subhalos group for snapshot # %d with "
                "%" STR_FMT " halos\n",
                snapshot_number + incr, Ngroups1);
        t_sectionstart = time(NULL);
        find_hierarchy_level(group1, Ngroups1, PARAMS.OUTPUT_DIR);
        t_sectionend = time(NULL);
        print_time(t_sectionstart, t_sectionend, "hierarchy level at next snapshot ");
        /* #endif */

#ifdef MAKE_LEAN
        fprintf(stderr, "freeing memory associated with particle positions \n");
        free_group_positions(group1, Ngroups1);
#endif

        t_sectionstart = time(NULL);
        Nparentsfound = findfofparents(group0, Ngroups0, group1, Ngroups1, PARAMS.OUTPUT_DIR);
        t_sectionend = time(NULL);
        print_time(t_sectionstart, t_sectionend, "Find FOF parents");

        fprintf(stderr, "Nparentsfound = %" STR_FMT " NFof0 = %" STR_FMT "\n", Nparentsfound, NFof0);
        notfound = NFof0;
        t_sectionstart = time(NULL);
        for (int64 i = 0; i < Ngroups0; i++)
        {
            if (group0[i].ParentId >= 0 && group0[i].isFof == 1)
                notfound--;
        }
        t_sectionend = time(NULL);
        print_time(t_sectionstart, t_sectionend, "Fill Fof ranks");

        if (notfound == 0)
        {
            /* done finding FOF parents */
            fprintf(stderr,
                    "\nAll FOF halos at %4d snapshot have been assigned to other FOF "
                    "halos at snapshot %4d \n",
                    snapshot_number, snapshot_number + incr);
        }
        else
        {
            /* some of the FOF halos at snapshot have not been assigned to other FOF
             * halos at snapshot+1 */
            fprintf(stderr, "\nOut of the %6" STR_FMT " FOF halos, %6" STR_FMT " halos could not be assigned.\n", NFof0,
                    NFof0 - Nparentsfound);
        }

        notfound = Ngroups0 - Nparentsfound;
        incr = 1;
        t_bigsectionstart = time(NULL);
        while (notfound > 0 && snapshot_number + incr <= PARAMS.MAX_SNAPSHOT_NUM && incr <= PARAMS.MAX_INCR)
        {
            if (abs(incr) != 1)
            {
                fprintf(stderr,
                        "\n\n %" STR_FMT " halos could not be assigned parents. Trying "
                        "to check the snapshot %d now \n",
                        notfound, snapshot_number + incr);
                free_group(group1, Ngroups1);

                Ngroups1 = returnNhalo(&PARAMS, snapshot_number + incr, fof_only);

                fprintf(stderr, "\nNow looking for subhalo parents \n");
                group1 = allocate_group(Ngroups1);
                loadgroups(&PARAMS, snapshot_number + incr, group1);

                /* #if !defined(SUSSING_TREES) && !defined(AHF_INPUT) */
                fprintf(stderr,
                        "find hierarchy level for the subhalos group for snapshot # %d "
                        "with %" STR_FMT " halos\n",
                        snapshot_number + incr, Ngroups1);
                t_sectionstart = time(NULL);
                find_hierarchy_level(group1, Ngroups1, PARAMS.OUTPUT_DIR);
                t_sectionend = time(NULL);
                /* 			  fprintf(stderr," done ...\n\n"); */
                print_time(t_sectionstart, t_sectionend, "hierarchy level at next snapshot ");
                /* #endif	 */

#ifdef MAKE_LEAN
                fprintf(stderr, "freeing memory associated with particle positions \n");
                free_group_positions(group1, Ngroups1);
#endif
            }

            Nparentsfound = findallparents(group0, Ngroups0, group1, Ngroups1, (const int)snapshot_number + incr,
                                           PARAMS.OUTPUT_DIR);
            notfound = Ngroups0 - Nparentsfound;

            if (incr == 1)
            {
                t_sectionstart = time(NULL);
                check_fof_matches(group0, Ngroups0, group1, Ngroups1, (const int)snapshot_number + incr,
                                  PARAMS.OUTPUT_DIR); // in switchfof.c
                t_sectionend = time(NULL);
                print_time(t_sectionstart, t_sectionend, "finished check_fof_matches  ");

                /* Also, check for subhaloes at the next snapshot that do not have
                 * progenitors */
                t_sectionstart = time(NULL);
                find_progenitor(group1, Ngroups1, group0, Ngroups0, PARAMS.OUTPUT_DIR);
                t_sectionend = time(NULL);
                print_time(t_sectionstart, t_sectionend, "finished findprogenitor  ");
            }

            incr++;
            fprintf(stderr, "incr = %d\n", incr);
        }

        if (notfound > 0)
        {
            if (incr > PARAMS.MAX_INCR)
                fprintf(stderr,
                        "\n\n %" STR_FMT " halos could not be assigned even after "
                        "looking at snapshot %d ..Giving up\n\n",
                        notfound, snapshot_number + incr);

            if ((snapshot_number + incr) > PARAMS.MAX_SNAPSHOT_NUM)
            {
                fprintf(stderr, "\n %" STR_FMT " halos are missing parents ", notfound);
                fprintf(stderr, "\n But there are no further snapshots to "
                                "load..ignoring these halos\n");
            }
        }
        else
        {
            fprintf(stderr,
                    "\n All the halos with missing parents have now been assigned "
                    "parents by snapshot %d \n",
                    snapshot_number + incr);
        }
        t_sectionend = time(NULL);
        print_time(t_bigsectionstart, t_sectionend, "While Loop finding all parents with skipping ");

        /* All halos have been found or some limiting condition has been met. Output
         * the data */
        t_sectionstart = time(NULL);
        my_snprintf(outfname, MAXLEN, "%s/parents_%03d.txt", PARAMS.OUTPUT_DIR, snapshot_number);
        fprintf(stderr, "\n\n\nopening parents file `%s' \n\n\n", outfname);
        fd = my_fopen(outfname, "w");
        fprintf(fd, "##############################################################"
                    "##############################################################"
                    "###################################\n");
        fprintf(fd, "# Snapshot    ThisGroupId    ParentLevel   ContainerId      "
                    "Nsub      NpartinHalo         ParentId   ParentSnap      "
                    "NpartinParent      Rank          NCommon  \n");
        fprintf(fd, "#   i             l              i              l            "
                    "l           l                    l          i                "
                    "l               d              l     \n");
        fprintf(fd, "##############################################################"
                    "##############################################################"
                    "###################################\n");
        for (int64 i = 0; i < Ngroups0; i++)
        {
            fprintf(fd,
                    "%4d  %14" STR_FMT " %14d  %10" STR_FMT "   %14" STR_FMT "  %14" STR_FMT "  %14" STR_FMT
                    " %10d        %14" STR_FMT " %14.4lf  %12" STR_FMT " \n",
                    snapshot_number, i, group0[i].ParentLevel, group0[i].ContainerIndex, group0[i].Nsub, group0[i].N,
                    group0[i].ParentId, group0[i].ParentSnapshot, group0[i].NpartinParent, group0[i].Rank,
                    group0[i].Ncommon);
            if (group0[i].ParentId < 0 && group0[i].Nsub > 1)
            {
                fprintf(stderr, "WARNING: Missing parent halo properties: \n");
                fprintf(stderr,
                        "Snapshot = %d  groupid = %" STR_FMT "  isfof = %d   npartinhalo = %" STR_FMT
                        "  Nsub = %" STR_FMT " \n",
                        snapshot_number, i, group0[i].isFof, group0[i].N, group0[i].Nsub);
            }
        }
        fclose(fd);

        t_sectionend = time(NULL);
        fprintf(stderr, " done ...\n\n");
        print_time(t_sectionstart, t_sectionend, "Writing out parents file");

        free_group(group0, Ngroups0);
        free_group(group1, Ngroups1);

        fprintf(stderr, "\n\nCompleted assigning groups at snapshot # %d successfully\n\n", snapshot_number);
        t_codeend = time(NULL);
        print_time(t_codestart, t_codeend, "Entire code");
    }
    else
    {
        fprintf(stderr, "There are no groups to match ..exiting \n");
    }
    free(REDSHIFT);

    return EXIT_SUCCESS;
}

void print_makefile_options(void)
{

#if 0
#ifdef FOF_ONLY
    fprintf(stderr, "The code is going to read in FOF groups only\n");
#endif

#ifdef SUBFIND
    fprintf(stderr, "The code is going to read in Subfind groups \n");
#endif

#ifdef SUSSING_TREES
    fprintf(stderr, "The code will assume data for the SUSSING Mergertree "
                    "Comparison Project\n");
#endif

#ifdef ASCII_DATA
    fprintf(stderr, "The code will read in ASCII input data (only valid with "
                    "-DSUSSING_TREES; ignored otherwise) \n");
#endif
#endif

#ifdef GET_GROUPVEL
    fprintf(stderr, "The code is going to read in velocities from the group files\n");
#else
    fprintf(stderr, "The code is *NOT* going to read in velocities from the group files\n");
#endif

#ifdef LONGIDS
    fprintf(stderr, "Assumes that Gadget particle ids are 8 bytes\n");
#else
    fprintf(stderr, "Assumes that Gadget particle ids are 4 bytes\n");
#endif

#ifdef BIGSIM
    fprintf(stderr, "Assumes particle load is larger than INT_MAX (~2 billion)\n");
#else
    fprintf(stderr, "Assumes particle load is smaller than INT_MAX (~2 billion)\n");
#endif

#ifdef MAKE_LEAN
    fprintf(stderr, "The code will free particle positions after the hierarchies "
                    "are found\n");
#endif
}
