#pragma once

#include "hinge.h"
#include "io.h"
#include "loadparents.h"

struct node_data
{
    int64 nodeloc;
    int64 haloid;
    float z;
    short snapshot;

    float xcen;
    float ycen;
    float zcen;

    float vxcen;
    float vycen;
    float vzcen;

    float meanvel[3]; /*all velocities are assumed to be "Gadget" snapshot
                         velocities and converted to peculiar velocities by
                         multiplying with sqrt(a) */

#ifdef SUSSING_TREES
    long SUSSING_haloID;
#endif

    double Mtot;
    double BoundFofMtot;
    double InfallMass;
    short InfallSnapshot;
    double Mstar; // MS - added on 6th Oct, 2011. Assigned using SAM models
    float FormationRedshift;
    float DestructionRedshift;
    float RedshiftofLastMerger;
    int64 Nmergers;
    int64 TotNmergers;

    int64 NDisruptions;
    int64 TotNDisruptions;

    int64 NDissolutions;
    int64 TotNDissolutions;

    int64 NFlybys;
    int64 TotNFlybys;

    /* Variables for the density profile fit */
    double Rvir_anyl; /* analytically obtained Rvir for some overdensity given
                         Mtot and z (assuming spherical halo) */
    double Rvir;      /* Numerically obtained Rvir by computing the density profile and
                         locating the overdensity radius. */
    double Rhalf;     /* Numerically obtained value for half mass radius*/
    float MaxOverDensity;
    float OverDensityThresh;
    float Conc;
    double Rs;

    double RVmax;
    double Vmax; // km/s
#ifndef FOF_ONLY
    double Xoffset; /*offset between potential center and density center*/
    double Yoffset;
    double Zoffset;
#endif
    int64 ParentID; /* Group number for the parent group at a following snapshot */
    float ParentZ;
    short ParentSnapshot;

    short isFof;
    short ParentLevel; /* misnomer. means the hierarchy level */
    int64 ContainerId; /* Group number that contains this halo. For a Fof, this
                          will be the Fof itself */

#ifndef FOF_ONLY
    int64 Nsub; /* actual number of subhalos. Nsub equals to subhalos for the Fof,
                   and the number of (sub)_n-subs for subhalos */
#endif
    int64 Nchild; /* actual number of pointers that consider this halo to be a
                     parent */

    short VisitedForMassloss;
    short VisitedForCumulativeMerger;

    unsigned int CumulativeNmergers; /* Not sure why I have unsigned int ??!!
                                        int64 may be?*/
    unsigned int CumulativeNFlybys;

    struct node_data *FofHalo; /* FOF container at this snapshot */
    struct node_data *Sibling;
    struct node_data *Parent;   /* The halo this is going to go into in a future snapshot.
                                   Depends on the parent matching algorithm - mostly at the
                                   immediately next snapshot */
    struct node_data *BigChild; /* Progenitor at a previous snapshot. Similar to Parent (mostly
                                   the immediately previous snapshot but not necessarily so ) */

    struct node_data *ContainerHalo; /* container halo; for a FOF/sub-halo,
                                        container is the Fof; a sub-subhalo will
                                        have a container as a subhalo and so on*/

    /*  node** at some point make an array with all the children, i.e., all halos
     that point to this halo as a Parent. */

    /* Black hole related stuff */
    short Seeded;
    int Nbh;
    struct bh_data *BH;

    /*environment related variables - used in environment.c */
    short environment;
    int64 N_Neighbours;
    /*   struct node_data *Neighbours; */
    struct node_data *(*Neighbours)[]; // So this is a  pointer to array of pointers to struct
                                       // node_data (as explained by cdecl)
};

double assign_stellar_mass_from_mvir(struct node_data *const thisnode, int model);
void maketree(struct parent_data *allparents[], int64 *Ngroups, struct node_data *tree[]);
void assign_haloid(struct node_data *tree[], int64 *Ngroups);
void find_mergers(struct node_data *tree[], int64 *Ngroups);
void assign_parent(struct node_data *halo, int64 haloid, struct node_data *parenthalo, int64 parentid);
void find_subsub_mergers(struct node_data *tree[], int64 *Ngroups);
int get_tag(struct node_data *const node1, struct node_data *const node2, float *z, float *rnorm, int *destruction_snap,
            double *submtot, double *parentmtot, double *parentmstar);
void ensure_same_snapshot_for_halos(struct node_data **node1, struct node_data **node2);
float get_separation_between_centres(struct node_data *g, struct node_data *f);
int halo_is_always_fof(struct node_data *fof);
