#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

#include "defs.h"

    struct group_data
    {
        id64 *id; // the exact data type depends on Makefile option -DLONGIDS
        float *x;
        float *y;
        float *z;
        float xcen;
        float ycen;
        float zcen;
        long haloID;
        long fofID;
#ifdef SUSSING_TREES
        float *ParticleEnergy;
        float *vx;
        float *vy;
        float *vz;
#endif

        /*   float minpos[3]; */
        /*   float maxpos[3]; */
        int64 groupnum; // the exact data type depends on Makefile option -DBIGSIM
        float vxcen;
        float vycen;
        float vzcen;

        short snapshot;
        float redshift;
        int64 nodeloc;

        /*   short is_wrapped[3]; */
        /*   short is_straddled[3]; */
        /*   short is_lower[3]; */
        /*   short is_upper[3]; */
        /*   short bin_wrapped[3]; */

        /* Group specific data that need to be assigned from the groups file*/
        double Mtot;
        double Rvir;
        double Rvir_anyl;
        double Rmax;

        short isFof; /* Originally here to show FOF halo or subhalo.
                                                superceded by the next two fields. */

        /*
               the following two fields store the hierarchy for substructure.
               ParentLevel is 1 for FOF halo, 2 for subhalo, 3 for sub-sub halo
               and so on. ContainerIndex contains the index for the container
               group. FOF halos will have a container index same as their own
               index, subhalos will have the index of their FOF halos, sub-subs
               will have the index of their subhalos and so on..

               In case you are using a different kind of group-finder, it is
               imperative that you identify the ParentLevel and the ContainerIndex
               correctly -- they are very important for the parent matching bit.
               Look in findallparents.c for the sequence in which the ParentLevel
               comes into play.
         */

        short int ParentLevel; // hierarchy level of the halo
        int64 ContainerIndex;  //

        int64 N;       // number of particles in the group
        int64 Nsub;    // number of [sub-] subhalos [subhalos will have the number
                       // sub-subhalos listed, and so on]
        int64 FOFHalo; /* FOF halo container for this group*/

        /* Variables for parent matching at the group level */
        int64 ParentId;
        int ParentSnapshot;
        int64 NParents;
        int64 Ncommon;
        double Rank;
        int64 NpartinParent;
        short Switched;

        /* Variables for parent matching at the particle level */
        int64 *parentgroupforparticle;
        int *parentsnapshotforparticle;
        short N_per_wedge;

#ifndef FOF_ONLY
        double Xoffset; /*offset between potential center  and density center*/
        double Yoffset;
        double Zoffset;
#endif

        /* Variables for the density profile fit */
        double Rhalf;
        float MaxOverDensity;
        float OverDensityThresh;
        float Conc;
        double Rs;

        double RVmax;
        double Vmax; // kms/s
    };

    struct io_header
    {
        int32_t npart[6];               /*!< number of particles of each type in this file */
        double mass[6];                 /*!< mass of particles of each type. If 0, then the masses are
                                           explicitly stored in the mass-block of the snapshot file,
                                           otherwise they are omitted */
        double time;                    /*!< time of snapshot file */
        double redshift;                /*!< redshift of snapshot file */
        int32_t flag_sfr;               /*!< flags whether the simulation was including star
                                           formation */
        int32_t flag_feedback;          /*!< flags whether feedback was included (obsolete) */
        uint32_t npartTotal[6];         /*!< total number of particles of each type in this
                                           snapshot. This can be different from npart if one
                                           is dealing with a multi-file snapshot. */
        int32_t flag_cooling;           /*!< flags whether cooling was included  */
        int32_t num_files;              /*!< number of files in multi-file snapshot */
        double BoxSize;                 /*!< box-size of simulation in case periodic boundaries were
                                           used */
        double Omega0;                  /*!< matter density in units of critical density */
        double OmegaLambda;             /*!< cosmological constant parameter */
        double HubbleParam;             /*!< Hubble parameter in units of 100 km/sec/Mpc */
        int32_t flag_stellarage;        /*!< flags whether the file contains formation times
                                           of star particles */
        int32_t flag_metals;            /*!< flags whether the file contains metallicity values
                                           for gas and star particles */
        uint32_t npartTotalHighWord[6]; /*!< High word of the total number of
                                           particles of each type */
        int32_t flag_entropy_instead_u; /*!< flags that IC-file contains entropy
                                           instead of u */
        char fill[60];                  /*!< fills to 256 Bytes */
    };

#define MAXTAGLEN 50

    extern const char GROUP_FORMAT_NAMES[][MAXTAGLEN];
    extern const enum valid_group_formats GROUP_FORMAT_ENUMS[];
    extern const int nvalid_group_format_names;

#define SQR_PERIODIC(dx) (periodic(dx) * periodic(dx))

#ifdef __cplusplus
}
#endif
