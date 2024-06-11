#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

#include "defs.h"
#include "io.h"

struct hinge_halo
{
    int64_t halo_id;
    int64_t fof_id;
    int64_t nsub;
    int64_t npart;
    double Mvir;
    double Xc;
    double Yc;
    double Zc;
    double VXc;
    double VYc;
    double VZc;
    double Rvir;
};

struct hinge_catalog
{
    int64_t nallocated;
    int64_t nhalos;
    int64_t nfofs;
    int64_t totnpart;
    struct hinge_halo *halos;
};


    extern struct hinge_catalog* read_hinge_ascii_halo_catalog(const char *fname, const int fof_only) __attribute__((warn_unused_result));
    extern void free_hinge_halocat(struct hinge_catalog *halocat);

#ifdef __cplusplus
}
#endif