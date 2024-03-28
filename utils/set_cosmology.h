#pragma once

#include <stdio.h>
#include <stdlib.h>

#include "io.h"

struct cosmology_data
{
    double G;
    double H0;
    double h100;
    double Omega_0;
    double Omega_b;
    double Omega_dm;
    double Omega_m;
    double Omega_lambda;
    double Omega_rad;
    double Omega_k;
    float K_fit; /* Maccio et al. 2008 fitting parameters*/
    float F_fit;
};

/* functions in set_cosmology */
double agefunc(double z, void *params);
double epeebles(float z);
double get_age(float z);
void set_cosmology(struct cosmology_data *CP);
float get_scalefactor(float z);
float get_overdensity(float z, struct cosmology_data *CP);
double get_rhocrit_at_redshift(float z, struct cosmology_data *CP);
double get_hubble_at_redshift(float z, struct cosmology_data *CP);
double getrvir_anyl(const double mvir, const float z, struct cosmology_data *CP);
float getconc_anyl(const double mvir, const float z);
void getrvir_from_overdensity(struct group_data *group, int NBINS, const double RhoCrit, const float OverDensity);
