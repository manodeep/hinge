#pragma once

#include <stdio.h>
#include <stdlib.h>

#include "hinge.h" //for struct cosmology_data definition
#include "io.h"    //for struct group_data definition

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
