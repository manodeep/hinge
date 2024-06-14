#include "set_cosmology.h"
#include "hinge.h"
#include "io.h"
#include "read_param.h"
#include "utils.h"
#include "macros.h"
#include <gsl/gsl_integration.h>

void set_cosmology(struct cosmology_data *CP)
{

    fprintf(stderr, "Setting cosmology parameters ....");
    CP->G = 43007.1; /* km/s kpc 1e10 Msun*/

#ifdef WMAP5
    CP->H0 = 71.9 / 1e3; /* in km/s/kpc */
    CP->Omega_0 = 1.0;
    CP->Omega_b = 0.044;
    CP->Omega_dm = 0.214;
    CP->Omega_m = CP->Omega_b + CP->Omega_dm;
    CP->Omega_lambda = CP->Omega_0 - CP->Omega_m;
    CP->Omega_rad = 0.0;
    CP->Omega_k = 0.0;
    CP->K_fit = 3.6; /* fit for overdensity of 200 and all halos. not just relaxed
                        ones. */
    CP->F_fit = 0.01;
#endif

#ifdef WMAP3
    CP->H0 = 71.9 / 1e3; /* in km/s/kpc */
    CP->Omega_0 = 1.0;
    CP->Omega_b = 0.044;
    CP->Omega_dm = 0.214;
    CP->Omega_m = CP->Omega_b + CP->Omega_dm;
    CP->Omega_lambda = CP->Omega_0 - CP->Omega_m;
    CP->Omega_rad = 0.0;
    CP->Omega_k = 0.0;
    CP->K_fit = 3.4;
    CP->F_fit = 0.01;
#endif

#ifdef WMAP1
    CP->H0 = 71.9 / 1e3; /* in km/s/kpc */
    CP->Omega_0 = 1.0;
    CP->Omega_b = 0.044;
    CP->Omega_dm = 0.214;
    CP->Omega_m = CP->Omega_b + CP->Omega_dm;
    CP->Omega_lambda = CP->Omega_0 - CP->Omega_m;
    CP->Omega_rad = 0.0;
    CP->Omega_k = 0.0;
    CP->K_fit = 3.8;
    CP->F_fit = 0.01;
#endif
    CP->h100 = CP->H0 / (100 / 1e3); /*little h  = H0/(100 km/s/Mpc)*/

    fprintf(stderr, "....done\n");
}

double get_age(float z)
{
    XASSERT(PARAMS.COSMO != NULL, "Cosmology struct not set\n");
    const int NWORKSPACE = 1000;
    const double RECOMBINATION_REDSHIFT = 1e3;
    const double AGE_AT_RECOMBINATION = 0.37 * 1e-3; /*in Gyr ( recombination = 0.37 Myr)*/
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(NWORKSPACE);
    gsl_function F;
    double dummy = 0.0;
    double result = 0.0, error = 0.0;

    F.function = &agefunc;
    F.params = &dummy;
    gsl_integration_qags(&F, z, RECOMBINATION_REDSHIFT, 0, 1e-7, NWORKSPACE, w, &result, &error);
    result *= 9.77813 / PARAMS.COSMO->h100;

    result += AGE_AT_RECOMBINATION;

    gsl_integration_workspace_free(w);
    return result;
}

double agefunc(double z, void *params)
{
    return 1.0 / (epeebles(z) * (1.0 + z));
}

double epeebles(float z)
{
    double ez = sqrt(PARAMS.COSMO->Omega_m * (1.0 + z) * (1.0 + z) * (1.0 + z) +
                     PARAMS.COSMO->Omega_lambda); /*assumes flat Universe with only matter and lambda*/
    return ez;
}

float get_scalefactor(float z)
{
    const float a0 = 1.0;
    return a0 / (1 + z);
}

///

float get_overdensity(float z, struct cosmology_data *CP)
{
    /*
           From the Bryan and Norman paper ApJ, 1998, 495, 80
           In their notation, Omega_0 is my Omega_m. My Omega_k
           is their Omega_R

    */

    float omega, Ez, x;
    Ez = sqrt(CP->Omega_lambda + CP->Omega_k * SQR(1.0 + z) + CP->Omega_m * CUBE(1.0 + z) +
              CP->Omega_rad * SQR(SQR(1.0 + z)));
    omega = CP->Omega_m * CUBE(1.0 + z) / SQR(Ez);
    x = omega - 1.0;
    return (18.0 * SQR(PI) + 18.0 * x - 39.0 * SQR(x)) / omega; /* As explained in Bullock 2001 paper (for that factor
                                                                   of omega in the denom.) */
}

double get_rhocrit_at_redshift(float z, struct cosmology_data *CP)
{
    double rhocrit;
    double hubble;
    hubble = get_hubble_at_redshift(z, CP);
    rhocrit = 3.0 * (hubble * hubble / 1e6) / (8.0 * PI * CP->G);
    return rhocrit;
}

double get_hubble_at_redshift(float z, struct cosmology_data *CP)
{
    double hubble;
    double oneplusz = 1.0 + z;
    double oneplusz_sqr = (1.0 + z) * (1.0 + z);
    hubble = sqrt(CP->Omega_lambda + CP->Omega_k * oneplusz_sqr + CP->Omega_m * oneplusz * oneplusz_sqr +
                  CP->Omega_rad * oneplusz_sqr * oneplusz_sqr) *
             CP->H0 * 1e3;
    return hubble;
}

double getrvir_anyl(const double mvir, const float z, struct cosmology_data *CP)
{
    double rvir;
    float overdensity = get_overdensity(z, CP);
    double rhocrit = get_rhocrit_at_redshift(z, CP);

    rvir = pow(mvir / (4. / 3. * PI * overdensity * rhocrit), 1.0 / 3.0);
    return rvir * (1.0 + z); /*does this need to be converted to co-moving ? */
}

float getconc_anyl(const double mvir, const float z)
{
    float conc = -1.0;

    /* Going to follow Maccio et al 2008 mnras 391, 1940 since that actually has
       data for WMAP5. c_200 and not c_vir Neglect redshift evolution for the time
       being.
    */

#ifdef WMAP5
    conc = pow(10.0, 0.917 - 0.104 * log10(mvir / 1e2)); /*compute pretending to be z=0*/
#endif

#ifdef WMAP3
    conc = pow(10.0, 0.830 - 0.098 * log10(mvir / 1e2));
#endif

#ifdef WMAP1
    conc = pow(10.0, 0.769 - 0.083 * log10(mvir / 1e2));
#endif

    return conc;
}

void getrvir_from_overdensity(struct group_data *group, int NBINS, const double RhoCrit, const float OverDensity)
{
    /* Get an overdensity by making an average density profile in logarithmic bins
     */

    /* I am worried about using a constant number of bins. Might be a better idea
       to get a fixed number of particles/bin and allowing non-uniform binsizes
       etc. Or at the very least, vary the binsize according to the total number
       of particles.
    */

    double *r = NULL, *rho = NULL;
    float r_minus1, r1;
    int64 *numberdensity = NULL;
    float rmax, rbinsize, rbin;
    float xcen, ycen, zcen;
    int i;
    int index;
    float maxoverdensity = 0.0;
    const int Npart = group->N;
    const double halfmass = group->Mtot * 0.5;
    float rmin = PARAMS.BOXSIZE;
    int nbins = NBINS;

    if (Npart <= 1)
    {
        fprintf(stderr, "Number of particles needs to be greated than 1 ..exiting \n\n");
        exit(EXIT_FAILURE);
    }
    else if (Npart <= 100)
        nbins = 20;

    // Setting the global massarr variable (if not already set)
    if (PARAMS.MASSARR[DM_PART_TYPE] <= 0.0 && group->Mtot > 0.0 && Npart > 0)
    {
        PARAMS.MASSARR[DM_PART_TYPE] = group->Mtot / Npart;
    }

    r = (double *)my_malloc(sizeof(*r), Npart);
    numberdensity = (int64 *)my_malloc(sizeof(*numberdensity), nbins);
    rho = (double *)my_malloc(sizeof(*rho), nbins);

    if (group->Mtot > 0 && fabs(PARAMS.MASSARR[DM_PART_TYPE] * Npart - group->Mtot) / group->Mtot > 0.01)
    {
        fprintf(stderr,
                "Particle masses may not be set correctly. Npart = %d Mtot = %e "
                "Mpart*Npart = %e Mpart = %e \n",
                Npart, group->Mtot, PARAMS.MASSARR[DM_PART_TYPE] * Npart, PARAMS.MASSARR[DM_PART_TYPE]);
        exit(EXIT_FAILURE);
    }

    xcen = group->xcen; // or use group->x[0]??
    ycen = group->ycen;
    zcen = group->zcen;
    assert(xcen >= 0 && xcen <= PARAMS.BOXSIZE && "xcen must be in [0.0, BoxSize]");
    assert(ycen >= 0 && ycen <= PARAMS.BOXSIZE && "ycen must be in [0.0, BoxSize]");
    assert(zcen >= 0 && zcen <= PARAMS.BOXSIZE && "zcen must be in [0.0, BoxSize]");

    /*
           Needs proper handling of box-wrapping -- taken from Groupfinder.

    */

    for (i = 0; i < nbins; i++)
    {
        rho[i] = 0.0;
        numberdensity[i] = 0;
    }

    /* In principle, I could update xcen's until rho[nbins-1] is the highest.
           But no guarantee for convergence given all the halos and the numerical
       noise. Note that with my convention rho[0] is the outermost particle.

           Switching the convention so that r[0] -> innermost.
     */

    rmax = 0.0;
    for (i = 1; i < Npart; i++)
    {
        const float dx = periodic(group->x[i] - xcen);
        const float dy = periodic(group->y[i] - ycen);
        const float dz = periodic(group->z[i] - zcen);

        r[i] = sqrt(dx * dx + dy * dy + dz * dz);
        if (r[i] > rmax)
            rmax = r[i];

        if (r[i] < rmin && r[i] > 0.0)
            rmin = r[i];
    }

    assert(rmin > 0 && "Min. radius must be non-zero");
    assert(rmax > rmin && "Max. radius must be greater than minimum radius");
    assert(nbins > 0 && "Number of bins must be non-zero");
    rbinsize = (log10(rmax) - log10(rmin)) / nbins;
    /*   for(i=0;i<nbins;i++) */
    /* 	fprintf(stderr,"rmax = %f nbins = %d rbin[%d] = %f
     * \n",rmax,nbins,i,rmin*pow(10.0,i*rbinsize)); */

    /* particle 0 is in bin 0 */
    numberdensity[0] = 1;
    //   rho[0] += PARAMS.MASSARR[group->type[0]];
    rho[0] += PARAMS.MASSARR[DM_PART_TYPE];

    for (i = 1; i < Npart; i++)
    {
        index = 0;
        if (r[i] >= rmin && r[i] <= rmax)
            index = (int)floor((log10(r[i]) - log10(rmin)) / rbinsize);

        if (index >= nbins) /* Make sure there are seg. faults */
            index = nbins - 1;

        /* 	  fprintf(stderr,"in get rvir from overdensity. nbins = %4d index = %4d
         * rmax = %10.3f r[%d] = %10.3f rbinsize  = %10.4f
         * \n",nbins,index,rmax,i,r[i],rbinsize); */
        numberdensity[index] += 1; /* Should be used as error estimation for fitting */
                                   //   rho[index] +=  PARAMS.MASSARR[group->type[i]];
        rho[index] += PARAMS.MASSARR[DM_PART_TYPE];
    }

    /* Now get the cumulative mass, i.e. M( < r) */
    for (i = 1; i < nbins; i++)
        rho[i] += rho[i - 1];

    /*rho currently contains mass enclosed within i'th bin -> get half mass radius
     */
    i = 0;
    while (i < nbins && rho[i] < halfmass)
        i++;

    if (i < nbins && i > 1 && rho[i] > halfmass && rho[i - 1] < halfmass)
    {
        r_minus1 = pow(10.0, (i - 1) * rbinsize + log10(rmin));
        r1 = pow(10.0, i * rbinsize + log10(rmin));
        group->Rhalf = r1 * (halfmass - rho[i - 1]) +
                       r_minus1 * (rho[i] - halfmass); /* linear interpolation to half-mass radius */
        group->Rhalf /= (rho[i] - rho[i - 1]);
    }
    else
        group->Rhalf = (rmax + rmin) * 0.5;

    for (i = 0; i < nbins; i++)
    {
        rbin = pow(10.0, i * rbinsize + log10(rmin));
        rho[i] /= (4. / 3. * PI * pow(rbin, 3.0)); /* computes average density and NOT actual density */
        rho[i] /= RhoCrit;                         /* Convert rho to overdensity  -- have to set Cosmology first*/
    }

    /*
          Now lets first search through the entire density profile and see if we
       are actually crossing the overdensity limit. Its subject to numerical noise
       as well.
     */

    /*
           Should also check that drho/dr < 0 for all r. However, that is just the
           best case scenario -- might not be valid for unrelaxed halos.
    */

    i = 0;
    maxoverdensity = rho[0];
    while (i < (nbins - 1) && rho[i] > OverDensity)
    {
        if (maxoverdensity < rho[i])
            maxoverdensity = rho[i];
        i++;
    }

    group->MaxOverDensity = maxoverdensity;
    group->OverDensityThresh = OverDensity;
    group->Conc = getconc_anyl(group->Mtot, REDSHIFT[group->snapshot]);

    if (i > 0 && i < nbins)
    {
        group->Rvir = pow(10.0, i * rbinsize + log10(rmin));
    }
    else
    {
        if (i == nbins)
        {
            group->Rvir = rmax;
        }
        else
        {
            group->Rvir = 0.0;
        }
    }

    // this can not happen really
    if (group->Rvir > 0.0 && group->Rhalf > group->Rvir)
        group->Rhalf = 0.5 * group->Rvir;

    /*   my_free((void **) &r); */
    /*   my_free((void **) &rho); */
    /*   my_free((void **) &numberdensity); */

    free(r);
    free(rho);
    free(numberdensity);
}
