#include "impulse.h"
#include "defs.h"
#include "read_param.h"
//
//  a currently sparse set of functions designed to figure out 'physically
//  important' flyby encounters.
//

double get_impulse(struct node_data *g, struct node_data *f, float rsep, float vsep)
{
    double Impulse = 0;
    double gc;
    double Rvir = f->Rvir_anyl;
    double Mtot2 = g->Mtot;
    double Mtot1 = f->Mtot;
    float Conc = f->Conc;
    double G = PARAMS.COSMO->G;

    if (Conc < 1.0)
        Conc = 1.0;

    gc = 1.0 / (log(1.0 + Conc) - (Conc / (1.0 + Conc)));

    //  Impulse=(pow(f->Rvir,3.0)*pow(g->Mtot,2.0));
    ///
    //(3.0*pow(rsep,4.0)*pow(vsep,2.0)*Conc*pow(gc,2.0)*f->Mtot);

    Impulse =
        (G * pow(Rvir, 3.0) * pow(Mtot2, 2.0)) / (3.0 * pow(rsep, 4.0) * pow(vsep, 2.0) * Conc * pow(gc, 2.0) * Mtot1);
    //  fprintf(stderr,"in Impulse %f %f %f %f %f %f %f\n",gc, Mtot1, Mtot2, Rvir,
    //  Conc, rsep, vsep);
    return Impulse;
}

double get_external_delpot_prim(struct node_data *secondary, struct node_data *primary, float rsep, float vsep)
{
    double DelPotPrim = 0.0;
    double Rvir1 = primary->Rvir_anyl;
    double Mtot1 = primary->Mtot;
    double Mtot2 = secondary->Mtot;
    double Omega;
    double G = PARAMS.COSMO->G;

    if (secondary->haloid == primary->haloid && secondary->snapshot == primary->snapshot)
        DelPotPrim = 0.0;
    else
    {

        Omega = (vsep / rsep) * sqrt(pow(Rvir1, 3.0) / (G * Mtot1));

        if (Omega <= 0.5)
        {
            DelPotPrim = 1.3e-3 * (pow(Mtot2, 2.0) / pow(Mtot1, 2.0)) * (pow(Rvir1, 6.0) / pow(rsep, 6.0));
        }
        else if ((Omega > 0.5) && (Omega <= 1.0))
        {
            DelPotPrim = 7e-4 * (pow(Mtot2, 2.0) / pow(Mtot1, 2.0)) * (pow(Rvir1, 6.0) / pow(rsep, 6.0));
        }
        else if ((Omega > 1.0) && (Omega <= 2.5))
        {
            DelPotPrim = 4.7e-4 * (pow(Mtot2, 2.0) / pow(Mtot1, 2.0)) * (pow(Rvir1, 6.0) / pow(rsep, 6.0));
        }
        else if ((Omega > 2.5) && (Omega <= 5.0))
        {
            DelPotPrim = 2.8e-4 * (pow(Mtot2, 2.0) / pow(Mtot1, 2.0)) * (pow(Rvir1, 6.0) / pow(rsep, 6.0));
        }
        else if (Omega > 5.0)
        {
            DelPotPrim = 0.0;
        }
    }

    return DelPotPrim;
}

double get_external_delpot_sec(struct node_data *secondary, struct node_data *primary, float rsep, float vsep)
{
    double DelPotSec = 0.0;
    double Rvir2 = secondary->Rvir_anyl;
    double Mtot2 = secondary->Mtot;
    double Mtot1 = primary->Mtot;
    double Omega;
    double G = PARAMS.COSMO->G;
    if (secondary->haloid == primary->haloid && secondary->snapshot == primary->snapshot)
        DelPotSec = 0.0;
    else
    {
        Omega = (vsep / rsep) * sqrt(pow(Rvir2, 3.0) / (G * Mtot2));

        if (Omega <= 0.5)
        {
            DelPotSec = 1.3e-3 * (pow(Mtot1, 2.0) / pow(Mtot2, 2.0)) * (pow(Rvir2, 6.0) / pow(rsep, 6.0));
        }
        else if ((Omega > 0.5) && (Omega <= 1.0))
        {
            DelPotSec = 7e-4 * (pow(Mtot1, 2.0) / pow(Mtot2, 2.0)) * (pow(Rvir2, 6.0) / pow(rsep, 6.0));
        }
        else if ((Omega > 1.0) && (Omega <= 2.5))
        {
            DelPotSec = 4.7e-4 * (pow(Mtot1, 2.0) / pow(Mtot2, 2.0)) * (pow(Rvir2, 6.0) / pow(rsep, 6.0));
        }
        else if ((Omega > 2.5) && (Omega <= 5.0))
        {
            DelPotSec = 2.8e-4 * (pow(Mtot1, 2.0) / pow(Mtot2, 2.0)) * (pow(Rvir2, 6.0) / pow(rsep, 6.0));
        }
        else if (Omega > 5.0)
        {
            DelPotSec = 0.0;
        }
    }
    return DelPotSec;
}

double get_internal_delpot_prim(struct node_data *secondary, struct node_data *primary, float rsep, float vsep)
{
    /* When the perturber is inside. Table 2 of Vesperini & Weinberg 1999 */

    double K, DelPotPrim, alpha;
    /*   const float half_mass_factor = 0.5; */
    /*   double R_half = primary->Rvir_anyl*half_mass_factor; /\*very approx.
     * Should use eqn. 28 in Lokas 2001. However that  */
    /* 														  needs a proper c
     * -- which has not yet been implemented*\/ */
    double R_half = primary->Rhalf;
    double ratio = rsep / R_half;
    const float big_perturbation = 999999999.999999;

    float Mtot1 = primary->Mtot;
    float Mtot2 = secondary->Mtot;
    double temp = Mtot2 / (0.1 * Mtot1);
    /*   const int maxUlps = 5; */
    /*   if (float_almost_equal(rsep,0.0,maxUlps) == 0) /\* i.e, rsep != 0.0 *\/
     */

    if (secondary->haloid == primary->haloid && secondary->snapshot == primary->snapshot)
        DelPotPrim = 0.0;
    else
    {
        if (fabs(rsep) > 0.0)
        {
            if (rsep >= primary->Rvir_anyl || ratio > 2.0) /* quite far flung -> may be ~Rvir -> use external
                                                              prescription */
            {
                DelPotPrim = get_external_delpot_prim(secondary, primary, rsep, vsep);
            }
            else
            {
                if (ratio > 0.0 && ratio <= 1.0)
                {
                    K = 0.5;
                    alpha = 1.9;
                }
                else
                {
                    K = 1.0;
                    alpha = 1.96;
                }
                DelPotPrim = temp * temp * K / pow(vsep / 200.0, alpha);
            }
        }
        else
            DelPotPrim = (double)big_perturbation;
    }
    return DelPotPrim;
}

double get_internal_delpot_sec(struct node_data *secondary, struct node_data *primary, float rsep, float vsep)
{
    /* When the perturber is inside. Table 2 of Vesperini & Weinberg 1999 */

    double K, DelPotSec, alpha;
    /*   const float half_mass_factor = 0.5; */
    /*   double R_half = secondary->Rvir_anyl*half_mass_factor; /\*very approx.
     * Should use eqn. 28 in Lokas 2001. However that  */
    /* 														  needs a proper c
     * -- which has not yet been implemented*\/ */

    double R_half = secondary->Rhalf;
    double ratio = rsep / R_half;
    const float big_perturbation = 999999999.999999;

    float Mtot1 = primary->Mtot;
    float Mtot2 = secondary->Mtot;
    double temp = Mtot1 / (0.1 * Mtot2);
    /*   const int maxUlps = 5; */
    /*   if (float_almost_equal(rsep,0.0,maxUlps) == 0) /\* i.e, rsep != 0.0 *\/
     */

    if (secondary->haloid == primary->haloid && secondary->snapshot == primary->snapshot)
        DelPotSec = 0.0;
    else
    {
        if (fabs(rsep) > 0.0)
        {
            if (rsep >= secondary->Rvir_anyl || ratio > 2.0) /* quit far flung -> may be ~Rvir -> use external
                                                                prescription */
            {
                DelPotSec = get_external_delpot_sec(secondary, primary, rsep, vsep);
            }
            else
            {

                if (ratio > 0.0 && ratio <= 1.0)
                {
                    K = 0.5;
                    alpha = 1.9;
                }
                else
                {
                    K = 1.0;
                    alpha = 1.96;
                }
                DelPotSec = temp * temp * K / pow(vsep / 200.0, alpha);
            }
        }
        else
            DelPotSec = (double)big_perturbation; /*rsep should only be 0 only when the primary and
                                                     secondary are the same, i.e., dissolving in case
                                                     of a merger or a detachment in a flyby */
    }

    return DelPotSec;
}
