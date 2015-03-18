#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define MAXLEN    1000

const char modelnames[][MAXLEN] = {"Conroy & Wechsler (2009)","Moster et al (2010)","Behroozi et al (2010)","Behroozi et al (2013)"};
const int nmodels = sizeof(modelnames)/MAXLEN;
const float MinZ_for_Models[] = {2.0,2.0,2.0,8.0};

double behroozi_fx(const double x, const double alpha, const double delta, const double gamma)
{
  if(x == 0.0) {
    return -log10(2.0) + delta* pow(log10(2.0), gamma)/(1.0 + exp(1.0));
  } else {
    double first_term = -log10(pow(10.0,alpha*x) + 1.0);
    double second_term_numerator = pow(log10(1.0 + exp(x)), gamma);
    double second_term_denom   = 1.0 + exp(pow(10.0,-x));
    double second_term = delta * second_term_numerator/second_term_denom;
    return first_term + second_term;
  }
}  


double mhalo_to_mstar(const double Mvir, const float z, int model)
{
  double alpha,beta,gamma,m,M1,M2;//
  double mstar=0.0;

  if(model < nmodels) {
    if(z <= MinZ_for_Models[model]) {
      switch(model) {
      case 0: //Conroy & Wechsler (2009)
	{
	  M1 = pow(10.0,  0.056*z*z + 0.068*z + 9.5);
	  M2 = pow(10.0,  0.320*z*z + 0.018*z + 11.2);
	  alpha = 0.021*pow(z,4.86) + 3.39;
	  beta  = 0.085*z + 0.36;
	  
	  mstar = M1*pow(Mvir,alpha)*pow(M2,-beta)*pow(0.5*(M2+Mvir),  beta-alpha);
	  break;
	}
      case 1://Moster (2010)
	{
	  m  = 0.0282*pow(1.0+z, -0.72);
	  M1 = pow(10.0, 11.884*pow(1.0+z, 0.019));
	  beta  = 0.17*z + 1.06;
	  gamma = 0.556*pow(1.0+z, -0.26);
	
	  mstar = 2.0*Mvir*m/( pow(Mvir/M1, -beta) + pow(Mvir/M1,gamma) );
	  break;
	}
      case 2://Behroozi (2010)
	{
	  M1 = pow(10.0,0.03500*z*z - 0.19200*z + 10.199);
	  M2 = pow(10.0,0.00509*z*z + 0.00299*z + 11.824);
	  alpha =       -0.20760*z*z + 0.75200*z + 2.423;
	  beta  =        0.12000*z*z - 0.09940*z + 0.206;  
	  
	  mstar = M1*pow(Mvir,alpha)*pow(M2,-beta)*pow(0.5*(M2+Mvir),  beta-alpha);
	  break;
	}
      case 3://Behroozi (2013)
	{
	  double scale_factor = 1.0/(1.0 + (double) z);
	  double nu = exp(-4.0*scale_factor*scale_factor);
	  double log10_epsilon = -1.777 + (-0.006*(scale_factor-1.0) + (-0.0)*z ) * nu +
	    (-0.119 * (scale_factor-1.0));
	  
	  double log10M1 = 11.514 + (-1.793*(scale_factor-1.0) + (-0.251)*z ) * nu ;
	  alpha = -1.412 + (0.731*(scale_factor-1.0))*nu;
	  double delta = 3.508 + (2.608 *(scale_factor-1.0) + (-0.043)*z )*nu;
	  gamma = 0.316 + (1.319 *(scale_factor-1.0) + (0.279 )*z )*nu; 
	  
	  double first_term  = log10_epsilon + log10M1;
	  double log10Mh_over_M1 = log10(Mvir)-log10M1;
	  double second_term = behroozi_fx(log10Mh_over_M1,alpha,delta,gamma);
	  double third_term  = behroozi_fx(0.0, alpha,delta,gamma);
	  
	  mstar = pow(10.0,first_term + second_term - third_term);
	  break;
	}
      default:
	{
	  fprintf(stderr,"Mvir-Mstar model = %d not implemented\n The options for assigning stellar mass as a function of Mvir are :\n",model);
	  for(int i=0;i<nmodels;i++)
	    fprintf(stderr,"%s  [%d]\n",modelnames[i],i);
	  exit(EXIT_FAILURE);
	}
      }

      if(mstar <= 0.0 || mstar >= Mvir) {
	fprintf(stderr,"mstar has an unphysical value (with model = %s [option %d]). Mvir = %lf at z = %f with mstar = %lf\n",modelnames[model],model,Mvir,z,mstar);
	fprintf(stderr,"exiting..\n");
	exit(EXIT_FAILURE);
      }
    } else {
      mstar = 0.0;
    }
  } else {
    fprintf(stderr,"Mvir-Mstar model = %d not implemented\n The options for assigning stellar mass as a function of Mvir are :\n",model);
    for(int i=0;i<nmodels;i++)
      fprintf(stderr,"%s  [%d]\n",modelnames[i],i);
    exit(EXIT_FAILURE);
  }

  return mstar;
}

int main(int argc, char **argv)
{
  //Taken from data compiled by Stewart, K arxiv:1109.3207v1 Table 1
  double mstar=0.0;
  float z=0.0;
  int output_table = 0;
  if(argc < 4) {
    fprintf(stderr,"ERROR: Supply <Mvir, Msun/h>  <model> <Table output (> 0 if a table of Mhalo, Mstellar is desired> [z]\n");
    for(int imodel=0;imodel<nmodels;imodel++) {
      fprintf(stderr,"Model [%d]  -> %s \n",imodel, modelnames[imodel]);
    }
    return EXIT_FAILURE;
  }
  
  double Mvir = atof(argv[1]);
  int model = atoi(argv[2]);
  output_table = atoi(argv[3]);
  assert(model >=0 && model < nmodels && "model is within implemented list");
  assert(Mvir >= 1e6 && Mvir <= 1e16 && "Mass is in Msun/h units");
  if(argc > 4) {
    z = atof(argv[4]);
  } 

  if(output_table > 0) {
    const double minmass=1e9,maxmass=1e15;
    const int nbins=100000;
    const double binsize = (maxmass-minmass)/(nbins-1.0);
    Mvir = minmass;
    for(int i=0;i<nbins;i++) {
      mstar = mhalo_to_mstar(Mvir,z,model);
      fprintf(stdout," %e  %e %f\n",Mvir,mstar,z);
      Mvir += binsize;
    }
  } else {
    mstar = mhalo_to_mstar(Mvir,z,model);
    fprintf(stdout,"Mstar = %e Mvir = %e z = %f\n",mstar,Mvir,z);
  }

  return EXIT_SUCCESS;
}

