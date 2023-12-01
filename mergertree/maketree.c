#include "defs.h"
#include "maketree.h"
#include "read_param.h"
#include "utils.h"
#include "progressbar.h"
#include "impulse.h"
#include "loadfillprogenitors.h"

void increment_mergers(struct node_data *BaseNode,int64 igroup,struct node_data *FofHalo);
void increment_fof_nmergers(struct node_data *FofHalo,int incr);
void check_if_halo_disappears(struct node_data *node,const int numsnaps,const char *outpath,const int tag);
void print_header_disruption(FILE *fp);
double behroozi_fx(const double x, const double alpha, const double delta, const double gamma);

float get_separation_between_centres(struct node_data *g,struct node_data *f)
{
	if (g->haloid == f->haloid)
		return 0.0;

	float dx = periodic(g->xcen - f->xcen);
	float dy = periodic(g->ycen - f->ycen);
	float dz = periodic(g->zcen - f->zcen);

	float r2 = dx*dx + dy*dy + dz*dz ;

	return sqrtf(r2);
}


int halo_is_always_fof(struct node_data *fof)
{
	struct node_data *thisnode = fof;
	while(thisnode != NULL)
	{
		if(thisnode->isFof==0)
			return 0;

		thisnode=thisnode->BigChild;
	}

	return 1;
}



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


double assign_stellar_mass_from_mvir(struct node_data * const thisnode,int model)
{
  //Taken from data compiled by Stewart, K arxiv:1109.3207v1 Table 1

  double alpha,beta,gamma,m,M1,M2;//
  double mstar=0.0;
  float z = thisnode->z;
  double Mvir = thisnode->InfallMass*ActualMassUnits;

  const char modelnames[][MAXLEN] = {"Conroy & Wechsler (2009)","Moster et al (2010)","Behroozi et al (2010)","Behroozi et al (2013)"};
  const int nmodels = sizeof(modelnames)/MAXLEN;
  const float MinZ_for_Models[] = {2.0,2.0,2.0,8.0};

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
  return mstar/ActualMassUnits;
}



void print_header_disruption(FILE *fp)
{
  if(fp != NULL)
	{
	  fprintf(fp,"#########################################################################\n");
	  fprintf(fp,"#   Snapshot       Groupnum           Haloid            Mtot         Tag \n");
	  fprintf(fp,"#      i              l                  l                d           i  \n");
	  fprintf(fp,"#########################################################################\n");
	}

}


void check_if_halo_disappears(struct node_data *node,const int numsnaps,const char *outpath,const int tag)
{
  struct node_data *thisnode = node;
  int count = 0;
  char outfname[MAXLEN];
  FILE *fp=NULL;
  my_snprintf(outfname,MAXLEN,"%s/%s",outpath,"halos_that_disrupt.txt");
  fp = my_fopen_carefully(outfname,&print_header_disruption);

  while (thisnode != NULL && count < numsnaps) {
    if(thisnode->Parent != NULL) {
      if(thisnode->Parent->haloid == thisnode->haloid) {
	thisnode = thisnode->Parent;
      } else {
	thisnode = NULL; /* so that the following "if"  condition is not satisfied*/
	break;
      }
    } else {
      break;
    }
    count = thisnode->snapshot - node->snapshot;
  }

  if(thisnode != NULL && count < numsnaps && thisnode->snapshot < PARAMS.MAX_SNAPSHOT_NUM) {
    /* so after the encounter, the halo actaully vanished from the simulation.
       let's output data for this halo.
    */

    thisnode = node;
    while(thisnode != NULL && thisnode->haloid == node->haloid) {
      fprintf(fp,"%8d    %14"STR_FMT"    %14"STR_FMT"   %14.4g   %8d\n",thisnode->snapshot,thisnode->nodeloc,
	      thisnode->haloid,thisnode->Mtot,tag);
      thisnode = thisnode->Parent;
    }
  }

  fclose(fp);

}


void assign_parent(struct node_data *halo,int64 haloid, struct node_data *parenthalo,int64 parentid)
{
  struct node_data *tmp_loc=NULL;

  parenthalo[parentid].Nchild++;
  if(parenthalo[parentid].Nchild == 1) {
    parenthalo[parentid].BigChild = &halo[haloid];
    halo[haloid].Sibling = NULL;
  } else {
    /*
       sort the children halos by mass by changing the pointers. DO NOT change the order
       in which they appear in the node structure. Make sure that the biggest appear first.
       reset all the sibling pointers (they point toward less massive halos)
    */


    tmp_loc = parenthalo[parentid].BigChild;
    while(tmp_loc->Mtot > halo[haloid].Mtot && tmp_loc->Sibling !=NULL)
      tmp_loc = tmp_loc->Sibling;

    /* Is the new halo going to be the BigChild ? */
    if(tmp_loc == parenthalo[parentid].BigChild) {
      if(tmp_loc->Mtot < halo[haloid].Mtot) {
	halo[haloid].Sibling = tmp_loc;
	parenthalo[parentid].BigChild = &halo[haloid];
      } else {
	halo[haloid].Sibling = tmp_loc->Sibling;
	tmp_loc->Sibling = &halo[haloid];
      }
    } else {
      if (tmp_loc->Sibling == NULL) {
	/* The new halo is at the end*/
	tmp_loc->Sibling = &halo[haloid];
	halo[haloid].Sibling = NULL;
      } else {
	/* The new halo is somewhere in the middle*/
	halo[haloid].Sibling = tmp_loc->Sibling;
	tmp_loc->Sibling = &halo[haloid];
      }
    }
  }
  halo[haloid].Parent = &parenthalo[parentid];
}



void maketree(struct parent_data* allparents[],int64 *Ngroups,struct node_data * tree[])
{
  struct parent_data *p=NULL;
  struct node_data   *ParentNode=NULL,*BaseNode=NULL;
  int parentsnapshot;
  int64 parentid;
  int interrupted;
  fprintf(stderr,"\n\nAssigning parents.. \n");
  /* PRINTSTEP = (int)floor(0.1*(PARAMS.MAX_SNAPSHOT_NUM-PARAMS.MIN_SNAPSHOT_NUM)); */
  /* SMALLPRINTSTEP = ceil(0.01*(PARAMS.MAX_SNAPSHOT_NUM-PARAMS.MIN_SNAPSHOT_NUM)) > 1 ? ceil(0.01*(PARAMS.MAX_SNAPSHOT_NUM-PARAMS.MIN_SNAPSHOT_NUM)):1; */
  init_my_progressbar(PARAMS.MAX_SNAPSHOT_NUM-PARAMS.MIN_SNAPSHOT_NUM+1,&interrupted);
  for(int isnapshot=PARAMS.MIN_SNAPSHOT_NUM;isnapshot<PARAMS.MAX_SNAPSHOT_NUM;isnapshot++) {
    my_progressbar(isnapshot-PARAMS.MIN_SNAPSHOT_NUM,&interrupted);
    p = allparents[isnapshot];
    BaseNode = tree[isnapshot];

    for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup++) {
      parentsnapshot = p[igroup].parentsnapshot;
      parentid       = p[igroup].parentid;
      if(parentsnapshot >= PARAMS.MIN_SNAPSHOT_NUM && parentsnapshot <= PARAMS.MAX_SNAPSHOT_NUM && parentid < Ngroups[parentsnapshot]) {
	ParentNode = tree[parentsnapshot];
	assign_parent(BaseNode,igroup,ParentNode,parentid);
      }
    }
  }

  /* fprintf(stderr,"..done\n"); */
  finish_myprogressbar(&interrupted);
}

void assign_haloid(struct node_data *tree[],int64 *Ngroups)
{
  int64 haloid = 0;
  struct node_data *BaseNode=NULL;
  struct node_data *thisnode=NULL;
  struct node_data *tmp_node=NULL;
  float formationz=-1,destructionz=-1.0;
  double infallmass=0.0;
  short infallsnap = -1;

  const int model = 3;//0 - Conroy & Wechsler (2009), 1 - Moster (2010), 2 - Behroozi (2010), 3 - Behroozi (2013)

  fprintf(stderr,"\n\n Assigning halo ids to all halos \n\n");
  haloid = 0;

  /* reset all haloid's (in case the tree has been modified) */
  for(int isnapshot=PARAMS.MIN_SNAPSHOT_NUM;isnapshot<=PARAMS.MAX_SNAPSHOT_NUM;isnapshot++) {
    BaseNode = tree[isnapshot];
    if(Ngroups[isnapshot] > 0) {
      for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup++) {
	thisnode = &BaseNode[igroup];
	thisnode->haloid = -1;

	/* Make sure that thisnode->Nchild is correct */
	if(thisnode->BigChild != NULL) {
	  thisnode->Nchild = 1;
	  tmp_node = thisnode->BigChild->Sibling;
	  while(tmp_node != NULL) {
	    thisnode->Nchild++;
	    tmp_node = tmp_node->Sibling;
	  }
	}
      }
    }
  }

  for(int isnapshot=PARAMS.MAX_SNAPSHOT_NUM;isnapshot>=PARAMS.MIN_SNAPSHOT_NUM;isnapshot--) {
    BaseNode = tree[isnapshot];
    if(Ngroups[isnapshot] > 0) {
      for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup++) {
	thisnode = &BaseNode[igroup];
	if (thisnode->haloid < 0) {
	  destructionz = thisnode->z;
	  thisnode->haloid = haloid;

	  while(thisnode->BigChild != NULL) {
	    thisnode = thisnode->BigChild;
	    if (thisnode->haloid >= 0) {
	      fprintf(stderr,"\n This should not have happened.. found a haloid while assigning haloids -- exiting \n");
	      fprintf(stderr,"snapshot = %d this haloid = %"STR_FMT"  this nodeloc = %"STR_FMT" this parent haloid  = %"STR_FMT"  this parent nodeloc = %"STR_FMT" snapshot = %hd\n",
		      thisnode->snapshot,thisnode->haloid,thisnode->nodeloc,thisnode->Parent->haloid,thisnode->Parent->nodeloc,thisnode->Parent->snapshot);
	      fprintf(stderr,"thisnode->parent->nchild = %"STR_FMT"   thisnode->parent->bigchild->haloid = %"STR_FMT" at snapshot = %d \n",
		      thisnode->Parent->Nchild,thisnode->Parent->BigChild->haloid,thisnode->Parent->BigChild->snapshot);

	      exit(EXIT_FAILURE);
	    } else {
	      thisnode->haloid = haloid;
	      thisnode->DestructionRedshift = destructionz;
	      formationz = thisnode->z;
	    }
	  }

	  infallmass = thisnode->Mtot;
	  infallsnap    = thisnode->snapshot;
	  while(thisnode != NULL) {
	    if(thisnode->isFof == 1) {
	      infallmass = thisnode->Mtot;//BoundFofMtot could also be used in principle (however, that would double count the total stellar mass)
	      infallsnap    = thisnode->snapshot;
	    }

	    if(thisnode->haloid == haloid) {
	      thisnode->FormationRedshift = formationz;
	      thisnode->InfallMass = infallmass;
	      thisnode->InfallSnapshot = infallsnap;
	      thisnode->Mstar = assign_stellar_mass_from_mvir(thisnode,model);

	      //ensure that the stellar mass does not reduce.
	      if(thisnode->BigChild != NULL && thisnode->BigChild->Mstar > thisnode->Mstar)
		thisnode->Mstar = thisnode->BigChild->Mstar;

	    } else {
	      break;
	    }
	    thisnode = thisnode->Parent;
	  }
	  haloid++;
	}
      }
    }
  }
  MaxHaloId = haloid-1;
}

void increment_fof_nmergers(struct node_data *FofHalo,int incr)
{
  int64 haloid = FofHalo->haloid;
  struct node_data *thisnode;

  FofHalo->Nmergers += incr;
  FofHalo->TotNmergers += incr;

  thisnode=FofHalo->Parent;
  while(thisnode!=NULL && thisnode->haloid == haloid)
	{
	  thisnode->TotNmergers +=incr;
	  thisnode = thisnode->Parent;
	}
}


void find_mergers(struct node_data *tree[],int64 *Ngroups)
{
  int64 haloid = 0;
  float destructionz=-1.0,lastmergerz;
  float destructionradius=0.0;

  double destruction_subtmot,destruction_parentmtot,destruction_parentmstar;
  int destruction_snapshot;

  struct node_data *BaseNode=NULL;
  struct node_data *thisnode=NULL, *savenode=NULL,*savethisnode=NULL,*tmpnode=NULL;
  /*struct node_data ,*savethisnode,*fofnode;; */
  /*   struct node_data *thisfofhalo,*bigchild, *sibling; */
  /*   int tag_fofhaloid,tag_haloid; */

  int64 totnmergers,totndissolve,totndisrupt;
  FILE *fp=NULL,*fp2=NULL;
  char outfname[MAXLEN];
  int64 count_mergers;
  int tag;
  float rsep,vsep;

  double Impulse=0.0;
  double DelPotPrim=0.0, DelPotSec=0.0,MaxDelPotPrim=0.0,MaxDelPotSec=0.0,tmpmaxdelpotsec=0.0,tmpmaxdelpotprim=0.0;

  /* variables for getting min. rsep -- required for max. delpots*/
  double minrsep=0.0,tmprsep=0.0,tmpvsep=0.0,rvir_at_minrsep=0.0,fof_rvir_at_minrsep=0.0;
  int minrsep_snap=-1;

  /* for halos that disrupt after an encounter */
  int NUM_SNAPSHOTS_TO_CHECK = 10;

  my_snprintf(outfname,MAXLEN,"%s%s%s","rm -f ",PARAMS.OUTPUT_DIR,"/halos_that_disrupt*.txt");
  fprintf(stderr,"executing system command '%s' \n",outfname);
  system(outfname); /* will get the proper use at the next line */

  my_snprintf(outfname,MAXLEN,"%s/%s",PARAMS.OUTPUT_DIR,"fof_to_subhalo_mergers.txt");
  fprintf(stderr,"About to  open file %s\n",outfname);
  fp = my_fopen(outfname,"w");
  fprintf(fp,"############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################\n");
  fprintf(fp,"#   Snapshot       ParentGrpnum      ParentHaloid    ParentMtot       SubGroupNum     Subhaloid       SubMtot      rsep         vsep        Impulse      DelPotPrim     DelPotSec       Tag          DestructionZ   MinRsep      MaxDelPotPrim     MaxDelPotSec        NormalisedDestructionRadius           ParentRvir        MinRsepSnap       SubRviratMinRsep            ParentRviratMinRsep    DestructionSnap    SubmtotDestructionSnap   ParentmtotDestructionSnap SubStellarMassDestructionSnap   ParentStellarMassDestructionsnap \n");
  fprintf(fp,"#     i               l                   l              d                l               l              d          f             f            d            d              d              i             f             f               d                d                       d                                 f                i                      f                            f                   i                       d                         d                        d                                d                    \n");
  fprintf(fp,"###########################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################\n");


  my_snprintf(outfname,MAXLEN,"%s/%s",PARAMS.OUTPUT_DIR,"initial_positions_for_flyby_orbits.txt");
  fprintf(stderr,"About to  open file %s\n",outfname);
  fp2 = my_fopen(outfname,"w");
  fprintf(fp2,"# Mass in units of 1e10 Msun/h, distance (physical) kpc/h,  velocity (peculiar)  km/s. Primary is centered at (0,0,0) \n");
  fprintf(fp2,"###############################################################################################################################################\n");
  fprintf(fp2,"# redshift   BoundPrimMass SubInfallMass   xcen        ycen         zcen          vxcen        vycen        vzcen      PrimHaloId    SubHaloId \n");
  fprintf(fp2,"#   f               d            d           f           f           f              f            f            f            l            l      \n");
  fprintf(fp2,"###############################################################################################################################################\n");


  /*
	Tag:  0 -- Fof halo falls in. becomes its own subhalo and then disrupts at some later time inside the same container Fof halo
	Tag:  1 -- The subhalo disappears because the FOF container has now merged with some other FOF halo.
	Tag:  2 -- the subhalo remains a subhalo inside that same FOF container until z=0
	Tag:  3 -- Fof halo falls in, but later becomes a fof halo on its own again. Possible flyby. Destructionz then represents
	the z at which *this* subhalo becomes a FOF halo again. (Multiple flybys by the same object should show up as another
	entry at a later time -- with the possibility of a tag 0.

  */


  /*Before anything else: merge all the Mstar masses for subhalos with Nchild > 1.
	For FOF halos, the mass used to assign Mstar increases (mostly) with `a' but the
	subhalos only use the infallmass to get their Mstar. If the subhalo merges with
	another subhalo, then that stellar mass content needs to be added in as well.
	Also, this has to proceed from high-z to low-z so that all the masses are
	propagated to z=0.
  */

  for(int isnapshot=PARAMS.MIN_SNAPSHOT_NUM;isnapshot <= PARAMS.MAX_SNAPSHOT_NUM;isnapshot++) {
    BaseNode = tree[isnapshot];
    for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup++) {
      thisnode = &BaseNode[igroup];
      if(thisnode->Nchild > 1 && thisnode->isFof == 0)	{
				savenode  = thisnode;
				thisnode = savenode->BigChild->Sibling;
				while(thisnode != NULL)	{
					savenode->Mstar += thisnode->Mstar;
					thisnode = thisnode->Sibling;
				}
      }
    }
  }


  for(int isnapshot=PARAMS.MAX_SNAPSHOT_NUM;isnapshot>=PARAMS.MIN_SNAPSHOT_NUM;isnapshot--)	{
    BaseNode = tree[isnapshot];
    if(Ngroups[isnapshot] > 0) {
      for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup++) {
				thisnode = &BaseNode[igroup];
				haloid = thisnode->haloid;
				lastmergerz  = -1.0;
				totnmergers  = 0;
				totndissolve = 0;
				totndisrupt  = 0;

				if(thisnode->Nchild > 1) {
					/* Disruptions and Dissolutions */
					if(thisnode->isFof == 1) {
						savenode = thisnode;
						thisnode = thisnode->BigChild;
						while(thisnode->Sibling != NULL) {
							thisnode = thisnode->Sibling;
							if(thisnode->FofHalo->haloid == haloid) {
								savenode->NDisruptions +=1;
								totndisrupt +=1;
								lastmergerz = savenode->z;
							} else {
								savenode->NDissolutions +=1;
								totndissolve +=1;
								lastmergerz = savenode->z;
							}

						}

						thisnode = savenode; /* reset thisnode to actual node that we are dealing with */
						assert(thisnode->haloid == haloid);
					} else {
						/*thisnode is a subhalo*/

						/*
							 Subhalo merging. Its possible that a FOF halo falls in and directly points into
							 this subhalo. In principle that could be the same as a `dissolution' event but
							 I am treating it as a `merger' event.
						*/

						savenode = thisnode;
						count_mergers=thisnode->Nchild-1;
						thisnode->Nmergers = count_mergers;
						totnmergers += count_mergers;
						lastmergerz = thisnode->z;

						thisnode = thisnode->BigChild;
						while(thisnode !=NULL && count_mergers > 0) {
							if (thisnode->FofHalo->haloid  == savenode->FofHalo->haloid) {
								count_mergers--;
							} else {
								/* 							  fprintf(stderr,"*WARNING* hmmm subhalo has more than one child coming from a different fof halo *WARNING*\n"); */
								if (savenode->haloid == thisnode->haloid)	{
									/* Figure out the tag */
									savethisnode = thisnode;
									thisnode=thisnode->Parent; /*come back down to current snapshot*/

									destructionz=-1.0;
									destructionradius=-1.0;
									destruction_subtmot=-1.0;
									destruction_parentmtot=-1.0;
									destruction_snapshot=-1;
									destruction_parentmstar=-1.0;

									/* 								  tag = get_tag(thisnode,savenode->FofHalo,&destructionz,&destructionradius); */
									//MS 7th Dec, 2011 - updated to use the new tagging system that accounts for
									//subhalos missing snapshots

									//MS 30th Sep, 2013 - no longer necessary. ensure_same_snapshot() accounts for missing snasphots
									tag = get_tag(thisnode,savenode->FofHalo,&destructionz,&destructionradius,&destruction_snapshot,&destruction_subtmot,&destruction_parentmtot,&destruction_parentmstar);
									if (tag == -1 && thisnode->Parent == NULL)
										tag = 99;

									if (tag == 3 || tag == 5) {
										thisnode->NFlybys++;
										savenode->FofHalo->NFlybys++;
										fprintf(fp2,"%e %e %e %e %e %e %e %e %e %"STR_FMT" %"STR_FMT"\n",
														savenode->z,
														savenode->FofHalo->BoundFofMtot,
														thisnode->InfallMass,
														get_scalefactor(thisnode->z)*periodic(savenode->FofHalo->xcen-thisnode->xcen),
														get_scalefactor(thisnode->z)*periodic(savenode->FofHalo->ycen-thisnode->ycen),
														get_scalefactor(thisnode->z)*periodic(savenode->FofHalo->zcen-thisnode->zcen),
														savenode->FofHalo->meanvel[0]-thisnode->meanvel[0],
														savenode->FofHalo->meanvel[1]-thisnode->meanvel[1],
														savenode->FofHalo->meanvel[2]-thisnode->meanvel[2],
														savenode->FofHalo->haloid,
														thisnode->haloid);

									}


									check_if_halo_disappears(thisnode,NUM_SNAPSHOTS_TO_CHECK,PARAMS.OUTPUT_DIR,tag);

									rsep = get_separation_between_centres(thisnode,savenode->FofHalo);
									vsep = 0.0;
									for(int i=0;i<3;i++)
										vsep += pow( (thisnode->meanvel[i] - savenode->FofHalo->meanvel[i]),2.0);

									vsep = sqrt(vsep);

									Impulse=get_impulse(thisnode,savenode->FofHalo,rsep,vsep);
									DelPotPrim=get_external_delpot_prim(thisnode,savenode->FofHalo,rsep,vsep);
									DelPotSec=get_external_delpot_sec(thisnode,savenode->FofHalo,rsep,vsep);

									/* Find the max delta phi created -> smallest rsep while subhalo is still inside this same fof */
									minrsep=rsep;
									minrsep_snap=thisnode->snapshot;
									rvir_at_minrsep=thisnode->Rvir;
									fof_rvir_at_minrsep=savenode->FofHalo->Rvir;

									tmpnode=thisnode;
									MaxDelPotSec = DelPotSec;
									MaxDelPotPrim = DelPotPrim;
									/* 								  while(tmpnode->Parent !=NULL && tmpnode->Parent->haloid==thisnode->haloid && tmpnode->FofHalo->haloid==thisnode->FofHalo->haloid) */
									while(tmpnode->Parent !=NULL && tmpnode->Parent->haloid==thisnode->haloid && tmpnode->Parent->FofHalo->haloid==thisnode->FofHalo->haloid ) {
										tmpnode=tmpnode->Parent;
										tmprsep = get_separation_between_centres(tmpnode,tmpnode->FofHalo);;
										if (tmprsep < minrsep) {
											minrsep=tmprsep;
											minrsep_snap = tmpnode->snapshot;
											rvir_at_minrsep=tmpnode->Rvir;
											fof_rvir_at_minrsep=tmpnode->FofHalo->Rvir;

										}

										if (float_almost_equal(minrsep,0.0,5) == 1  && (tag==3 || tag==5)) {
											fprintf(stderr,"WARNING (UPPER): This should not happen\n Min. rsep is 0.0 for a flyby\n");
											fprintf(stderr,"thisnode->haloid = %"STR_FMT" snapshot= %d  tmpnode->haloid=%"STR_FMT" tmpnode->snap=%d savenode->Fofid = %"STR_FMT" \n",
															thisnode->haloid,thisnode->snapshot,tmpnode->haloid,tmpnode->snapshot,savenode->FofHalo->haloid);
											fprintf(stderr,"thisnode->fofhaloid = %"STR_FMT" tmpnode->id = %"STR_FMT" tmpnode->Fofhaloid = %"STR_FMT" \n",
															thisnode->FofHalo->haloid,tmpnode->haloid,tmpnode->FofHalo->haloid);

											fprintf(stderr,"tmpnode->xcen =%f tmpnode->ycen = %f tmpnode->zcen =%f\n",
															tmpnode->xcen,tmpnode->ycen,tmpnode->zcen);
											fprintf(stderr,"tmpfof->xcen =%f tmpfof->ycen = %f tmpfof->zcen =%f\n",
															tmpnode->FofHalo->xcen,tmpnode->FofHalo->ycen,tmpnode->FofHalo->zcen);

										}


										tmpvsep = 0.0;
										for(int i=0;i<3;i++)
											tmpvsep += pow( (tmpnode->meanvel[i] - tmpnode->FofHalo->meanvel[i]),2.0);

										tmpvsep = sqrt(tmpvsep);

										tmpmaxdelpotprim = get_internal_delpot_prim(tmpnode,tmpnode->FofHalo,minrsep,tmpvsep);
										MaxDelPotPrim = tmpmaxdelpotprim > MaxDelPotPrim ? tmpmaxdelpotprim:MaxDelPotPrim;

										tmpmaxdelpotsec = get_internal_delpot_sec(tmpnode,tmpnode->FofHalo,minrsep,tmpvsep);
										MaxDelPotSec  = tmpmaxdelpotsec > MaxDelPotSec ?  tmpmaxdelpotsec:MaxDelPotSec;

									}
									//MS 08/24/2011 -- Changed Submtot to InfallMass

									/* 								  fprintf(fp," %10d     %12"STR_FMT"   %12"STR_FMT"      %12.4g    %12"STR_FMT"   %12"STR_FMT"    %12.6g   %10.4f   %10.4f  %12.4g  %12.4g  %12.4g   %6d  %18.4f   %12.4f   %12.4g %12.4g       %12.6e       %12.4g      %10d       %12.4g         %12.4g \n", */
									/* 										  savenode->FofHalo->snapshot,savenode->FofHalo->nodeloc,savenode->FofHalo->haloid,savenode->FofHalo->Mtot, */
									/* 										  thisnode->nodeloc,thisnode->haloid,thisnode->Mtot,rsep, vsep,Impulse,DelPotPrim,DelPotSec,tag,destructionz, */
									/* 										  minrsep,MaxDelPotPrim,MaxDelPotSec,destructionradius,savenode->FofHalo->Rvir,minrsep_snap, rvir_at_minrsep,fof_rvir_at_minrsep); */



									fprintf(fp," %10d     %12"STR_FMT"   %12"STR_FMT"      %12.4g    %12"STR_FMT"   %12"STR_FMT"    %12.6g   %10.4f   %10.4f  %12.4g  %12.4g  %12.4g   %6d  %18.4f   %12.4f   %12.4g %12.4g       %12.6e       %12.4g      %10d       %12.4g         %12.4g  %10d %12.4g %12.4g  %16.4g  %16.4g\n",
													savenode->FofHalo->snapshot,savenode->FofHalo->nodeloc,savenode->FofHalo->haloid,savenode->FofHalo->Mtot,
													thisnode->nodeloc,thisnode->haloid,thisnode->InfallMass,rsep, vsep,Impulse,DelPotPrim,DelPotSec,tag,destructionz,
													minrsep,MaxDelPotPrim,MaxDelPotSec,destructionradius,savenode->FofHalo->Rvir,minrsep_snap, rvir_at_minrsep,fof_rvir_at_minrsep,
													destruction_snapshot,
													destruction_subtmot,
													destruction_parentmtot,
													thisnode->Mstar,
													destruction_parentmstar);


		  thisnode = savethisnode; /* go back to previous snapshot */
		} else {
		  /*so something fell from the outside and has now disappeared inside this subhalo*/
		  tag = 0;
		  destructionz = savenode->z;
		  rsep = 0.0;
		  vsep = 0.0;
		  Impulse = 0.0;
		  DelPotPrim = 0.0;
		  DelPotSec  = 0.0;
		  minrsep=0.0;
		  minrsep_snap = savenode->snapshot;
		  rvir_at_minrsep = 0.0;
		  fof_rvir_at_minrsep = savenode->FofHalo->Rvir;
		  MaxDelPotSec=0.0;
		  MaxDelPotPrim=0.0;
		  destruction_snapshot=-1;
		  destructionradius=-1.0;
		  destruction_parentmtot=-1.0;
		  destruction_parentmstar=-1.0;

		  //MS -- 08/24/2011 - changed Submtot to Infallmass
		  /* 								  fprintf(fp," %10d     %12"STR_FMT"   %12"STR_FMT"      %12.4g    %12"STR_FMT"   %12"STR_FMT"    %12.6g   %10.4f   %10.4f  %12.4g  %12.4g  %12.4g   %6d  %18.4f  %12.4f   %12.4g %12.4g       %12.6e        %12.4g      %10d       %12.4g         %12.4g \n", */
		  /* 										  savenode->FofHalo->snapshot,savenode->FofHalo->nodeloc,savenode->FofHalo->haloid,savenode->FofHalo->Mtot, */
		  /* 										  thisnode->nodeloc,thisnode->haloid,thisnode->Mtot,rsep, vsep,Impulse,DelPotPrim,DelPotSec,tag,destructionz, */
		  /* 										  minrsep,MaxDelPotPrim,MaxDelPotSec,destructionradius,savenode->FofHalo->Rvir,minrsep_snap, rvir_at_minrsep,fof_rvir_at_minrsep); */


		  fprintf(fp," %10d     %12"STR_FMT"   %12"STR_FMT"      %12.4g    %12"STR_FMT"   %12"STR_FMT"    %12.6g   %10.4f   %10.4f  %12.4g  %12.4g  %12.4g   %6d  %18.4f  %12.4f   %12.4g %12.4g       %12.6e        %12.4g      %10d       %12.4g         %12.4g  %10d  %12.4g  %12.4g  %16.4g  %16.4g \n",
			  savenode->FofHalo->snapshot,savenode->FofHalo->nodeloc,savenode->FofHalo->haloid,savenode->FofHalo->Mtot,
			  thisnode->nodeloc,thisnode->haloid,thisnode->InfallMass,rsep, vsep,Impulse,DelPotPrim,DelPotSec,tag,destructionz,
			  minrsep,MaxDelPotPrim,MaxDelPotSec,destructionradius,savenode->FofHalo->Rvir,minrsep_snap, rvir_at_minrsep,fof_rvir_at_minrsep,
			  destruction_snapshot,
			  destruction_subtmot,
			  destruction_parentmtot,
			  thisnode->Mstar,
			  destruction_parentmstar);



		}

	      }
	      thisnode = thisnode->Sibling;
	    }

	    thisnode  = savenode;
	    assert(thisnode->haloid == haloid);
	    increment_fof_nmergers(thisnode->FofHalo,count_mergers);
	  }
	} else {
	  /*So, I have at most one child. Now if the Fofhalo thisnode and bigchild are the same then it's
	    normal growth (without mergers. However, if the FofHalo is different, then this subhalo has
	    been created by an infalling FOFHalo.
	  */

	  if(thisnode->Nchild == 1) { /*Nchild could be 0 as well at this section in the code -> hence choose Nchild==1*/
	    savenode = thisnode;//MS Dec 8, 2011 - this line was not there. However, that was
	    //causing incorrect merger tags to appear. And looking at how the tags are done
	    //I believe this line is correct/required.
	    if(thisnode->isFof==0) {
	      if(thisnode->FofHalo->haloid != thisnode->BigChild->FofHalo->haloid) {
		count_mergers = 1;
		thisnode->Nmergers += count_mergers;
		totnmergers +=count_mergers;
		lastmergerz = thisnode->z;
		increment_fof_nmergers(thisnode->FofHalo,count_mergers);

		destructionz=-1.0;
		destructionradius=-1.0;
		destruction_subtmot=-1.0;
		destruction_parentmtot=-1.0;
		destruction_snapshot=-1;
		destruction_parentmstar=-1.0;

		/* Figure out the tag */
		/* 							  tag = get_tag(thisnode,thisnode->FofHalo,&destructionz,&destructionradius); */
		//MS 7th Dec, 2011 - updated to use the new tagging system that accounts for
		//subhalos missing snapshots

		//MS 30th Sep, 2013 - no longer necessary. ensure_same_snapshot() accounts for missing snasphots
		tag = get_tag(thisnode,thisnode->FofHalo,&destructionz,&destructionradius,&destruction_snapshot,&destruction_subtmot,&destruction_parentmtot, &destruction_parentmstar);
		if (tag == -1 && thisnode->Parent == NULL)
		  tag = 99;


		rsep = get_separation_between_centres(thisnode,thisnode->FofHalo);
		vsep = 0.0;
		for(int i=0;i<3;i++)
		  vsep += pow( (thisnode->meanvel[i] - thisnode->FofHalo->meanvel[i]),2.0);

		vsep = sqrt(vsep);


		if (tag == 3 || tag == 5) {
		  thisnode->NFlybys++;
		  thisnode->FofHalo->NFlybys++;

		  fprintf(fp2,"%e %e %e %e %e %e %e %e %e %"STR_FMT"  %"STR_FMT"\n",
			  thisnode->z,
			  thisnode->FofHalo->BoundFofMtot,
			  thisnode->InfallMass,
			  get_scalefactor(thisnode->z)*periodic(thisnode->FofHalo->xcen-thisnode->xcen),
			  get_scalefactor(thisnode->z)*periodic(thisnode->FofHalo->ycen-thisnode->ycen),
			  get_scalefactor(thisnode->z)*periodic(thisnode->FofHalo->zcen-thisnode->zcen),
			  thisnode->FofHalo->meanvel[0]-thisnode->meanvel[0],
			  thisnode->FofHalo->meanvel[1]-thisnode->meanvel[1],
			  thisnode->FofHalo->meanvel[2]-thisnode->meanvel[2],
			  thisnode->FofHalo->haloid,
			  thisnode->haloid);

		}


		/* 							  if (tag == 3 || tag == 5) */
		check_if_halo_disappears(thisnode,NUM_SNAPSHOTS_TO_CHECK,PARAMS.OUTPUT_DIR,tag);


		Impulse=get_impulse(thisnode,thisnode->FofHalo,rsep,vsep);
		DelPotPrim=get_external_delpot_prim(thisnode,thisnode->FofHalo,rsep,vsep);
		DelPotSec=get_external_delpot_sec(thisnode,thisnode->FofHalo,rsep,vsep);

		minrsep=rsep;
		minrsep_snap=thisnode->snapshot;
		tmpnode=thisnode;
		MaxDelPotSec = DelPotSec;
		MaxDelPotPrim = DelPotPrim;
		while(tmpnode->Parent !=NULL && tmpnode->Parent->haloid==thisnode->haloid && tmpnode->Parent->FofHalo->haloid==thisnode->FofHalo->haloid ) {
		  tmpnode=tmpnode->Parent;
		  tmprsep = get_separation_between_centres(tmpnode,tmpnode->FofHalo);;
		  if (tmprsep < minrsep) {
		    minrsep=tmprsep;
		    minrsep_snap = tmpnode->snapshot;
		    rvir_at_minrsep=tmpnode->Rvir;
		    fof_rvir_at_minrsep=tmpnode->FofHalo->Rvir;

		  }

		  if (float_almost_equal(minrsep,0.0,5) == 1  && (tag==3 || tag==5)) {
		    fprintf(stderr,"WARNING (LOWER): This should not happen\n Min. rsep is 0.0 for a flyby\n");
		    fprintf(stderr,"thisnode->haloid = %"STR_FMT" snapshot= %d  Fofid = %"STR_FMT" \n",
			    thisnode->haloid,thisnode->snapshot,savenode->FofHalo->haloid);
		    fprintf(stderr,"thisnode->fofhaloid = %"STR_FMT" tmpnode->id = %"STR_FMT" tmpnode->Fofhaloid = %"STR_FMT" \n",
			    thisnode->FofHalo->haloid,tmpnode->haloid,tmpnode->FofHalo->haloid);
		    fprintf(stderr,"tmpnode->xcen =%f tmpnode->ycen = %f tmpnode->zcen =%f\n",
			    tmpnode->xcen,tmpnode->ycen,tmpnode->zcen);
		    fprintf(stderr,"tmpfof->xcen =%f tmpfof->ycen = %f tmpfof->zcen =%f\n",
			    tmpnode->FofHalo->xcen,tmpnode->FofHalo->ycen,tmpnode->FofHalo->zcen);

		  }



		  tmpvsep = 0.0;
		  for(int i=0;i<3;i++)
		    tmpvsep += pow( (tmpnode->meanvel[i] - tmpnode->FofHalo->meanvel[i]),2.0);

		  tmpvsep = sqrt(tmpvsep);

		  tmpmaxdelpotprim = get_internal_delpot_prim(tmpnode,tmpnode->FofHalo,minrsep,tmpvsep);
		  MaxDelPotPrim = tmpmaxdelpotprim > MaxDelPotPrim ? tmpmaxdelpotprim:MaxDelPotPrim;

		  tmpmaxdelpotsec = get_internal_delpot_sec(tmpnode,tmpnode->FofHalo,minrsep,tmpvsep);
		  MaxDelPotSec  = tmpmaxdelpotsec > MaxDelPotSec ?  tmpmaxdelpotsec:MaxDelPotSec;
		}


		//MS -- 08/24/2011 - Changed submtot to Infallmass

		/* 							  fprintf(fp," %10d     %12"STR_FMT"   %12"STR_FMT"      %12.4g    %12"STR_FMT"   %12"STR_FMT"    %12.6g   %10.4f   %10.4f  %12.4g  %12.4g  %12.4g   %6d  %18.4f  %12.4f   %12.4g %12.4g       %12.6e         %12.4g      %10d       %12.4g         %12.4g \n", */
					/* 									  savenode->FofHalo->snapshot,savenode->FofHalo->nodeloc,savenode->FofHalo->haloid,savenode->FofHalo->Mtot, */
					/* 									  thisnode->nodeloc,thisnode->haloid,thisnode->Mtot,rsep, vsep,Impulse,DelPotPrim,DelPotSec,tag,destructionz, */
					/* 									  minrsep,MaxDelPotPrim,MaxDelPotSec,destructionradius,savenode->FofHalo->Rvir,minrsep_snap, rvir_at_minrsep,fof_rvir_at_minrsep); */


		fprintf(fp," %10d     %12"STR_FMT"   %12"STR_FMT"      %12.4g    %12"STR_FMT"   %12"STR_FMT"    %12.6g   %10.4f   %10.4f  %12.4g  %12.4g  %12.4g   %6d  %18.4f  %12.4f   %12.4g %12.4g       %12.6e         %12.4g      %10d       %12.4g         %12.4g  %10d  %12.4g  %12.4g   %16.4g  %16.4g\n",
			savenode->FofHalo->snapshot,savenode->FofHalo->nodeloc,savenode->FofHalo->haloid,savenode->FofHalo->Mtot,
			thisnode->nodeloc,thisnode->haloid,thisnode->InfallMass,rsep, vsep,Impulse,DelPotPrim,DelPotSec,tag,destructionz,
			minrsep,MaxDelPotPrim,MaxDelPotSec,destructionradius,savenode->FofHalo->Rvir,minrsep_snap, rvir_at_minrsep,fof_rvir_at_minrsep,
			destruction_snapshot,
			destruction_subtmot,
			destruction_parentmtot,
			thisnode->Mstar,
			destruction_parentmstar);

	      }
	    }
	  }

	}

	if (lastmergerz >= 0.0)
	  thisnode->RedshiftofLastMerger = lastmergerz;

	thisnode->TotNDissolutions = totndissolve;
	thisnode->TotNDisruptions  = totndisrupt;

	/*explicitly showing that TotNmergers might have already been updated at this point */
	if(thisnode->TotNmergers == 0)
	  thisnode->TotNmergers      = totnmergers;
	else
	  thisnode->TotNmergers      += totnmergers;/* in case increment_fof_mergers has been called on thisnode */

      }

    }
  }
  fclose(fp2);
  fclose(fp);

}






void find_subsub_mergers(struct node_data *tree[],int64 *Ngroups)
{
  struct node_data *BaseNode=NULL;
  struct node_data *thisnode=NULL;
  FILE *fp=NULL;
  char outfname[MAXLEN];

  my_snprintf(outfname,MAXLEN,"%s/%s",PARAMS.OUTPUT_DIR,"subhalo_to_subhalo_mergers.txt");
  fp = my_fopen(outfname,"w");

  fprintf(fp,"###################################################################################################################################################\n");
  fprintf(fp,"#  Snapshot        GroupNum        HaloId       ContainerHaloId      ParentSnapshot   ParentGrpNum      ParentHaloId   HaloidofParentContainer     \n");
  fprintf(fp,"#     i               l               l               l                  i                l                l                  l                    \n");
  fprintf(fp,"###################################################################################################################################################\n");



  for(int isnapshot=PARAMS.MAX_SNAPSHOT_NUM;isnapshot>=PARAMS.MIN_SNAPSHOT_NUM;isnapshot--) {
    BaseNode = tree[isnapshot];
    if(Ngroups[isnapshot] > 0) {
      for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup++) {
	thisnode = &BaseNode[igroup];
	if(thisnode->isFof ==0 && thisnode->Parent != NULL && thisnode->Parent->haloid != thisnode->haloid && thisnode->Parent->isFof==0
	   && thisnode->ContainerHalo->haloid != thisnode->Parent->haloid && thisnode->Parent->ContainerHalo->haloid == thisnode->ContainerHalo->haloid) {
	  //subhalo-subhalo merger; the container halo takes care of the dissolving in the same host (so the
	  //cases that should trigger here are only the ones which are (sub)-subhalos [i.e., it could be any level
	  // in the hierarchy] that merge into other (sub)-subhalos within the same containerhalo
	  fprintf(fp,"%8hd   %14"STR_FMT " %14"STR_FMT"   %14"STR_FMT"          %8hd   %14"STR_FMT "   %14"STR_FMT"     %14"STR_FMT"  \n",
		  thisnode->snapshot,thisnode->nodeloc,thisnode->haloid,thisnode->ContainerHalo->haloid,
		  thisnode->Parent->snapshot,thisnode->Parent->nodeloc,thisnode->Parent->haloid,thisnode->Parent->ContainerHalo->haloid);

	}
      }
    }
  }
  fclose(fp);

}


/*
  This sequence of evaluation is important!!
  First check that the subhalo is not a FOF halo, because
  if it is, it's irrelevant whether or not it as at the same
  snapshot as it's original FOF container.

  Everything that follows has the original FOF and subhalo
  linked together -> so it's important that they be at the
  same snapshot. So the "unknown" tag check is here.

  Then check if the FOF has fallen into another FOF because
  in that case the subhalo will lose it's id (unless subs
  of subs are loaded in). Then check for disruption and
  survivors.

*/

int get_tag(struct node_data * const node1,struct node_data * const node2,float *z,float *rnorm, int *destruction_snap, double *submtot, double *parentmtot, double *parentmstar)
{
  struct node_data *thisnode = node1;
  struct node_data *fofnode  = node2;
  struct node_data *new_fofnode=NULL;
  struct node_data *parentnode=NULL;
  struct node_data *tmpnode=NULL;
  int64 tag_fofhaloid,tag_haloid;
  int tag=-1;
  double vdotv = 0.0;
  float vxcen=0.0,vycen=0.0,vzcen=0.0;

  tag_fofhaloid = fofnode->haloid;
  tag_haloid = thisnode->haloid;

  while(thisnode != NULL && fofnode != NULL) {
      if(thisnode->isFof == 1 && thisnode->FofHalo->haloid != fofnode->FofHalo->haloid) {
          /*thisnode was a FOF when it fell in*/
        if(node1->BigChild !=NULL && node1->BigChild->isFof==1) {

            /* compute the vdotv when the other halo fell in. this should be negative */
            parentnode = node1->BigChild; /* bad naming. parentnode here is actually in the past..against the convention in the code*/
            new_fofnode  = node2->BigChild; /* since node1 (subhalo of node2) has a progenitor, node2 must have one too */
            if(new_fofnode != NULL) {
                vxcen = 0.5*(parentnode->meanvel[0]+new_fofnode->meanvel[0]);
                vycen = 0.5*(parentnode->meanvel[1]+new_fofnode->meanvel[1]);
                vzcen = 0.5*(parentnode->meanvel[2]+new_fofnode->meanvel[2]);

                vdotv = (parentnode->vxcen - vxcen)* (new_fofnode->vxcen - vxcen) + (parentnode->vycen - vycen) * (new_fofnode->vycen - vycen)+
                    (parentnode->vzcen - vzcen) * (new_fofnode->vzcen - vzcen);

                if (vdotv < 0.0)
                    tag = 3;/* fly by */
                else
                    tag = 5;


            } else {
                tag = -2; /* so somehow fofnode->BigChild is NULL but thisnode->BigChild is not.
                                                            Only way this can happen is if a new Fof (fofnode) formed and thisnode->BigChild (itself a fof)
                                                            fell into it and became a subhalo (thisnode)
                            */
            }
        } else {
            tag = 4; /* split ..
                    Note added 1st Dec, 2023 (MS): Think this means that the original subhalo, after
                    being followed into  future, is now located in a different FOF halo than the FOF halo
                    it fell into originally. (basically if this subhalo became a FOF halo,
                    then it would a classed as a flyby but it is still a subhalo but just in a different FOF halo)
                    */
        }

        *z = thisnode->z;
        *destruction_snap = thisnode->snapshot;
        *submtot = thisnode->Mtot;
        *parentmtot = fofnode->Mtot;
        *parentmstar = fofnode->Mstar;
        if(thisnode != NULL && thisnode->BigChild != NULL)
            *rnorm = get_separation_between_centres(thisnode->BigChild,thisnode->BigChild->FofHalo)/thisnode->BigChild->FofHalo->Rvir; /* more like a separation radius*/
        else
            *rnorm = -1.0;

        break;
    }


    if(fofnode->FofHalo->haloid != tag_fofhaloid  && thisnode->FofHalo->haloid == fofnode->FofHalo->haloid) {
      tag = 1; /* the Fof container has now fallen into a different FOF halo-> need subs of subs to proceed */
      /* With the hierarchy levels, I have that info.
	 But the algorithm is too complex right now, some day I will figure out how to handle this */
      *z = fofnode->z;
      *rnorm = -1.0;
      break;
    }

    if (thisnode->haloid != tag_haloid && thisnode->FofHalo->haloid == tag_fofhaloid) {
      tag = 0; /*the subhalo has disappeared, thisnode is not the same haloid
		 as the original infalling subhalo */

      tmpnode=node1;
      while(tmpnode->Parent != NULL && tmpnode->Parent->haloid == tag_haloid)
	tmpnode = tmpnode->Parent;

      if (tmpnode != NULL && thisnode->FofHalo != NULL && thisnode->FofHalo->BigChild != NULL) {
	*z = tmpnode->z;
	*rnorm = get_separation_between_centres(tmpnode,thisnode->FofHalo->BigChild)/thisnode->FofHalo->BigChild->Rvir; /* last known position of the subhalo*/
      }

      break;
    }


    if((thisnode->Parent==NULL || fofnode->Parent==NULL) &&  (thisnode->snapshot == (PARAMS.MAX_SNAPSHOT_NUM))) {
      tag = 2;/* so we have reached the end of the simulation and hence at least one pointer has disappeared */
      *z = REDSHIFT[PARAMS.MAX_SNAPSHOT_NUM];
      *rnorm = get_separation_between_centres(thisnode,thisnode->FofHalo)/thisnode->FofHalo->Rvir;

      break;
    }


    thisnode = thisnode->Parent;
    fofnode  = fofnode->Parent;

    //The following part ensures that both the halos are at the snapshot -- which
    //is non-trivial in the general case.
    if(thisnode != NULL && fofnode!=NULL && thisnode->snapshot != fofnode->snapshot)
      ensure_same_snapshot_for_halos(&thisnode,&fofnode);

  }

  return tag;
}

void ensure_same_snapshot_for_halos(struct node_data **node1,struct node_data **node2)
{

  struct node_data *thisnode;
  struct node_data *fofnode;
  int64 thisnode_haloid,fofnode_haloid;

  if( (*node1) != NULL && (*node2) !=NULL) {
    if((*node1)->snapshot < (*node2)->snapshot)	{
      fofnode  = *node1;
      thisnode = *node2;
    } else {
      fofnode  = *node2;
      thisnode = *node1;
    }

    thisnode_haloid = thisnode->haloid;
    fofnode_haloid = fofnode->haloid;

/* 	  fprintf(stderr,"in ensure_same_snapshot: thisnode.snapshot =%8hd fofnode.snapshot = %8hd\n",thisnode->snapshot,fofnode->snapshot); */
/* 	  fprintf(stderr,"fofnode->haloid = %14"STR_FMT" fofnode->FofHalo->haloid = %14"STR_FMT" fofnode->nodeloc=%14"STR_FMT" \n", */
/* 			  fofnode->haloid,fofnode->FofHalo->haloid,fofnode->nodeloc); */

/* 	  fprintf(stderr,"thisnode->haloid = %14"STR_FMT" thisnode->FofHalo->haloid = %14"STR_FMT" thisnode->nodeloc=%14"STR_FMT"\n", */
/* 			  thisnode->haloid,thisnode->FofHalo->haloid,thisnode->nodeloc); */

    while(1) {
      while(fofnode != NULL && fofnode->snapshot < thisnode->snapshot && fofnode->haloid == fofnode_haloid)
	fofnode = fofnode->Parent;

      if((fofnode == NULL) || (fofnode != NULL && fofnode->haloid != fofnode_haloid)) {
	fofnode = NULL;
	break;
      } else {
	if(thisnode->snapshot == fofnode->snapshot) {
	  //thisnode and fofnode are now the same snapshot..break;
	  break;
	} else {
	  //now fofnode is at a later snapshot..increment thisnode until it reaches fofnode
	  while(thisnode != NULL && thisnode->snapshot < fofnode->snapshot && thisnode->haloid == thisnode_haloid) {
	    thisnode = thisnode->Parent;
	  }

	  if(thisnode==NULL || (thisnode != NULL && thisnode->haloid != thisnode_haloid)) {
	    thisnode=NULL;
	    break;
	  }

	  if(thisnode != NULL && thisnode->snapshot == fofnode->snapshot)
	    break;
	}
      }

      if(thisnode == NULL || fofnode == NULL) {
	fprintf(stderr,"in ensure_same_snapshot: failed to get to same snapshots. (at least) one of the node pointers is now NULL \n");
	break;
      }
    }

    /* 	  if(thisnode != NULL && fofnode != NULL && thisnode->snapshot == fofnode->snapshot) */
    /* 		fprintf(stderr,"success: thisnode.snapshot is equal to fofnode.snapshot  = %8hd\n",thisnode->snapshot); */

  }

}
