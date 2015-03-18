#include "defs.h"
#include "genplotdata.h"
#include "read_param.h"
#include "utils.h"
#include "impulse.h"

void print_header_for_massloss(FILE *fp)
{
  fprintf(fp,"#############################################################################################################################################################\n");
  fprintf(fp,"# Snapshot         GroupNum          Haloid         Mtot      ParentLevel       ContainerNum        Nsub         FofHaloId       FofMtot       Rcontainersep \n");
  fprintf(fp,"#    i                l                 l             d           i                   l               l               l              d                 d     \n");
  fprintf(fp,"#############################################################################################################################################################\n");
}



void output_plot_data(struct node_data * tree[],int64 *Ngroups)
{

  FILE *fp=NULL;
  FILE *fp1=NULL;
  FILE *fp2=NULL;
  FILE *fp3=NULL;
  char outfname[MAXLEN];
  double Mthresh = 100.0; /* 1d12 Msun */
  float Mratio  = 0.3;
  const float MinMassRatio = 0.01;
  struct node_data *thisnode,*BaseNode;
  int64 haloid=0;
  float separation;
  double fofmtot;
  int64 Nmissing=0;
  int tag=0;

  const double Mmin = 1e-2;
  const double log10_Mmin=log10(Mmin);
  const double Mmax = 3e4;
  const int Nbins1 = 100;

  double binsize = (log10(Mmax) - log10_Mmin)/(double)Nbins1;
  int bins_noprev_fof[Nbins1],bins_noprev_sub[Nbins1],bins_nonext_fof[Nbins1],bins_nonext_sub[Nbins1],bins_none_fof[Nbins1],bins_none_sub[Nbins1];
  int64 bins_ngroups_sub[Nbins1],bins_ngroups_fof[Nbins1];
  int index;
  int64 nonext,noprev,none;

  int64 fofnum;
  struct node_data *fofnode,*firstnode,*secondnode;
  float rsep,vsep;
  double Impulse=0.0;
  double DelPotPrim=0.0, DelPotSec=0.0;

  int64 jgroup=0;
  
  my_snprintf(outfname,MAXLEN,"%s/MW_massloss.txt",PARAMS.OUTPUT_DIR);
  fp = my_fopen(outfname,"w");
  fprintf(fp,"#   Mthresh = %12.4g   Mratio = %f \n",Mthresh,Mratio);
  print_header_for_massloss(fp);

  my_snprintf(outfname,MAXLEN,"%s/progenitorless_subhalos.txt",PARAMS.OUTPUT_DIR);
  fp1 = my_fopen(outfname,"w");
  fprintf(fp1,"########################################\n");
  fprintf(fp1,"# snapshot        N_without_progenitors \n");
  fprintf(fp1,"#     i                   l             \n");
  fprintf(fp1,"########################################\n");

  my_snprintf(outfname,MAXLEN,"%s/spurioushalos.txt",PARAMS.OUTPUT_DIR);
  fp2 = my_fopen(outfname,"w");
  fprintf(fp2,"#################################################################\n");
  fprintf(fp2,"# snapshot        noprev        nonext       none       Ngroups  \n");
  fprintf(fp2,"#    i               l             l          l            l     \n");
  fprintf(fp2,"#################################################################\n");

  my_snprintf(outfname,MAXLEN,"%s/spurioushalos_bins.txt",PARAMS.OUTPUT_DIR);
  fp3 = my_fopen(outfname,"w");
  fprintf(fp3,"###########################################################################################################################################################\n");
  fprintf(fp3,"# snapshot       log_mass          noprev_fof         noprev_sub         nonext_fof       nonext_sub       none_fof    none_sub     totn_fof       totn_sub\n");
  fprintf(fp3,"#    i              d                  l                   l                 l                 l               l           l            l              l   \n");
  fprintf(fp3,"###########################################################################################################################################################\n");

  for(int isnapshot=PARAMS.MIN_SNAPSHOT_NUM;isnapshot<=PARAMS.MAX_SNAPSHOT_NUM;isnapshot++) {
    Nmissing = 0;
    nonext = 0;
    noprev = 0;
    none   = 0;
    
    if(Ngroups[isnapshot] > 0) {
      for(int i=0;i<Nbins1;i++)	{
				bins_noprev_fof[i] = 0;
				bins_noprev_sub[i] = 0;
				bins_nonext_fof[i] = 0;
				bins_nonext_sub[i] = 0;
				bins_none_fof[i] = 0;
				bins_none_sub[i] = 0;
				bins_ngroups_sub[i]=0;
				bins_ngroups_fof[i]=0;
      }
      
      BaseNode = tree[isnapshot];
      for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup++) {
				thisnode = &BaseNode[igroup];
				if(thisnode->VisitedForMassloss == 0 && (thisnode->Mtot/Mthresh) > Mratio && (thisnode->Mtot/Mthresh) < (1.0+Mratio)) {
					haloid = thisnode->haloid;
					while(thisnode != NULL && thisnode->haloid == haloid) {
						separation = get_separation_between_centres(thisnode,thisnode->ContainerHalo);
	    
						fofmtot = thisnode->FofHalo->Mtot;
						for(int64 gnum=thisnode->FofHalo->nodeloc+1;gnum < Ngroups[isnapshot];gnum++) {
							if(BaseNode[gnum].isFof == 1)
								break;
	      
							fofmtot += BaseNode[gnum].Mtot;
						}
	    
						fprintf(fp,"%6d     %14"STR_FMT"   %14"STR_FMT"   %12.4g   %10hd    %14"STR_FMT"  %12"STR_FMT"    %14"STR_FMT "  %12.4g    %14.4g  \n",
										thisnode->snapshot,thisnode->nodeloc,thisnode->haloid,thisnode->Mtot,thisnode->ParentLevel,thisnode->ContainerId,
										thisnode->Nsub,thisnode->FofHalo->haloid,fofmtot,separation);
	    
						thisnode->VisitedForMassloss = 1;
						thisnode = thisnode->Parent;
					}
				  
				}
			 
				thisnode  = &BaseNode[igroup];
				if(thisnode->isFof != 1 && thisnode->BigChild == NULL )
					Nmissing++;
	
				index = (int) (  (log10(thisnode->Mtot) - log10_Mmin)/binsize);
				if (index < Nbins1 && index >=0) {
					if(thisnode->Parent == NULL && thisnode->BigChild != NULL && thisnode->snapshot != (PARAMS.MAX_SNAPSHOT_NUM)) {
						nonext++;
						if(thisnode->isFof==1)
							bins_nonext_fof[index]++;
						else
							bins_nonext_sub[index]++;
					}
	  
					if(thisnode->BigChild == NULL && thisnode->Parent != NULL && thisnode->snapshot != PARAMS.MIN_SNAPSHOT_NUM) {
						noprev++;
						if(thisnode->isFof==1)
							bins_noprev_fof[index]++;
						else
							bins_noprev_sub[index]++;
					}
	  
					if(thisnode->BigChild == NULL && thisnode->Parent == NULL && (thisnode->snapshot != PARAMS.MIN_SNAPSHOT_NUM || thisnode->snapshot != (PARAMS.MAX_SNAPSHOT_NUM))) {
						none++;
						if(thisnode->isFof==1)
							bins_none_fof[index]++;
						else
							bins_none_sub[index]++;
					}
	  
					if(thisnode->isFof==1)
						bins_ngroups_fof[index]++;
					else
						bins_ngroups_sub[index]++;
				}
      }
      for(int i=0;i<Nbins1;i++)
				fprintf(fp3,"%12hd    %12.4g   %12d     %12d     %12d     %12d    %12d    %12d   %12"STR_FMT"   %12"STR_FMT"\n",isnapshot,i*binsize+log10_Mmin,
								bins_noprev_fof[i],bins_noprev_sub[i],bins_nonext_fof[i],bins_nonext_sub[i],bins_none_fof[i],bins_none_sub[i],
								bins_ngroups_fof[i],bins_ngroups_sub[i]);
    }
    fprintf(fp1,"%10d    %14"STR_FMT" \n",isnapshot,Nmissing);
    fprintf(fp2,"%10d  %12"STR_FMT"    %12"STR_FMT"   %12"STR_FMT"   %12"STR_FMT" \n",isnapshot,noprev,nonext,none,Ngroups[isnapshot]);
    
  }
  
  fclose(fp);
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);


/*   my_snprintf(outfname,MAXLEN,"%s%s%s","rm -f ",PARAMS.OUTPUT_DIR,"/halo_subhalo_flybys.txt"); */
/*   fprintf(stderr,"executing system command `%s' \n",outfname); */
/*   system(outfname);  */

  my_snprintf(outfname,MAXLEN,"%s/halo_subhalo_flybys.txt",PARAMS.OUTPUT_DIR);
  fp = my_fopen(outfname,"w");
  print_header_for_subsub_flyby(fp);

  /* Now I am going to locate the subhalos that have undergone flybys */
  for(int isnapshot=PARAMS.MIN_SNAPSHOT_NUM;isnapshot<PARAMS.MAX_SNAPSHOT_NUM;isnapshot++) {
    if(Ngroups[isnapshot] > 0) {
      BaseNode = tree[isnapshot];
      fofnum = 0;
      fofnode = &BaseNode[fofnum];
      for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup+=fofnode->Nsub) {
				fofnode = &BaseNode[igroup];
				fofnum = fofnode->nodeloc;
				if(fofnode->isFof != 1) {
					fprintf(stderr,"Error: While looking for sub-sub flybys -- could not locate fof halo\n");
					fprintf(stderr,"snapshot = %d  fofnodeloc = %"STR_FMT"   \n",isnapshot,fofnum);
					exit(EXIT_FAILURE);
				}
	
				for(int64 first=fofnum;first<fofnum+fofnode->Nsub-1;first++) {
					for(int64 second=first+1;second<fofnum+fofnode->Nsub;second++) {
						firstnode  = &BaseNode[first];
						secondnode = &BaseNode[second];
						rsep = 0.0;
						vsep = 0.0;
						if(secondnode->Mtot/firstnode->Mtot >= MinMassRatio) {
							if(check_for_flyby(firstnode,secondnode,&rsep,&vsep) == 1) {
								tag = -1;
								/* sub-> sub-sub->sub flyby */
								if(secondnode->BigChild != NULL && secondnode->BigChild->ParentLevel >= 2 && secondnode->BigChild->ParentLevel < secondnode->ParentLevel &&
									 secondnode->Parent != NULL && secondnode->Parent->haloid == secondnode->haloid && secondnode->Parent->ParentLevel < secondnode->ParentLevel &&
									 firstnode->isFof == 0)
									tag = 3;
								else {
									/* potentially no change of hierarchy at this snapshot */
									if(secondnode->BigChild != NULL && secondnode->BigChild->ParentLevel >= 2 && secondnode->BigChild->ParentLevel <= secondnode->ParentLevel &&
										 secondnode->Parent != NULL && secondnode->Parent->haloid == secondnode->haloid && secondnode->Parent->ParentLevel <= secondnode->ParentLevel &&
										 firstnode->isFof == 0)
										tag = 5;
		  
									if(secondnode->ContainerHalo == firstnode && firstnode->isFof==1)
										tag = 2;
		  
								}
		
								/* Previously impulse.c assumed that the first arg. was the more massive one. 
									 This assumption might no longer be valid.
								*/
		
								Impulse=get_impulse(secondnode,firstnode,rsep,vsep);
								DelPotPrim=get_internal_delpot_prim(secondnode,firstnode,rsep,vsep);
								DelPotSec=get_internal_delpot_sec(secondnode,firstnode,rsep,vsep);
		
								fprintf(fp,"%10d   %14"STR_FMT"  %14"STR_FMT"  %14.6g    %14.6g     %12.4f    %12.4f   %14"STR_FMT"   %14"STR_FMT "   %10.4g    %10.4g   %10.4g  %12d\n",
												firstnode->snapshot,firstnode->nodeloc,secondnode->nodeloc,firstnode->Mtot,secondnode->Mtot,rsep,vsep,
												firstnode->ContainerId,secondnode->ContainerId,Impulse,DelPotPrim,DelPotSec,tag);
							}
						}
					}
				}
      }
    }
  }
  if(fp != NULL)
    fclose(fp);

  
  /* Now search for all nearby FOFs */
  /*   struct node_data *nextfofnode=NULL; */
  /*   float pre_fac = 3.0; */
  /*   float r_thresh = 0.0; */
  /*   double vdotv=0.0; */
  /*   double MinMassThresh = 1.0; */
  /*   my_snprintf(outfname,MAXLEN,"%s/fof_and_fof_flybys.txt",PARAMS.OUTPUT_DIR); */
  /*   fp = my_fopen(outfname,"w"); */
  /*   print_header_for_fof_fof_flyby(fp); */
  
  /*   /\* Now I am going to locate the subhalos that have undergone flybys *\/ */
  /*   for(int isnapshot=PARAMS.MIN_SNAPSHOT_NUM;isnapshot<PARAMS.MAX_SNAPSHOT_NUM-1;isnapshot++) */
  /* 	{ */
  /* 	  if(Ngroups[isnapshot] > 0) */
  /* 		{ */
  /* 		  BaseNode = tree[isnapshot]; */
  /* 		  fofnum = 0; */
  /* 		  fofnode = &BaseNode[fofnum]; */
  
  /* 		  for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup+=fofnode->Nsub) */
  /* 			{ */
  /* 			  fofnode = &BaseNode[igroup]; */
  /* 			  fofnum = fofnode->nodeloc; */
  /* 			  jgroup = igroup + fofnode->Nsub; */
  /* 			  while(jgroup < Ngroups[isnapshot]) */
  /* 				{ */
  /* 				  nextfofnode = &BaseNode[jgroup]; */
  /* 				  if(fofnode->Mtot > MinMassThresh || nextfofnode->Mtot > MinMassThresh) */
  /* 					{ */
  /* 					  if(fofnode->BigChild != NULL && fofnode->Parent != NULL && fofnode->Parent->haloid == fofnode->haloid && */
  /* 						 nextfofnode->BigChild != NULL && nextfofnode->Parent !=NULL && nextfofnode->Parent->haloid == nextfofnode->haloid) */
  /* 						{ */
/* 						  rsep=get_separation_between_centres(fofnode,nextfofnode); */
/* 						  r_thresh = pre_fac * ((fofnode->Rvir > nextfofnode->Rvir ) ? fofnode->Rvir:nextfofnode->Rvir); */
  
/* 						  if(rsep < r_thresh && find_fof_flyby_future(fofnode,nextfofnode) == 1) */
/* 							{ */
  
/* 							  vsep = 0.0; */
/* #ifdef GET_MEANVEL */
/* 							  for(int i=0;i<3;i++) */
/* 								vsep += pow( (fofnode->meanvel[i] - nextfofnode->meanvel[i]),2.0); */
/* #else */
/* #ifdef GET_GROUPVEL */
  
/* #endif */
							  
/* #endif */
							  
/* 							  vsep = sqrt(vsep); */
/* 							  vdotv = get_vdotv(fofnode,nextfofnode); */
							  
/* 							  Impulse=get_impulse(fofnode,nextfofnode,rsep,vsep); */

/* 							  /\*for get_delpot functions, first node should be secondary => nextfofnode should appear first *\/ */
/* 							  DelPotPrim=get_external_delpot_prim(nextfofnode,fofnode,rsep,vsep); */
/* 							  DelPotSec=get_external_delpot_sec(nextfofnode,fofnode,rsep,vsep); */
							  
/* 							  fprintf(fp,"%10d   %14"STR_FMT"   %14"STR_FMT"   %14.6g   %14.6g    %12.4f     %12.4f    %14.6g  %14.6lg   %12.4lg   %12.4lg  %12.4lg %12.4g\n", */
/* 									  isnapshot,fofnum,jgroup,fofnode->Mtot,nextfofnode->Mtot,rsep,vsep,fofnode->Rvir,nextfofnode->Rvir,Impulse, */
/* 									  DelPotPrim,DelPotSec,vdotv); */

/* 							  fflush(fp); */
/* 							} */
/* 						} */
/* 					} */
/* 				  jgroup += nextfofnode->Nsub; */
/* 				} */
/* 			} */
/* 		} */
/* 	} */

/*   fclose(fp); */



  /*
		Find half mass assembly history for halos that exist at some supplied set of 
		snapshots. 
	
		Oddly enough I am going to define all the required variables here.
	
  */

/*   int NOUTPUTS = 6; */
/*   int Snapshots[] = {107,60,41,26}; */
/*   int Snapshots[] = {111,64,45,30};//for the 50 Mpc 1024^3 */
/*   int Snapshots[] = {60,41,26}; */
/*   int Snapshots[] = {44,30,20,10}; */
/*   int Snapshots[] = {21,15,10,5}; */
/*   int Snapshots[] = {14,10,6,3}; */
/*   int Snapshots[] = {10,7,5,2}; */
/*   int Snapshots[] = {2,1}; */
/*   int Snapshots[] = {8,6,4,2}; */
/*   int Snapshots[] = {53,30,20,12}; */
  int Snapshots[] = {141,131,111,83,63,25};//for the new 512^3 phase sims.
/*   int Snapshots[] = {79,65,55,41,31,12};//for the new 512^3 phase sims - skip one */
/*   int Snapshots[] = {53,43,36,27,21,8};//for the new 512^3 phase sims. */
/*   int Snapshots[]  = {61,54,47,41,26};//for the 10Mpc 2lpt sims */
  int Nbins2[] = {40,40,40,40,40,20}; /* these choices produce at least one halo from the min to the max (of the halo bins that
																				 are occupied. Nbins-1 itself and lower ones might be 0. */
  
/*   double mass_limits[] = {6.7207024e9,2e14}; */
/*   double log_mass_limits[] = {-1.07567532,4.3010}; /\* in units of 1d10 for the 1024^3, 50 Mpc*\/ */
  double log_mass_limits[] = {log10(Mmin),log10(Mmax)}; /* in units of 1d10 for the 1024^3, 50 Mpc*/
/*   double log_mass_limits[] = {-2.97190, 0.904237}; /\* in units of 1d10*\/ */
  double log_mass_bin = 0.0;
/*   int64 *h=NULL; */
  int64 *nhalos=NULL;
  double *halfz=NULL,*halfz_sigma=NULL;
  int NOUTPUTS =(int) (sizeof(Snapshots)/sizeof(int));

  for (int i=0;i<NOUTPUTS;i++)	{
    if(Snapshots[i] >= PARAMS.MIN_SNAPSHOT_NUM && Snapshots[i] <= PARAMS.MAX_SNAPSHOT_NUM) {
      log_mass_bin = (log_mass_limits[1]-log_mass_limits[0])/(double)Nbins2[i];
      halfz = my_malloc(sizeof(*halfz),Nbins2[i]+1);
      halfz_sigma = my_malloc(sizeof(*halfz_sigma),Nbins2[i]+1);
      nhalos = my_malloc(sizeof(*nhalos),Nbins2[i]+1);
      for(int j=0;j<=Nbins2[i];j++) {
				/* 		  h[j] = 0; */
				halfz[j] = 0.0;
				halfz_sigma[j] = 0.0;
				nhalos[j] = 0;
      }
      
      BaseNode = tree[Snapshots[i]];
      for(int64 igroup=0;igroup<Ngroups[Snapshots[i]];igroup++) {
				thisnode = &BaseNode[igroup];
				firstnode = &BaseNode[igroup];
				if(log10(thisnode->Mtot) >= log_mass_limits[0] && log10(thisnode->Mtot) <= log_mass_limits[1]) {
					while(thisnode != NULL && thisnode->Mtot > 0.5*firstnode->Mtot) {
						thisnode = thisnode->BigChild;
					}
	  
					if(thisnode != NULL) {
						index = (int) ((log10(firstnode->Mtot)-log_mass_limits[0])/log_mass_bin);
						if (index <= Nbins2[i]) {
							nhalos[index]++;
							halfz[index] += thisnode->z;
							halfz_sigma[index] = halfz_sigma[index] + thisnode->z*thisnode->z;
						} else {
							fprintf(stderr,"index  = %d is outside nbins = %d \n",index,Nbins2[i]);
							fprintf(stderr,"firstnode->mtot = %g  log_mass_bin = %g \n",firstnode->Mtot,log_mass_bin);
						}
					}
	  
				}
      }
      
      /* output the data */
      my_snprintf(outfname,MAXLEN,"%s/halfmasshistory_%03d.txt",PARAMS.OUTPUT_DIR,Snapshots[i]);
      fp = my_fopen(outfname,"w");
      print_header_for_halfmass(fp);
      
      for(int j=0;j<=Nbins2[i];j++) {
				if(nhalos[j] > 0) {
					const double mean_halfz = halfz[j]/(double) nhalos[j];
					const double mean_sqr_halfz = halfz_sigma[j]/(double)nhalos[j];
					fprintf(fp,"%8d   %10d   %10"STR_FMT" %14.4g  %14.4g   %14.4g  \n",Snapshots[i],j,nhalos[j],
									log_mass_limits[0] + log_mass_bin*j,
									mean_halfz, sqrt(mean_sqr_halfz - mean_halfz*mean_halfz));
				}
      }
      
      fclose(fp);
      
      /* free memory */
      /* 	  my_free((void **) h); */
      my_free((void **) &halfz);
		  my_free((void **) &halfz_sigma);
		  my_free((void **) &nhalos);
    }
  }
  
  
  
  /*Get the t0.04 for the halos following Zhao et al 2009. Really, all that is different
    between this section and the previous is that the inner while loop here
    has a multiplicative constant of 0.04 (instead of 0.5)*/

  for (int i=0;i<NOUTPUTS;i++)	{
    if(Snapshots[i] <= PARAMS.MAX_SNAPSHOT_NUM && Snapshots[i] >= PARAMS.MIN_SNAPSHOT_NUM) {
      log_mass_bin = (log_mass_limits[1]-log_mass_limits[0])/(double)Nbins2[i];
      halfz = my_malloc(sizeof(*halfz),Nbins2[i]+1);
      halfz_sigma = my_malloc(sizeof(*halfz_sigma),Nbins2[i]+1);
      nhalos = my_malloc(sizeof(*nhalos),Nbins2[i]+1);
      for(int j=0;j<=Nbins2[i];j++) {
				/* 		  h[j] = 0; */
				halfz[j] = 0.0;
				halfz_sigma[j] = 0.0;
				nhalos[j] = 0;
      }
      
      BaseNode = tree[Snapshots[i]];
      for(int64 igroup=0;igroup<Ngroups[Snapshots[i]];igroup++) {
				thisnode = &BaseNode[igroup];
				firstnode = &BaseNode[igroup];
				if(log10(thisnode->Mtot) >= log_mass_limits[0] && log10(thisnode->Mtot) <= log_mass_limits[1]) {
					while(thisnode != NULL && thisnode->Mtot > 0.04*firstnode->Mtot) {
						thisnode = thisnode->BigChild;
					}
	  
					if(thisnode != NULL) {
						index = (int) ((log10(firstnode->Mtot)-log_mass_limits[0])/log_mass_bin);
						if (index <= Nbins2[i]) {
							nhalos[index]++;
							halfz[index] += thisnode->z;
							halfz_sigma[index] = halfz_sigma[index] + thisnode->z*thisnode->z;
						} else {
							fprintf(stderr,"index  = %d is outside nbins = %d \n",index,Nbins2[i]);
							fprintf(stderr,"firstnode->mtot = %g  log_mass_bin = %g \n",firstnode->Mtot,log_mass_bin);
						}
					}
	  
				}
      }
		  
      /* output the data */
      my_snprintf(outfname,MAXLEN,"%s/mah_t0.04_zhao09_%03d.txt",PARAMS.OUTPUT_DIR,Snapshots[i]);
      fp = my_fopen(outfname,"w");
      print_header_for_halfmass(fp);
      
      for(int j=0;j<=Nbins2[i];j++) {
				if(nhalos[j] > 0) {
					fprintf(fp,"%8d   %10d   %10"STR_FMT" %14.4g  %14.4g   %14.4g  \n",Snapshots[i],j,nhalos[j],
									log_mass_limits[0] + log_mass_bin*j,
									halfz[j]/(double)nhalos[j],sqrt(halfz_sigma[j]/(double)nhalos[j] - halfz[j]/(double)nhalos[j]*halfz[j]/(double)nhalos[j]));
				}
      }
      
      fclose(fp);
      
      /* free memory */
      /* 	  my_free((void **) h); */
      my_free((void **) &halfz);
      my_free((void **) &halfz_sigma);
      my_free((void **) &nhalos);
      
    }
  }


  /* Generate the data for the Fakhouri & Ma plot */
  my_snprintf(outfname,MAXLEN,"%s/fof_and_fof_fma_mergers.txt",PARAMS.OUTPUT_DIR);
  fp = my_fopen(outfname,"w");
  print_header_for_fma(fp);
  float destructionz,destructionradius;
  for(int isnapshot=PARAMS.MIN_SNAPSHOT_NUM;isnapshot<=PARAMS.MAX_SNAPSHOT_NUM;isnapshot++) {
    if(Ngroups[isnapshot] > 0) {
      BaseNode = tree[isnapshot];
      fofnum = 0;
      fofnode = &BaseNode[fofnum];
      for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup+=fofnode->Nsub) {
				fofnode = &BaseNode[igroup];
				jgroup = igroup + 1;
				while((jgroup-igroup) < fofnode->Nsub && jgroup < Ngroups[isnapshot]) {
					thisnode = &BaseNode[jgroup];
					if(thisnode->BigChild != NULL && thisnode->BigChild->isFof==1 && thisnode->BigChild->haloid != fofnode->haloid
						 && fofnode->BigChild != NULL && fofnode->BigChild->isFof==1 && fofnode->BigChild->snapshot == thisnode->BigChild->snapshot
						 && fofnode->BigChild->snapshot == (fofnode->snapshot-1) && thisnode->FofHalo == fofnode) /*the last check is to protect against my stupidity*/
					{
						/* 					  tag = get_tag(thisnode,thisnode->FofHalo,&destructionz,&destructionradius); */
						//MS 7th Dec, 2011 - updated to use the new tagging system that accounts for 
						//subhalos missing snapshots
						double destruction_subtmot=-1.0;
						double destruction_parentmtot=-1.0;
						int destruction_snapshot=-1;
						double destruction_mstar=-1.0;

						tag = get_tag(thisnode,thisnode->FofHalo,&destructionz,&destructionradius,&destruction_snapshot,&destruction_subtmot,&destruction_parentmtot,&destruction_mstar);
						if(tag==0 || tag==2) /*merger tags*/
						{
							fprintf(fp," %6d   %14.4e  %14"STR_FMT "  %14"STR_FMT"  %14"STR_FMT"  %14"STR_FMT"   %14.4lf    %14.4lf   %14.4lf\n",
											thisnode->BigChild->snapshot,REDSHIFT[thisnode->BigChild->snapshot],
											fofnode->BigChild->nodeloc,thisnode->BigChild->nodeloc,
											fofnode->BigChild->haloid,thisnode->BigChild->haloid,
											fofnode->BigChild->BoundFofMtot,thisnode->BigChild->BoundFofMtot,
											fofnode->BoundFofMtot);
						}
					}
					jgroup++;
				}
      }
    }
  }
  fclose(fp);

  /* my_snprintf(outfname,MAXLEN,"%s/halo_distribution_from_list.txt",PARAMS.OUTPUT_DIR); */
  /* fp1 = my_fopen(outfname,"w"); */
  /* fprintf(fp1,"#########################################################################################################################################################################################################################################\n"); */
  /* fprintf(fp1,"# Snapshot        FofHaloId     Sub-ContainerHaloId       SubhaloId      FofMtot      BoundFofMtot    ContainerMtot      SubMtot     R_rel      V_rel      FofRhalf   FofRvir    FofRvir_anyl    rp      vproj    fof_mstar   sub_mstar \n"); */
  /* fprintf(fp1,"#    i                l                   l                   l             d             d                d               d          d          d            d         d             d          d         d          d          d      \n"); */
  /* fprintf(fp1,"########################################################################################################################################################################################################################################\n"); */

	  struct node_data *subnode;
  
  /* my_snprintf(outfname,MAXLEN,"%s/list_of_halos_inside_search_radius.txt",PARAMS.OUTPUT_DIR); */
  /* fp = my_fopen(outfname,"r"); */
  /* char buffer[MAXLINESIZE]; */
  /* int nread,start_snapshot; */
  /* int64 primary_haloid,sub_haloid; */
  /* struct node_data *subnode; */
  /* while(fgets(buffer,MAXLINESIZE,fp) != NULL) { */
  /*   if(buffer[0] == '#') */
  /*     continue; */

  /*   nread = sscanf(buffer, "%d %"STR_FMT" %"STR_FMT" ",&start_snapshot,&primary_haloid,&sub_haloid); */
  /*   if(nread == 3) { */
  /*     thisnode = NULL; */
  /*     if(start_snapshot <= PARAMS.MAX_SNAPSHOT_NUM && start_snapshot >= PARAMS.MIN_SNAPSHOT_NUM) { */
	/* 			BaseNode = tree[start_snapshot]; */
	/* 			for(int64 igroup=0;igroup < Ngroups[start_snapshot];igroup++) { */
	/* 				thisnode = &BaseNode[igroup]; */
	/* 				if(thisnode->haloid == primary_haloid) */
	/* 					break; */
	/* 			} */

	/* 			for(int64 igroup=0;igroup < Ngroups[start_snapshot];igroup++) { */
	/* 				subnode = &BaseNode[igroup]; */
	/* 				if(subnode->haloid == sub_haloid) */
	/* 					break; */
	/* 			} */
	

	/* 			if(thisnode != NULL && thisnode->haloid == primary_haloid && subnode != NULL && subnode->haloid == sub_haloid) { */
	/* 				//Ok - so found both the primary and the secondary halos */
	/* 				while(thisnode != NULL && subnode != NULL) { */
	/* 					//output data about the relative separations + relative velocities */
	/* 					double vel=0.0; */
	/* 					double rp=0.0,vproj=0.0; */
	/* 					double scale_factor=1.0/(1.0+thisnode->z); */
	/* 					double dist=get_separation_between_centres(thisnode,subnode); //dist is in co-moving units */
	    
	/* 					for(int k=0;k<3;k++)  */
	/* 						vel += (subnode->meanvel[k] - thisnode->meanvel[k])*(subnode->meanvel[k] - thisnode->meanvel[k]); */

	/* 					rp += periodic(subnode->xcen - thisnode->xcen)*periodic(subnode->xcen - thisnode->xcen); */
	/* 					rp += periodic(subnode->ycen - thisnode->ycen)*periodic(subnode->ycen - thisnode->ycen); */
	/* 					rp = sqrt(rp); */
	    
	/* 					for(int k=0;k<2;k++) */
	/* 						vproj += (subnode->meanvel[k] - thisnode->meanvel[k])*(subnode->meanvel[k] - thisnode->meanvel[k]); */

	/* 					vproj = sqrt(vproj); */
	/* 					vproj += PARAMS.COSMO->H0 * epeebles(thisnode->z) * rp * scale_factor;//now adding the Hubble flow -> vel is physical velocity */
	    

	/* 					vel = sqrt(vel);//Only peculiar velocity -> Hubble flow not added in  */
	/* 					vel += PARAMS.COSMO->H0 * epeebles(thisnode->z) * dist * scale_factor;//now adding the Hubble flow -> vel is physical velocity */

	/* 					fprintf(fp1,"%10d   %10"STR_FMT"  %16"STR_FMT"  %20"STR_FMT"    %12.6g  %12.6g  %14.6g  %14.6g  %8.2lf  %8.2lf  %10.4g  %10.4g %10.4g  %10.4g  %10.4g  %10.4g  %10.4g\n", */
	/* 									thisnode->snapshot,thisnode->haloid,subnode->ContainerHalo->haloid,sub_haloid, */
	/* 									thisnode->Mtot,thisnode->BoundFofMtot,subnode->ContainerHalo->Mtot,subnode->Mtot,dist,vel, */
	/* 									thisnode->Rhalf,thisnode->Rvir,thisnode->Rvir_anyl,rp,vproj,thisnode->Mstar,subnode->Mstar); */
	/* 					fflush(fp1); */
	    
	/* 					thisnode = thisnode->BigChild; */
	/* 					subnode = subnode->BigChild; */
	/* 					while(thisnode != NULL && subnode != NULL && thisnode->snapshot != subnode->snapshot) { */
	/* 						if(subnode->snapshot < thisnode->snapshot) { */
	/* 							thisnode = thisnode->BigChild; */
	/* 						} else { */
	/* 							subnode = subnode->BigChild; */
	/* 						} */
	/* 					} */
	/* 				} */
	/* 			} */
  /*     } else { */
	/* 			fprintf(stderr,"WARNING: start_snapshot = %d not inside the snapshot limits [%d,%d] ..ignoring\n",start_snapshot,PARAMS.MIN_SNAPSHOT_NUM, PARAMS.MAX_SNAPSHOT_NUM); */
  /*     } */
  /*   } else { */
  /*     fprintf(stderr,"WARNING: Could not parse three values for [start_snapshot, primary_haloid, sub_haloid] \n"); */
  /*   } */
  
  /* } */
  /* fclose(fp1); */
  /* fclose(fp); */


  


  //Now the plot to make the rsep/vsep probability at the some pre-determined point in the simulation.
  //The basic idea is to target halos of a certain mass and see the distribution of first infall vs flybys
  //in two-D plots (rsep/rvir , vsep/vvir). 
  const int target_snapshots[]={112,87,74,64};
  int N_target_redshifts=sizeof(target_snapshots)/sizeof(int);
  const double Min_Target_Halo_Mass[]={1e2,1e3,1e4};//10^13 Msun/h
  const double Max_Target_Halo_Mass[]={2e2,2e3,3e4};//10^13 Msun/h
  int N_target_mass=sizeof(Min_Target_Halo_Mass)/sizeof(double);
  struct node_data *tmp_fofnode=NULL;
  int64 subhaloid;
  const double Time_Window_Factor=1.0;//Look for interactions in the window: [t_now, Time_Window_Factor*t_now]. Time_Window_Factor=1.0 implies *all* interactions regardless of when they occurred
  for(int itarget=0;itarget<N_target_redshifts;itarget++) {
    if(target_snapshots[itarget] >= PARAMS.MIN_SNAPSHOT_NUM && target_snapshots[itarget] <= PARAMS.MAX_SNAPSHOT_NUM) {
      int last_snapshot=target_snapshots[itarget];
      int earliest_snapshot;
      double Time_Window=Time_Window_Factor*PARAMS.Age[last_snapshot];
      double Earliest_Age=PARAMS.Age[last_snapshot] - Time_Window;

      for(earliest_snapshot=PARAMS.MIN_SNAPSHOT_NUM;earliest_snapshot<=PARAMS.MAX_SNAPSHOT_NUM;earliest_snapshot++) {
				if(PARAMS.Age[earliest_snapshot] >= Earliest_Age) {
					break;
				}
      }      

      if(earliest_snapshot > PARAMS.MIN_SNAPSHOT_NUM) {
				earliest_snapshot--;
      }

      fprintf(stderr,"Age[%d] = %e Time_window=%lf Earliest_Age=%lf Age[%d] = %lf \n",last_snapshot,PARAMS.Age[last_snapshot],Time_Window,Earliest_Age,earliest_snapshot,PARAMS.Age[earliest_snapshot]);
      
      //Find the earliest snapshot that has an age lower than Earliest_age
      for(int imass=0;imass<N_target_mass;imass++) {
				//Mass-bins
				my_snprintf(outfname,MAXLEN,"%s/halo_distribution_near_bighalos_imass_%02d_%03d.txt",PARAMS.OUTPUT_DIR,imass,last_snapshot);
				fp = my_fopen(outfname,"w");
				fprintf(fp,"# Min mass = %lf max. mass = %lf target redshift = %f \n",Min_Target_Halo_Mass[imass],Max_Target_Halo_Mass[imass],REDSHIFT[last_snapshot]);
				fprintf(fp,"###################################################################################################################################################################################################################################################################################################################\n");
				fprintf(fp,"# Snapshot        FofHaloId     ContainerHaloId       SubhaloId      FofMtot      BoundFofMtot    ContainerMtot      SubMtot     R_rel      V_rel      Tag     DestructionZ   FofRhalf   FofRvir    FofRvir_anyl       rp      vproj        SubInfallMass    InfallSnapshot    Sub_StellarMass      FofStellarMass \n");
				fprintf(fp,"#    i                l                l                   l             d             d                d               d          d          d         i         d               d         d             d            d         d               d                   i               d                    d        \n");
				fprintf(fp,"###################################################################################################################################################################################################################################################################################################################\n");

				BaseNode = tree[last_snapshot];
				for(int64 igroup=0;igroup<Ngroups[last_snapshot];igroup++) {
					thisnode = &BaseNode[igroup];
					if(thisnode->isFof==1 && thisnode->BoundFofMtot >= Min_Target_Halo_Mass[imass] && thisnode->BoundFofMtot < Max_Target_Halo_Mass[imass]) {
						//Okay so found a halo at the target redshift that satisfies the mass criteria
						while(thisnode->BigChild != NULL) {
							thisnode = thisnode->BigChild;
							if(thisnode->snapshot <= earliest_snapshot) {
								break;
							}
						}

						//Now at the earliest snapshot that encompasses the desired TimeWindow.
						while(thisnode != NULL && thisnode->snapshot <= last_snapshot) {
							if(thisnode->isFof==1) {
								double scale_factor=1.0/(1.0+thisnode->z);
								/* double hubble = get_hubble_at_redshift(thisnode->z,COSMO); */
								for(jgroup=1;jgroup<thisnode->Nsub;jgroup++) {
									tmp_fofnode = thisnode;
									subnode = thisnode + jgroup;
									subhaloid=subnode->haloid;
									double destruction_subtmot=-1.0;
									double destruction_parentmtot=-1.0;
									int destruction_snapshot=-1;
									double destruction_mstar=-1.0;
		  
									tag = get_tag(subnode,subnode->FofHalo,&destructionz,&destructionradius,&destruction_snapshot,&destruction_subtmot,&destruction_parentmtot,&destruction_mstar);

									//check if the interaction is classified as a split only because we are
									//not starting at the beginning of the interaction
									{
										if(tag == 4) {
											struct node_data *tmp1;
											tmp1 = subnode;
		      
											while(tmp1 != NULL && tmp1->isFof==0 && tmp1->FofHalo->haloid == tmp_fofnode->haloid)
												tmp1 = tmp1->BigChild;
		      
											if(tmp1 != NULL && tmp1->isFof==1)
												tag = 3;//the subhalo was a Fof at some point in the past. 
										}
									}

									while(subnode != NULL && subnode->haloid==subhaloid && subnode->snapshot <= last_snapshot) {
										double vel=0.0;
										double dist=0.0;
										double rp=0.0;
										dist  = periodic(subnode->xcen - tmp_fofnode->xcen)*periodic(subnode->xcen - tmp_fofnode->xcen);
										dist += periodic(subnode->ycen - tmp_fofnode->ycen)*periodic(subnode->ycen - tmp_fofnode->ycen);

										rp = sqrt(dist);
		    
										dist += periodic(subnode->zcen - tmp_fofnode->zcen)*periodic(subnode->zcen - tmp_fofnode->zcen);
										dist = sqrt(dist);
		    
										for(int k=0;k<2;k++) 
											vel += (subnode->meanvel[k] - tmp_fofnode->meanvel[k])*(subnode->meanvel[k] - tmp_fofnode->meanvel[k]);

										double vproj = sqrt(vel);//projected vel
										for(int k=2;k<3;k++) 
											vel += (subnode->meanvel[k] - tmp_fofnode->meanvel[k])*(subnode->meanvel[k] - tmp_fofnode->meanvel[k]);

										vel = sqrt(vel);

										//Add in scale_factor * dist * H(z) to get to physical velocity
										vel   += PARAMS.COSMO->H0 * epeebles(thisnode->z) * dist * scale_factor;//now adding the Hubble flow -> vel is physical velocity
										vproj += PARAMS.COSMO->H0 * epeebles(thisnode->z) * rp   * scale_factor;//now adding the Hubble flow -> vproj is projected physical velocity
		    
										fprintf(fp,"%10d   %10"STR_FMT"  %16"STR_FMT"  %20"STR_FMT"    %12.6g  %12.6g  %14.6g  %14.6g  %8.2lf  %8.2lf  %8d   %10.4lf  %10.4g  %10.4g %10.4g %12.4g  %12.4g  %12.4g %10d  %12.4g %12.4g\n",
														subnode->snapshot,tmp_fofnode->haloid,subnode->ContainerHalo->haloid,subhaloid,
														tmp_fofnode->Mtot,tmp_fofnode->BoundFofMtot,subnode->ContainerHalo->Mtot,subnode->Mtot,dist,vel,tag,destructionz,
														tmp_fofnode->Rhalf,tmp_fofnode->Rvir,tmp_fofnode->Rvir_anyl,rp,vproj,subnode->InfallMass,subnode->InfallSnapshot,subnode->Mstar,tmp_fofnode->Mstar);

										subnode=subnode->Parent;
										tmp_fofnode=tmp_fofnode->Parent;

										if(subnode != NULL && tmp_fofnode!=NULL && subnode->snapshot != tmp_fofnode->snapshot)
											ensure_same_snapshot_for_halos(&subnode,&tmp_fofnode);

									}
								}
							} /* else { */
							/* 	//Now the halo is not a FOF */
							/* 	fofnode=thisnode->FofHalo; */
							/* 	for(jgroup=1;jgroup<fofnode->Nsub;jgroup++) { */
							/* 	  subnode=fofnode+jgroup; */
							/* 	  tmp_fofnode = thisnode; */
							/* 	  tag = get_tag(subnode,fofnode,&destructionz,&destructionradius); */
							/* 	  if(subnode->ContainerHalo == thisnode) { */
							/* 	    subhaloid=subnode->haloid; */
							/* 	    while(subnode != NULL && subnode->haloid==subhaloid && subnode->snapshot <= last_snapshot) { */
							/* 	      double vel=0.0; */
							/* 	      double dist=0.0; */
							/* 	      dist  = periodic(subnode->xcen - tmp_fofnode->xcen)*periodic(subnode->xcen - tmp_fofnode->xcen); */
							/* 	      dist += periodic(subnode->ycen - tmp_fofnode->ycen)*periodic(subnode->ycen - tmp_fofnode->ycen); */
							/* 	      dist += periodic(subnode->zcen - tmp_fofnode->zcen)*periodic(subnode->zcen - tmp_fofnode->zcen); */
							/* 	      dist = sqrt(dist); */

							/* 	      for(int k=0;k<3;k++)  */
							/* 		vel += (subnode->meanvel[k] - tmp_fofnode->meanvel[k])*(subnode->meanvel[k] - tmp_fofnode->meanvel[k]); */

							/* 	      vel = sqrt(vel);//Only peculiar velocity -> Hubble flow not added in  */

							/* 	      fprintf(fp,"%10d   %10"STR_FMT"  %16"STR_FMT"  %20"STR_FMT"    %12.6g  %12.6g  %14.6g  %14.6g  %8.2lf  %8.2lf  %8d   %10.4lf  %10.4g  %10.4g %10.4g \n", */
							/* 		      subnode->snapshot,tmp_fofnode->haloid,subnode->ContainerHalo->haloid,subhaloid, */
							/* 		      tmp_fofnode->Mtot,tmp_fofnode->BoundFofMtot,subnode->ContainerHalo->Mtot,subnode->Mtot,dist,vel,tag,destructionz, */
							/* 		      tmp_fofnode->Rhalf,tmp_fofnode->Rvir,tmp_fofnode->Rvir_anyl); */
		      
							/* 	      /\* fprintf(fp,"%10d   %14"STR_FMT"  %20"STR_FMT"  %20"STR_FMT"    %16.6g  %16.6g  %20.6g  %16.6g  %12.4g  %12.4g  %10d  %10.4lf  %20.4g  %20.4g %20.4g \n", *\/ */
							/* 	      /\* 	      subnode->snapshot,tmp_fofnode->haloid,subnode->ContainerHalo->haloid,subhaloid, *\/ */
							/* 	      /\* 	      tmp_fofnode->Mtot,tmp_fofnode->BoundFofMtot,subnode->ContainerHalo->Mtot,subnode->Mtot,dist,vel,tag,destructionz, *\/ */
							/* 	      /\* 	      tmp_fofnode->Rhalf,tmp_fofnode->Rvir,tmp_fofnode->Rvir_anyl); *\/ */

							/* 	      subnode = subnode->Parent; */
							/* 	      tmp_fofnode=tmp_fofnode->Parent; */

							/* 	      if(subnode != NULL && tmp_fofnode!=NULL && subnode->snapshot != tmp_fofnode->snapshot) */
							/* 		ensure_same_snapshot_for_halos(&subnode,&tmp_fofnode); */
		      
							/* 	    } */
							/* 	  } */
							/* 	} */
							/* } */
							thisnode=thisnode->Parent;
						}
					}//if clause -> target halo at target redshift. All of the code lies inside these braces

				}//for loop over all the halos at the (current) target snapshot
				fclose(fp);
      } // for loop over target halo mass
    }// This redshift is actually a valid one and covered in the redshift-range of the simulation
  }//for loop over redshift



  


  
}


void print_header_for_fma(FILE *fp)
{
  fprintf(fp,"####################################################################################################################################################\n");
  fprintf(fp,"#  snapshot      redshift          halonum1       halonum2       haloid1         haloid2         boundfofmtot1     boundfofmtot2       finalfofmtot \n");
  fprintf(fp,"#     i             f                 l              l             l                l                 d                 d                  d        \n");
  fprintf(fp,"####################################################################################################################################################\n");
}




void print_header_for_halfmass(FILE *fp)
{
  fprintf(fp,"##################################################################################\n");
  fprintf(fp,"#  Snapshot      bin       nhalos    bin_mass     meanhalfmassz    halfmass_sigma \n");
  fprintf(fp,"#     i           i          l          d              d                 d        \n");
  fprintf(fp,"##################################################################################\n");

}


int find_fof_flyby_future(struct node_data * const n1,struct node_data * const n2)
{
  /* This is somewhat akin to the get_tag future. Except 
     that this will only return:
     0 -> if this pair will have an encounter that will be catalogued by get_tag/maketree
     1 -> if the pair really flyby and will never come together later.
     
     The const declaration for n1 and n2 ensures that I can't modify the pointer in here. Not
     that I should/would. JIC.
  */

  if(n1==NULL || n2==NULL)
	{
	  fprintf(stderr,"Null pointers encountered in finding the fof future function...... exiting .. \n");
	  exit(EXIT_FAILURE);
	}

  int64 id1  = n1->haloid;
  int64 id2  = n2->haloid;

  struct node_data *node1  = n1;
  struct node_data *node2  = n2;

  /* See if one falls into the other in the future */
  while(node1!=NULL && node2!=NULL && node1->haloid == id1 && node2->haloid == id2)
	{
	  if(node1->ContainerHalo==node2 || node2->ContainerHalo==node1 || node1->FofHalo == node2 || node2->FofHalo == node1)
			return 0;

	  node1 = node1->Parent;
	  node2 = node2->Parent;
	}

  if(node1==NULL || node2==NULL)
		return 1;

  if(node1->haloid != id1 || node2->haloid != id2)
		return 1;
  else
		if(node1->snapshot >= (PARAMS.MAX_SNAPSHOT_NUM))
			return 1;

  return 0;
}





void print_header_for_fof_fof_flyby(FILE *fp)
{

  fprintf(fp,"####################################################################################################################################################################################################################\n");
  fprintf(fp,"#   snapshot           fofnum1          fofnum2         m1               m2              rsep             vsep             rvir1           rvir2       impulse       delpotprim      delpotsec           vdotv     \n");
  fprintf(fp,"#      i                 l                l             d                d                f                f                 f               f            d               d               d                d       \n");
  fprintf(fp,"###################################################################################################################################################################################################################\n");

}




void print_header_for_subsub_flyby(FILE *fp)
{

  fprintf(fp,"###############################################################################################################################################################################################\n");
  fprintf(fp,"# snapshot         groupnum1          groupnum2           m1             m2           rsep           vsep       containernum1    containernum2    impulse       delpotprim      delpotsec  tag \n");
  fprintf(fp,"#    i                 l                  l               d               d             f              f            l                 l             d               d              d        i  \n");
  fprintf(fp,"###############################################################################################################################################################################################\n");

}


int check_for_flyby(struct node_data *first,struct node_data *second,float *rsep,float *vsep)
{
  float v;
  float pre_fac = 3.0;
  *rsep=get_separation_between_centres(first,second);
  
  float r_thresh = 0.0;
  if((first->isFof == 1 && second->ContainerHalo == first) || (second->isFof==1 && first->ContainerHalo==second))
		r_thresh = pre_fac * PARAMS.PhysicalLinkLength;
  else
	{
	  if(first->isFof == 0 && second->isFof == 0 && second->ContainerHalo != first && first->ContainerHalo != second)
			r_thresh = pre_fac * ((first->Rvir > second->Rvir) ? first->Rvir:second->Rvir );
	  else
		{
		  if(second->ContainerHalo == first || first->ContainerHalo == second)
				r_thresh = pre_fac * PARAMS.PhysicalLinkLength;

		}
	}
  

  if(*rsep <= r_thresh)
	{
	  v = 0.0;
	  for(int i=0;i<3;i++)
			v += pow( (first->meanvel[i] - second->meanvel[i]),2.0);

	  v = sqrt(v);
	  *vsep = v;
	  return 1;
	}
  

  return 0;

}
