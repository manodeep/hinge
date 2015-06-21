/*

Log started: July 2009: MS. 

07/01/2009:   started writing the codes. Complete with a Makefile.

07/19/2009:   loadgroups now written and compiles. Needs to be cross-matched
              with known-good IDL routines to make sure that the data is read 
			  in correctly. 

08/24/2009:   No more parent checking. Just read in the parents_xxx.txt and
              assign them. Added 3 categories for mergers:


08/03/2011:   Started the code for computing vmax

NDissolution -> Different FOF parents between subhalos/Fof halos that now reside
            in one FOF halo e.g., one FOF halo comes in and gets destroyed inside
			a different FOF halo.


NMergers     -> A subhalo/FOF halo comes inside a FOF halo and *RETAINS* its identity.
                This halo will continue to exist for multiple snapshots and get tidally
				stripped until it reaches the point where it disappears.


NDisruptions -> A subhalo disappears inside a FOF halo as a result of tidal stripping etc.
                A merger counted in Nmerger is complete when a corresponding entry is seen
				in NDisruptions.

*/

//check 64 bit
#include "defs.h"
#include "proto.h"
#include "read_param.h"
#include "set_cosmology.h"
#include "loadgroups.h"
#include "loadsnapshot.h"
#include "loadparents.h"
#include "loadfillprogenitors.h"
#include "utils.h"
#include "progressbar.h"
#include "genplotdata.h"

int main(int argc, char **argv)
{
  char outfname[MAXLEN];
  /* FILE *fp=NULL; */
  int64 Ngroups0=0;
  int NUM_SNAPSHOTS;
#if (defined(GET_GROUPVEL) == 0 && defined(GET_MEANVEL) == 0)
  struct particle_data *P=NULL;
#endif
  /* char snapshotname[MAXLEN]; */
  struct group_data *group0=NULL;
  struct node_data  *node=NULL;
  struct parent_data *parent=NULL;
  struct cosmology_data COSMO;
  char answer[MAXLEN];

  /* char str_line[MAXLEN]; */
  /* int line=0; */
  time_t t_codestart,t_codeend,t_sectionstart,t_sectionend,t_bigsectionstart;
  
  t_codestart = time(NULL);

  /* Check compilation options */
#if ((defined(WMAP1)  + defined(WMAP3) + defined(WMAP5))  > 1)
#error  Please select only one of the cosmologies from various WMAP options in the Makefile
#endif


#ifdef BIGSIM
  if(sizeof(size_t) != 8)
    { 
      fprintf(stderr,"Error: Code needs to be compiled in 64 bit mode \n");
      fprintf(stderr,"Please add -m64 to the options in the Makefile ..(or find a 64 bit compiler)\n");
      fprintf(stderr,"Exiting..\n");
      exit(EXIT_FAILURE);
    }
#endif


  // Check command line for the parameter file name.
  if (argc !=2)
    { 
      fprintf(stderr,"Usage: %s <parameterfile>\n", argv[0]);
      fprintf(stderr,"Code was compiled with \n\n");
      print_makefile_options();
      exit(EXIT_FAILURE);
    }
  
  my_snprintf(outfname,MAXLEN,"%s",argv[1]);//parameter file

  //read in the parameter file
  fprintf(stderr,"reading parameter file `%s'...",outfname);
  read_params(outfname,&PARAMS);
  fprintf(stderr,"..done\n");

  fprintf(stderr,"sanity checking param values ...");
  sanity_check_params(&PARAMS);
  fprintf(stderr,"..done\n");

  //fill in the config parameters
  fill_config_params(&PARAMS);

  //the declarations that depend on MAX_SNAPSHOT_NUM
  //If MAX_SNAPSHOT_NUM is too large, consider mallocing these
  NUM_SNAPSHOTS = PARAMS.MAX_SNAPSHOT_NUM+1;    
  int64 Ngroups[NUM_SNAPSHOTS];//there are groups corresponding to MAX_SNAPSHOT_NUM
  struct parent_data *  allparents[PARAMS.MAX_SNAPSHOT_NUM];//parents_??? files only go up to MAX_SNAPSHOT_NUM-1
  struct node_data *    tree[NUM_SNAPSHOTS];//there are groups in MAX_SNAPSHOT_NUM

  //read in boxsize and massarr from Gadget header
#ifdef SUBFIND
  my_snprintf(snapshotname,MAXLEN,"%s/%s_%03d",PARAMS.SNAPSHOT_DIR,PARAMS.SNAPSHOT_BASE,PARAMS.MAX_SNAPSHOT_NUM);
  struct io_header header = get_gadget_header(snapshotname);
  NUMPART = get_Numpart(&header);
  //copy Boxsize and massarr to the params structure
  PARAMS.BOXSIZE = header.BoxSize;
  for(int i=0;i<6;i++) {
		if(header.npart[i] > 0 && header.mass[i] <= 0.0) { //needed to make Rvir calculation
			fprintf(stderr,"ERROR: Gadget snapshot has individual masses. This code can not handle that yet\n");//loadgroups with fof_only will not be able to handle this.
			exit(EXIT_FAILURE);
		}
		PARAMS.MASSARR[i] = header.mass[i];
	}
#endif

	NUMPART = 1536*(int64_t) 1536 * (int64_t) 1536;
	
  // output the parameter file 
  my_snprintf(outfname,MAXLEN,"%s/complete_mergertree.params",PARAMS.OUTPUT_DIR);
  fprintf(stderr,"output parameter file to `%s'...",outfname);
  output_params(outfname,&PARAMS);
  fprintf(stderr,"..done\n");


  REDSHIFT = my_malloc(sizeof(*REDSHIFT),NUM_SNAPSHOTS);

#ifndef SUSSING_TREES
	my_snprintf(outfname,MAXLEN,"%s/redshift",PARAMS.GROUP_DIR);
#else
	my_snprintf(outfname,MAXLEN,"%s/redshifts.list",PARAMS.GROUP_DIR);
#endif

	{
		fprintf(stderr,"Reading redshifts from file `%s'\n",outfname);
		FILE *fd = my_fopen(outfname,"rt");
		int line=0;
		char buffer[MAXLINESIZE];
		while(line < NUM_SNAPSHOTS) {
			if(fgets(buffer, MAXLINESIZE,fd) !=NULL) {
				int nread = sscanf(buffer," %f ",&REDSHIFT[line]);
				if(nread == 1) {
					fprintf(stderr,"REDSHIFT[%d] = %g \n",line,REDSHIFT[line]);
					line++;
				}
			} else {
				fprintf(stderr,"WARNING: DID not find enough redshifts (expected %d, found %d) in the redshift file `%s'\n",NUM_SNAPSHOTS,line,outfname);
				break;
			}
		}
		fclose(fd);
	}


	set_cosmology(&COSMO);
  PARAMS.COSMO = &COSMO;
  set_simulation_params(&PARAMS);

  for(int isnapshot=PARAMS.MIN_SNAPSHOT_NUM;isnapshot <=PARAMS.MAX_SNAPSHOT_NUM;isnapshot++)
    Ngroups[isnapshot] = 0;


  /* WARNING: the definition of Ngroups changes with compilation option FOF_ONLY */
  for(int isnapshot = PARAMS.MIN_SNAPSHOT_NUM;isnapshot <= PARAMS.MAX_SNAPSHOT_NUM;isnapshot++)	
    {
      t_bigsectionstart = time(NULL);
      fprintf(stderr,"\n\n Now working on snapshot# %4d \n",isnapshot);
			/* read in Ngroups.  */
#ifdef SUBFIND
			my_snprintf(outfname, MAXLEN,"%s/groups_%03d.fofcat", PARAMS.GROUP_DIR,isnapshot);
			Ngroups0 = returnNhalo(outfname);
#endif



#ifdef SUSSING_TREES
#define RETURN_ONLY_FOFS 0
			/* my_snprintf(outfname,MAXLEN,"%s/%s_%03d.z%5.3f.AHF_halos", PARAMS.GROUP_DIR, PARAMS.GROUP_BASE,isnapshot,REDSHIFT[isnapshot]); */
			my_snprintf(outfname,MAXLEN,"%s/%s%05d.z%5.3f.AHF_halos", PARAMS.GROUP_DIR, PARAMS.GROUP_BASE,isnapshot,REDSHIFT[isnapshot]);
			Ngroups0 = returnNhalo_SUSSING(outfname,RETURN_ONLY_FOFS);
#undef RETURN_ONLY_FOFS
#endif

#ifdef BGC2
#define RETURN_ONLY_FOFS 0
			my_snprintf(outfname,MAXLEN,"%s/halos_%03d.0.bgc2",  PARAMS.GROUP_DIR, isnapshot);
			Ngroups0 = returnNhalo_bgc2(outfname,RETURN_ONLY_FOFS);
#undef RETURN_ONLY_FOFS			
#endif
			
      if (Ngroups0 > 0)
	{

          if(isnapshot < PARAMS.MAX_SNAPSHOT_NUM)
            { 
	      parent = (struct parent_data *)my_malloc(sizeof(*parent),Ngroups0);
	      my_snprintf(outfname,MAXLEN,"%s/%s_%03d%s",PARAMS.GROUP_DIR,"parents",isnapshot,".txt");
	      fprintf(stderr,"Reading in %s Ngroups0 = %"STR_FMT"\n",outfname,Ngroups0);
	      t_sectionstart = time(NULL);
	      parent  = loadparents(outfname,parent,Ngroups0);
	      t_sectionend = time(NULL);
	      fprintf(stderr," done ...\n\n");
	      print_time(t_sectionstart,t_sectionend,"loadparents");
	      allparents[isnapshot] = parent;
	    }
	  
	  group0=allocate_group(Ngroups0);
	  fprintf(stderr,"loading group for snapshot # %d with %"STR_FMT" halos ",isnapshot,Ngroups0);
	  group_init(group0,Ngroups0);
	  
	  t_sectionstart = time(NULL);	  
	  loadgroups(isnapshot,group0);
	  /* loadgroups_no_particles(isnapshot,group0); */
	  fprintf(stderr," done ...\n\n");
	  t_sectionend = time(NULL);
	  print_time(t_sectionstart,t_sectionend,"loadgroups");
	  Ngroups[isnapshot] = Ngroups0;
	  
          t_sectionstart = time(NULL);
          my_snprintf(outfname,MAXLEN,"%s/subhalolevel_%03d.txt",PARAMS.GROUP_DIR,isnapshot);
          readsubhalo_hierarchy_levels(outfname,group0);
          fprintf(stderr," done ...\n\n");
          t_sectionend = time(NULL);
          print_time(t_sectionstart,t_sectionend,"subhalo hierarchy levels");

	  /*Read in the actual snapshot */
#if (defined(GET_GROUPVEL) == 0 && defined(GET_MEANVEL) == 0)
	  my_snprintf(snapshotname, MAXLEN,"%s/%s_%03d",PARAMS.SNAPSHOT_DIR,PARAMS.SNAPSHOT_BASE,isnapshot);
	  t_sectionstart = time(NULL);
	  header = get_gadget_header(snapshotname);
	  P = loadsnapshot(snapshotname,&header);
	  t_sectionend = time(NULL);
	  print_time(t_sectionstart,t_sectionend,"loadsnapshot");
	  
	  t_sectionstart = time(NULL);
	  assign_vxcm(group0,Ngroups0,P,REDSHIFT[isnapshot]);
	  t_sectionend = time(NULL);
	  print_time(t_sectionstart,t_sectionend,"assign meanvel from snapshot");
	  my_free((void **) &P);
#endif

	  /* assign the groups for their snapshot number to the correct locations*/
	  node = (struct node_data *) my_malloc(sizeof(*node),Ngroups0);
	  t_sectionstart = time(NULL);
	  assign_node(group0,Ngroups0,parent,node,isnapshot);
	  t_sectionend = time(NULL);
	  print_time(t_sectionstart,t_sectionend,"assign_node");
	  
	  /* Now store the node data in the tree */
	  tree[isnapshot] = node;
	  free_group(group0,Ngroups0);
	  
	}
	  
      my_snprintf(outfname,MAXLEN,"%s %d ","Reading groups and assigning to node for i = ",isnapshot);/*use outfname as temporary message variable. Gets reset..*/
      t_sectionend = time(NULL);
      print_time(t_bigsectionstart,t_sectionend,outfname);
    }
  
  
  /* All the groups have now been assigned in the node/tree structure. Make the parent/child pointer associations */
  t_sectionstart = time(NULL);
  maketree(allparents,Ngroups,tree);
  t_sectionend = time(NULL);
  print_time(t_sectionstart,t_sectionend,"make_fof_tree");
  
  /* Give the halos unique ids starting at redshift 0 and tracing them all the way back. */
  t_sectionstart = time(NULL);
  assign_haloid(tree,Ngroups);
  t_sectionend = time(NULL);
  print_time(t_sectionstart,t_sectionend,"assign_haloid");

#ifndef FOF_ONLY
  my_snprintf(outfname,MAXLEN,"%s/found_progenitors.txt",PARAMS.GROUP_DIR);
  t_sectionstart = time(NULL);
  load_found_progenitors(tree,Ngroups,outfname);
  assign_haloid(tree,Ngroups);
  t_sectionend = time(NULL);
  print_time(t_sectionstart,t_sectionend,"load_found_progenitors");
#endif

  t_sectionstart = time(NULL);
  find_mergers(tree,Ngroups);
  t_sectionend = time(NULL);
  print_time(t_sectionstart,t_sectionend,"find_mergers");


  t_sectionstart = time(NULL);
  find_subsub_mergers(tree,Ngroups);
  t_sectionend = time(NULL);
  print_time(t_sectionstart,t_sectionend,"find subhalo-subhalo mergers");

  /* output the translation table from group number to group id at all snapshots*/
  t_sectionstart = time(NULL);
  writeids(tree,Ngroups);
  t_sectionend = time(NULL);
  print_time(t_sectionstart,t_sectionend,"writeids");


  /* output the children data */
  t_sectionstart = time(NULL);
  output_children(tree,Ngroups);
  t_sectionend = time(NULL);
  print_time(t_sectionstart,t_sectionend,"output_children");

  /* output the new list of parent data */
  t_sectionstart = time(NULL);
  output_parents(tree,Ngroups);
  t_sectionend = time(NULL);
  print_time(t_sectionstart,t_sectionend,"output_parents");

  /* output the merger data */
  t_sectionstart = time(NULL);
  output_mergers(tree,Ngroups);
  t_sectionend = time(NULL);
  print_time(t_sectionstart,t_sectionend,"output_mergers");

  t_sectionstart = time(NULL);
  cumulative_merger_history(tree,Ngroups);
  t_sectionend = time(NULL);
  print_time(t_sectionstart,t_sectionend,"output cumulative mergers for all snapshots");

  t_sectionstart = time(NULL);
  output_interesting_halos(tree,Ngroups);
  t_sectionend = time(NULL);
  print_time(t_sectionstart,t_sectionend,"output haloids for `interesting' halos");		  

  t_sectionstart = time(NULL);
  output_suspicious_halos(tree,Ngroups);
  t_sectionend = time(NULL);
  print_time(t_sectionstart,t_sectionend,"output haloids for `suspicious' halos");		  

  t_sectionstart = time(NULL);
  output_milkyway_halos(tree,Ngroups);
  t_sectionend = time(NULL);
  print_time(t_sectionstart,t_sectionend,"output all MW type halos at z = 0.0");		  

  t_sectionstart = time(NULL);
  output_flyby_futures(tree,Ngroups);
  t_sectionend = time(NULL);
  print_time(t_sectionstart,t_sectionend,"output all flyby futures ");		  

  {
    if(PARAMS.INDIVIDUAL_MERGERS == 1)
      {
	int64 output_haloid;
	fprintf(stderr,"Do you want to output mergers for a specific haloid (y/n) ? \n");
	scanf("%s",answer);
	while(answer[0] == 'y')
	  {
	    fprintf(stderr,"Enter haloid that you want the merger history for \n");
	    fscanf(stdin,"%"STR_FMT,&output_haloid);
	    t_sectionstart = time(NULL);
	    output_mergers_for_haloid(tree,Ngroups,output_haloid,0.1);
	    t_sectionend = time(NULL);
	    print_time(t_sectionstart,t_sectionend,"output_mergers_for_haloid");
	    fprintf(stderr,"Do you want to output mergers for another haloid (y/n) ? \n");
	    scanf("%s",answer);
	  }  
      }
  }

  t_sectionstart = time(NULL);
  output_all_haloids(tree,Ngroups);
  t_sectionend = time(NULL);
  print_time(t_sectionstart,t_sectionend,"output data for haloids");		  

  t_sectionstart = time(NULL);
  output_plot_data(tree,Ngroups);
  t_sectionend = time(NULL);
  print_time(t_sectionstart,t_sectionend,"output data for making plots");		  


  /* free up the memory */
  for(int isnapshot=PARAMS.MIN_SNAPSHOT_NUM;isnapshot<=PARAMS.MAX_SNAPSHOT_NUM;isnapshot++)
    { 
      if (Ngroups[isnapshot] > 0 )
        { 
          if(isnapshot < PARAMS.MAX_SNAPSHOT_NUM)
            free(allparents[isnapshot]);

          free(tree[isnapshot]);
        }
    }

  free(REDSHIFT);
  free(PARAMS.Age);

  fprintf(stderr,"\n\n Done...\n\n");
  
  t_codeend = time(NULL);
  print_time(t_codestart,t_codeend,"Entire code");
  
  exit(EXIT_SUCCESS);
  
}	



void print_makefile_options()
{
  
#ifdef FOF_ONLY
  fprintf(stderr,"The code was compiled to read in FOF groups only\n");
#else
  fprintf(stderr,"The code was compiled to read in Subfind groups \n");
#endif

#ifdef WMAP1
  fprintf(stderr,"The code was compiled for WMAP 1 parameters \n");
#endif


#ifdef WMAP3
  fprintf(stderr,"The code was compiled for WMAP 3 parameters \n");
#endif

#ifdef WMAP5
  fprintf(stderr,"The code was compiled for WMAP 5 parameters \n");
#endif

#ifdef GET_GROUPVEL
  fprintf(stderr,"The code was compiled to load groups with velocity info\n");
#else
  fprintf(stderr,"The code was compiled to *NOT* read in velocities from the group files \n");
#endif


#ifdef TRACK_BH_BY_ID
  fprintf(stderr,"The code was compiled to track black holes by the id \n");
#endif 

#ifdef TRACK_BH_BY_HALO
  fprintf(stderr,"The code was compiled to track black holes by the centres of halos they were assigned to \n");
#endif


#ifdef LONGIDS
  fprintf(stderr,"Assumes that Gadget particle ids are 8 bytes\n");
#else
  fprintf(stderr,"Assumes that Gadget particle ids are 4 bytes\n");
#endif

#ifdef BIGSIM
  fprintf(stderr,"Assumes particle load is larger than INT_MAX (2 billion)\n");
#else
  fprintf(stderr,"Assumes particle load is smaller than INT_MAX (2 billion)\n");
#endif

  fprintf(stderr,"\n");
 
}



void readsubhalo_hierarchy_levels(const char *fname,struct group_data *group0)
{
  FILE *fp = NULL;
  char str_line[MAXLINESIZE];
  const char comment='#';
  int64 i=0,dummy;
  short dummy1;
  
  fp = my_fopen(fname,"r");
  i=0;
  while(1)
	{
	  if(fgets(str_line, MAXLINESIZE,fp)!=NULL)
		{
		  if(str_line[0] !=comment)
			{
			  /*                    sscanf(str_line,"%*"STR_FMT" %"STR_FMT" %hd   %"STR_FMT"  %"STR_FMT"     %hd",  */
			  sscanf(str_line,"%*d  %"STR_FMT" %hd   %"STR_FMT"  %"STR_FMT"     %hd",
					 &dummy,&(group0[i].ParentLevel),&(group0[i].ContainerIndex),&(group0[i].Nsub),&dummy1); /*discarding the first field, Fofid */
			  
			  if(dummy != i)
				{
				  fprintf(stderr,"Error: While reading in file `%s' for subhalo hierarchy level \n ",fname);
				  fprintf(stderr,"expected igroup = %"STR_FMT" instead got %"STR_FMT" ..exiting \n\n",i,dummy);
				  exit(EXIT_FAILURE);
				  }
			  i++;
			}
		}
	  else
		  break;
	}
  fclose(fp);
  
}


void output_parents(struct node_data * tree[],int64 *Ngroups)
{
  FILE *fd=NULL;
  char outfname[MAXLEN];
  struct node_data *BaseNode=NULL,*thisnode=NULL;
  int interrupted;
  
  fprintf(stderr,"\n\nupdating parentlist.txt files\n");
  init_my_progressbar(PARAMS.MAX_SNAPSHOT_NUM-PARAMS.MIN_SNAPSHOT_NUM+1,&interrupted);
  //there will be no parents at the last snapshot
  for(int isnapshot=PARAMS.MIN_SNAPSHOT_NUM;isnapshot<PARAMS.MAX_SNAPSHOT_NUM;isnapshot++) {
    my_progressbar(isnapshot-PARAMS.MIN_SNAPSHOT_NUM,&interrupted);
    if(Ngroups[isnapshot] > 0) {
      BaseNode = tree[isnapshot];
      my_snprintf(outfname,MAXLEN,"%s/parentlist_%03d.txt",PARAMS.OUTPUT_DIR,isnapshot);
      fd = my_fopen(outfname,"w");
      fprintf(fd,"####################################################################################################################\n");
      fprintf(fd,"# Snapshot      Haloid         GroupNum      isFOF        ParentSnap    ParentId        ParentNum      ParentisFOF\n");
      fprintf(fd,"#    i            l               l              i              i              l               l              i     \n");
      fprintf(fd,"####################################################################################################################\n");
      for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup++) {
	thisnode = &BaseNode[igroup];
	if(thisnode->Parent != NULL)
	  fprintf(fd,"%6d  %14"STR_FMT " %14"STR_FMT"  %12d    %12d  %14"STR_FMT"   %14"STR_FMT "  %12d\n",
		  thisnode->snapshot,thisnode->haloid,thisnode->nodeloc,thisnode->isFof,
		  thisnode->Parent->snapshot,thisnode->Parent->haloid,thisnode->Parent->nodeloc,thisnode->Parent->isFof);
	else
	  fprintf(fd,"%6d  %14"STR_FMT " %14"STR_FMT"  %12d    %12d  %14d   %14d  %12d\n",
		  thisnode->snapshot,thisnode->haloid,thisnode->nodeloc,thisnode->isFof,-1,-1,-1,-1);
	
      }
      fclose(fd);
    }
  }
  /* fprintf(stderr,"..done\n"); */
  finish_myprogressbar(&interrupted);
}



void output_all_haloids(struct node_data * tree[],int64 *Ngroups)
{
  FILE *fp=NULL;
  char fname[MAXLEN];
  struct node_data *BaseNode=NULL,*thisnode=NULL,*tmp=NULL,*first_as_subhalo=NULL;
  int64 haloid=0;
  double maxmass=0.0;

  my_snprintf(fname,MAXLEN,"%s/halolifetimes.txt",PARAMS.OUTPUT_DIR);
  fp = my_fopen(fname,"w");

  fprintf(fp,"######################################################################################################################################################################\n");
  fprintf(fp,"#     Haloid      HaloIDAssignedSnap     FormationZ        FormationSnapshot     DestructionZ    DestructionSnap  HaloMass     EarliestSnap_as_Subhalo    MaxFOFMass  \n");
  fprintf(fp,"#        l              i                    f                   i                   f                 i               d                 i                   d        \n");
  fprintf(fp,"######################################################################################################################################################################\n");


  for(int isnapshot=PARAMS.MAX_SNAPSHOT_NUM;isnapshot>=PARAMS.MIN_SNAPSHOT_NUM;isnapshot--) {
    if(Ngroups[isnapshot] > 0) {
      BaseNode = tree[isnapshot];
      for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup++) {
	first_as_subhalo = NULL;
	maxmass = 0.0;
	thisnode=&BaseNode[igroup];
	if(thisnode->haloid == haloid) {
	  tmp=thisnode;
	  while(tmp->BigChild != NULL) {
	    if(tmp->isFof==0)
	      first_as_subhalo = tmp;
	    else
	      if(tmp->Mtot > maxmass)
		maxmass = tmp->Mtot;
	    
	    tmp=tmp->BigChild;
	  }
	  
	  fprintf(fp,"%14"STR_FMT"   %14d     %14.4g      %12d        %14.4g    %12d  %14.4g",
		  thisnode->haloid,isnapshot,thisnode->FormationRedshift,tmp->snapshot,thisnode->DestructionRedshift,
		  thisnode->snapshot,thisnode->Mtot);
	  
	  if(first_as_subhalo != NULL)
	    fprintf(fp,"  %20d",first_as_subhalo->snapshot);
	  else 
	    fprintf(fp,"  %20d",-1);
	  
	  if(maxmass > 0.0)
	    fprintf(fp," %18.4g",maxmass);
	  else
	    fprintf(fp," %18.4g",-1.0);
	  
	  fprintf(fp,"\n");
	  
	  haloid++;
	}
      }
    }
  }
  

  fclose(fp);
}



void assign_vxcm(struct group_data *group,int64 Ngroups,struct particle_data *P,float redshift)
{
  int64 i,j,index;
  double sumx,sumy,sumz,summ;
  double sumvx,sumvy,sumvz;
  float mass = 0.0;
  for(i=0;i<Ngroups;i++) {
    sumx=0.0;
    sumy=0.0;
    sumz=0.0;
    sumvx=0.0;
    sumvy=0.0;
    sumvz=0.0;
    summ=0.0;
    
#ifdef FOF_ONLY
    group[i].Mtot=0.0;
#endif
    for(j=0;j<group[i].N;j++) {
      index = group[i].id[j];
      mass  = P[index].Mass;
      sumvx += (P[index].Vel[0]*mass);
      sumvy += (P[index].Vel[1]*mass);
      sumvz += (P[index].Vel[2]*mass);
      summ += mass;
#ifndef FOF_ONLY
      sumx += (periodic(P[index].Pos[0]-group[i].xcen)*mass); /* computing the centre of mass -> needs to account for box wrapping. */
      sumy += (periodic(P[index].Pos[1]-group[i].ycen)*mass);
      sumz += (periodic(P[index].Pos[2]-group[i].zcen)*mass);
#else 
      sumx += (periodic(P[index].Pos[0]-group[i].xcen)*mass); /* Make sure that group[i].xcen is initiliased by loadgroups under the FOF_ONLY flag. */
      sumy += (periodic(P[index].Pos[1]-group[i].ycen)*mass);
      sumz += (periodic(P[index].Pos[2]-group[i].zcen)*mass);
#endif 
      
#ifdef FOF_ONLY
      group[i].Mtot += mass;
#endif
    }
    
#ifdef FOF_ONLY
    group[i].xcen = periodic_wrap(sumx/summ);
    group[i].ycen = periodic_wrap(sumy/summ);
    group[i].zcen = periodic_wrap(sumz/summ);
#else
    group[i].Xoffset = group[i].xcen - periodic_wrap(sumx/summ);
    group[i].Yoffset = group[i].ycen - periodic_wrap(sumy/summ);
    group[i].Zoffset = group[i].zcen - periodic_wrap(sumz/summ);
#endif
    
    group[i].vxcen = sumvx/summ;
    group[i].vycen = sumvy/summ;
    group[i].vzcen = sumvz/summ;
    
    group[i].vxcen *= sqrt(get_scalefactor(redshift));
    group[i].vycen *= sqrt(get_scalefactor(redshift));
    group[i].vzcen *= sqrt(get_scalefactor(redshift));
    
  }
  
}


void assign_node(struct group_data *group0,int64 Ngroups0,struct parent_data *parent,struct node_data *node,int isnapshot)
{

#ifdef GET_MEANVEL
  FILE *fp=NULL;
  int64 tmp_grpnum=0;
  char str_line[MAXLINESIZE];
  char fname[MAXLEN];
  const char comment='#';

  my_snprintf(fname,MAXLEN,"%s/allgroupvels_%03d.txt",PARAMS.GROUP_DIR,isnapshot);
  fp = my_fopen(fname,"rt");
  while(1) {
    if(fgets(str_line, MAXLINESIZE,fp)!=NULL)
      if(str_line[0] !=comment)
	break;
  }
#endif

  double rhocrit;
  float overdensity; /*Warning: Make sure that the definitions dont confuse between \rho_u and \rho_{crit}. Overdensity of 200 is usually defined wrt \rho_{crit} */
  struct node_data *FOF_Parent=NULL;
  float scale_factor = 1.0/(1.0+REDSHIFT[isnapshot]);
  float root_a = sqrt(scale_factor);
/*   double boundfofmtot=0.0; */
  
  rhocrit = get_rhocrit_at_redshift(REDSHIFT[isnapshot],PARAMS.COSMO);
  rhocrit *= (scale_factor*scale_factor*scale_factor);  /* in co-moving co-ordinates */
  /*   overdensity = get_overdensity(REDSHIFT[isnapshot],PARAMS.COSMO); */
  overdensity = 200.0; /* Using 200 instead of the Bryan Norman overdensity */
  
  for(int64 igroup = 0;igroup<Ngroups0;igroup++) {
/*       if(FOF_Parent!=NULL && group0[igroup].isFof==1) */
/* 		{ */
/* 		  /\* this will assign the sum of all the subhalo masses as the bound fof mtot. however, this */
/* 			 will miss the last fof halo if it has only itself as the bound subhalo -- */
/* 			 need to check for that at the end of the loop *\/ */
/* 		  FOF_Parent->BoundFofMtot = boundfofmtot; */
/* 		} */
      
    if(group0[igroup].isFof == 1) {
      FOF_Parent = &node[igroup];
      /* 		  boundfofmtot = group0[igroup].Mtot; */
      FOF_Parent->BoundFofMtot = group0[igroup].Mtot;
    } else {
      FOF_Parent->BoundFofMtot += group0[igroup].Mtot;
      node[igroup].BoundFofMtot = group0[igroup].Mtot; /* Not really proper -- subhalos have their 'normal' mass as the bound fof mass */
    }
    
    node[igroup].isFof = group0[igroup].isFof;
    /* 	  node[igroup].Nsub     = group0[igroup].Nsub; */
    node[igroup].haloid = -1;

#ifdef SUSSING_TREES
    node[igroup].SUSSING_haloID = group0[igroup].haloID;
#endif
    node[igroup].nodeloc = igroup;
    node[igroup].z = REDSHIFT[isnapshot];
    node[igroup].snapshot = isnapshot;
    node[igroup].xcen = group0[igroup].xcen;
    node[igroup].ycen = group0[igroup].ycen;
    node[igroup].zcen = group0[igroup].zcen;
      
#ifdef GET_GROUPVEL
    node[igroup].vxcen = group0[igroup].vxcen*root_a;
    node[igroup].vycen = group0[igroup].vycen*root_a;
    node[igroup].vzcen = group0[igroup].vzcen*root_a;
#else 
#ifndef GET_MEANVEL
    node[igroup].vxcen = group0[igroup].vxcen;/* group vels have been loaded in from snapshots. Copy them to the nodes. */
    node[igroup].vycen = group0[igroup].vycen;
    node[igroup].vzcen = group0[igroup].vzcen;
#endif
#endif
    
    node[igroup].Mtot = group0[igroup].Mtot;
    node[igroup].InfallMass = group0[igroup].Mtot;/* initialization -- will get the `real' value from maketree later on*/
    node[igroup].InfallSnapshot    = -1;
    node[igroup].Mstar = 0.0;
    
    /* 	  node[igroup].FormationRedshift = REDSHIFT[isnapshot]; */
    node[igroup].FormationRedshift = 0.0;
    node[igroup].RedshiftofLastMerger = -1.0;
    
    node[igroup].TotNmergers = 0;
    node[igroup].Nmergers = 0;
    
    node[igroup].NDisruptions = 0;
    node[igroup].TotNDisruptions = 0;
    
    node[igroup].NDissolutions = 0;
    node[igroup].TotNDissolutions = 0;
    
    node[igroup].NFlybys = 0;
    node[igroup].TotNFlybys = 0;
    
    node[igroup].Xoffset = group0[igroup].Xoffset;
    node[igroup].Yoffset = group0[igroup].Yoffset; 
    node[igroup].Zoffset = group0[igroup].Zoffset; 
    
      /* Set the group level density parameters */
    node[igroup].Rvir_anyl = getrvir_anyl(group0[igroup].Mtot,REDSHIFT[isnapshot],PARAMS.COSMO);
    getrvir_from_overdensity(&group0[igroup],Nbins,rhocrit,overdensity); /* Follows Bullock et al. 2001 MNRAS 321 559*/
    /* compute_vmax(&group0[igroup],REDSHIFT[isnapshot]); */ //Not implemented yet
    
    /*       if(group0[igroup].N > 10000 ) */
    /* 		fprintf(stderr,"DEBUG: After fitting profile\n  Rvir = %lf  Rvir_anyl=%lf  Rhalf = %lf Maxoverdensity=%e  group->N = %"STR_FMT" group->Mtot = %lf rhocrit=%e\n", */
    /* 		group0[igroup].Rvir,node[igroup].Rvir_anyl,group0[igroup].Rhalf,group0[igroup].MaxOverDensity,group0[igroup].N,group0[igroup].Mtot,rhocrit); */
      
    node[igroup].Rhalf     = group0[igroup].Rhalf;
    
    if (group0[igroup].Rvir > 0.0) {
      node[igroup].Rvir      = group0[igroup].Rvir;
      node[igroup].Conc      = group0[igroup].Conc;
    } else {
      node[igroup].Rhalf = node[igroup].Rvir_anyl*0.5;/* hack -> half mass radius = 1/2 Rvir. Better fit in Lokas 2001  */
      node[igroup].Rvir  = node[igroup].Rvir_anyl;
      node[igroup].Conc  = 0.0; /* use concentration and mass to fit concentration  (Bullock 2001) */
    }
    
    node[igroup].MaxOverDensity    = group0[igroup].MaxOverDensity;
    node[igroup].OverDensityThresh = group0[igroup].OverDensityThresh;
    node[igroup].Vmax              = group0[igroup].Vmax;
    node[igroup].RVmax             = group0[igroup].RVmax;
    
    
    //The parent structure will not exist for snapshot PARAMS.MAX_SNAPSHOT_NUM. initialize from
    // the group directory.
    
/* 	  node[igroup].ContainerId     = parent[igroup].containerid; */
/*    node[igroup].ParentLevel     = parent[igroup].parentlevel; */
/*    node[igroup].Nsub            = parent[igroup].nsub; */
/*    node[igroup].ContainerHalo   = &(node[node[igroup].ContainerId]); */
    
    node[igroup].ContainerId     = group0[igroup].ContainerIndex;
    node[igroup].ParentLevel     = group0[igroup].ParentLevel;
    node[igroup].Nsub            = group0[igroup].Nsub;
    node[igroup].ContainerHalo   = &(node[node[igroup].ContainerId]);
    
    node[igroup].VisitedForMassloss = 0;
    node[igroup].VisitedForCumulativeMerger = 0;
    node[igroup].CumulativeNmergers = 0;
    node[igroup].CumulativeNFlybys  = 0;
    node[igroup].FormationRedshift= node[igroup].z;
    node[igroup].DestructionRedshift = node[igroup].z;
    
    if((node[igroup].ContainerHalo)->Nsub == 0) {
      fprintf(stderr,"This should not have happened. The container claims to have no Nsubs. snapshot = %hd igroup = %"STR_FMT" Container id = %"STR_FMT" ..exiting\n",
	      isnapshot,igroup,node[igroup].ContainerId);
      exit(EXIT_FAILURE);
    }
    

    /* group0 is going to be freed once this function returns. No need to fill it up. */
    /* 	  group0[igroup].ContainerId   = parent[igroup].containerid; */
    /* 	  group0[igroup].ParentLevel   = parent[igroup].parentlevel; */
    
    if(isnapshot < PARAMS.MAX_SNAPSHOT_NUM && parent[igroup].parentsnapshot <= PARAMS.MAX_SNAPSHOT_NUM && parent[igroup].parentsnapshot >= PARAMS.MIN_SNAPSHOT_NUM && parent[igroup].parentid >=0 ) {
      /* 			  fprintf(stderr,"isnapshot = %d igroup = %d parent snapshot = %d\n",isnapshot,igroup,parent[igroup].parentsnapshot); */
      node[igroup].ParentID        = parent[igroup].parentid;
      node[igroup].ParentZ         = REDSHIFT[parent[igroup].parentsnapshot];
      node[igroup].ParentSnapshot  = parent[igroup].parentsnapshot;
    } else {
      node[igroup].ParentID        = -1;
      node[igroup].ParentZ         = -1.0;
      node[igroup].ParentSnapshot  = -1;
    }
    node[igroup].FofHalo  = FOF_Parent;
    node[igroup].Sibling  = NULL;
    node[igroup].Parent   = NULL;
    node[igroup].BigChild = NULL;
    node[igroup].Nchild = 0;
    node[igroup].Nbh = 0;
    node[igroup].BH  = NULL;
#ifdef GET_MEANVEL
    sscanf(str_line,"%"STR_FMT" %f  %f  %f\n",&tmp_grpnum,&(node[igroup].meanvel[0]),&(node[igroup].meanvel[1]),&(node[igroup].meanvel[2]));
    if(tmp_grpnum != igroup) {
      fprintf(stderr,"Error: While reading in velocites for the groups\n");
      fprintf(stderr,"expected groupnum = %"STR_FMT"  got groupnum = %"STR_FMT" ..exiting\n",igroup,tmp_grpnum);
      exit(EXIT_FAILURE);
    }
    
    if(fgets(str_line, MAXLINESIZE,fp)==NULL  && igroup != (Ngroups0-1)) {
      fprintf(stderr,"Error (eof): expected more lines to read in from group velocities file `%s' \n",fname);
      fprintf(stderr,"currently at igroup = %"STR_FMT" ..exiting\n",igroup);
      exit(EXIT_FAILURE);
    }
    
      /* MS: 5th Aug, 2010. Took out the sqrt(a) dependence of v */

      for(int k=0;k<3;k++)
	node[igroup].meanvel[k] *= root_a;

#else
#ifdef GET_GROUPVEL
      node[igroup].meanvel[0] = node[igroup].vxcen;
      node[igroup].meanvel[1] = node[igroup].vycen;
      node[igroup].meanvel[2] = node[igroup].vzcen;
#endif
#endif

      node[igroup].Seeded=0; /* This halo has not 'yet' been seeded with a black hole */
      
    }
#ifdef GET_MEANVEL
  fclose(fp);
#endif

}

void writeids(struct node_data * tree[],int64 *Ngroups)
{
  char outfname[MAXLEN];
  /* int PRINTSTEP,SMALLPRINTSTEP; */
  struct node_data *BaseNode;
  struct node_data *thisnode;
  FILE *fp=NULL;
  int interrupted=0;
  fprintf(stderr,"\n\nwriting out table.txt files\n");
  /* PRINTSTEP = (int)floor(0.1*(PARAMS.MAX_SNAPSHOT_NUM-PARAMS.MIN_SNAPSHOT_NUM)); */
  /* SMALLPRINTSTEP = ceil(0.01*(PARAMS.MAX_SNAPSHOT_NUM-PARAMS.MIN_SNAPSHOT_NUM)) > 1 ? ceil(0.01*(PARAMS.MAX_SNAPSHOT_NUM-PARAMS.MIN_SNAPSHOT_NUM)):1; */

  init_my_progressbar(PARAMS.MAX_SNAPSHOT_NUM-PARAMS.MIN_SNAPSHOT_NUM+1,&interrupted);
  for(int isnapshot = PARAMS.MIN_SNAPSHOT_NUM;isnapshot<= PARAMS.MAX_SNAPSHOT_NUM;isnapshot++) {
    /* if(PRINTSTEP > 0 && SMALLPRINTSTEP > 0) { */
    /*   if((isnapshot-PARAMS.MIN_SNAPSHOT_NUM)%PRINTSTEP == 0) */
    /* 	fprintf(stderr,"%d%%",(int)ceil((isnapshot-PARAMS.MIN_SNAPSHOT_NUM)/PRINTSTEP)*10); */
    /*   else */
    /* 	if((isnapshot-PARAMS.MIN_SNAPSHOT_NUM)%SMALLPRINTSTEP==0) */
    /* 	  fprintf(stderr,"."); */
    /* } */

    my_progressbar(isnapshot-PARAMS.MIN_SNAPSHOT_NUM,&interrupted);
    
    BaseNode = tree[isnapshot];
    my_snprintf(outfname,MAXLEN,"%s/table_%03d.txt",PARAMS.OUTPUT_DIR,isnapshot);
    
    /*       fprintf(stderr,"writing to file `%s' \n",outfname); */
    fp = my_fopen(outfname,"w");
    fprintf(fp,"# Mass is in %e Msun. Redshift = %f  \n",ActualMassUnits,REDSHIFT[isnapshot]);
    fprintf(fp,"#  \n");
    fprintf(fp,"# Snapshot        GroupNum        Haloid           Mtot         LastMergerRed     Formationz      Nmergers    TotNmergers     ParentLev       Nsub     FofHaloid    ContainerNum       Rvir          Rvir_anyl         Rhalf       Nflybys       TotNflybys     Mstar    InfallMass       xcen \n");
    fprintf(fp,"#   i                l               l               d                 f               f               l           l              i              l         l             l                d             d                 d           l               l            d          d            f   \n");
    fprintf(fp,"###############################################################################################################################################################################################################################################################################################\n");
    
    for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup++) {
      thisnode = &BaseNode[igroup];
      fprintf(fp,"%4d     %12"STR_FMT "     %12"STR_FMT"    %13.4f    %12.3f      %12.3f    %10"STR_FMT"    %10"STR_FMT"  %10d  %12"STR_FMT"  %12"STR_FMT"    %12"STR_FMT "   %12.4g    %12.4g    %12.4g     %12"STR_FMT"    %12"STR_FMT"   %14.4lf   %14.4lf %14.6g \n",
	      isnapshot,igroup,thisnode->haloid,
	      thisnode->Mtot,thisnode->RedshiftofLastMerger,thisnode->FormationRedshift,thisnode->Nmergers,
	      thisnode->TotNmergers,thisnode->ParentLevel,thisnode->Nsub,thisnode->FofHalo->haloid,thisnode->ContainerId,
	      thisnode->Rvir,thisnode->Rvir_anyl,thisnode->Rhalf,thisnode->NFlybys,thisnode->TotNFlybys,thisnode->Mstar,
							thisnode->InfallMass,thisnode->xcen);
    }
    
    fclose(fp);
  }
  /* fprintf(stderr,"..done\n"); */

  finish_myprogressbar(&interrupted);
}


void output_children(struct node_data * tree[],int64 *Ngroups)
{
  char outfname[300];
  FILE *fd=NULL;
  struct node_data *BaseNode;
  struct node_data *childnode;
  /* int PRINTSTEP,SMALLPRINTSTEP; */
  int interrupted;
  
  fprintf(stderr,"\n\nwriting children.txt files\n");
  /* PRINTSTEP = (int)floor(0.1*(PARAMS.MAX_SNAPSHOT_NUM-PARAMS.MIN_SNAPSHOT_NUM)); */
  /* SMALLPRINTSTEP = ceil(0.01*(PARAMS.MAX_SNAPSHOT_NUM-PARAMS.MIN_SNAPSHOT_NUM)) > 1 ? ceil(0.01*(PARAMS.MAX_SNAPSHOT_NUM-PARAMS.MIN_SNAPSHOT_NUM)):1; */
  init_my_progressbar(PARAMS.MAX_SNAPSHOT_NUM-PARAMS.MIN_SNAPSHOT_NUM,&interrupted);

  for(int isnapshot=PARAMS.MAX_SNAPSHOT_NUM;isnapshot>PARAMS.MIN_SNAPSHOT_NUM;isnapshot--) {
    my_progressbar(PARAMS.MAX_SNAPSHOT_NUM-isnapshot,&interrupted);
    if(Ngroups[isnapshot] > 0) {
      BaseNode = tree[isnapshot];
      my_snprintf(outfname,MAXLEN,"%s/children_%03d.txt",PARAMS.OUTPUT_DIR,isnapshot);
      /* 	  fprintf(stderr,"\n opening children file `%s' \n",outfname); */
      fd = my_fopen(outfname,"w");
      fprintf(fd,"###################################################################################################################################################################\n");
      fprintf(fd,"# Snapshot     GroupNum     HaloId         HierarLevel    Nsub          Mtot            ChildNum         ChildId        ChildSnap         Mchild   ChildHierarlevel\n");
      fprintf(fd,"#    i            l           l                i            l            d                 l                l              i                 d             i       \n");
      fprintf(fd,"###################################################################################################################################################################\n");
      for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup++) {
	childnode = BaseNode[igroup].BigChild;
	if(childnode != NULL) {
	  while(childnode !=NULL) {
	    fprintf(fd,"%4d  %14"STR_FMT"  %10"STR_FMT" %14d   %14"STR_FMT"   %14.4g  %14"STR_FMT"  %14"STR_FMT"  %14d   %16.4g %14d\n",
		    isnapshot,igroup,BaseNode[igroup].haloid,BaseNode[igroup].ParentLevel,BaseNode[igroup].Nsub,
		    BaseNode[igroup].Mtot,childnode->nodeloc,childnode->haloid,childnode->snapshot,childnode->Mtot,childnode->ParentLevel);
	    childnode = childnode->Sibling;
	  }
	} else {
	  fprintf(fd,"%4d  %14"STR_FMT"  %10"STR_FMT" %14d   %14"STR_FMT"   %14.4g  %14d  %14d  %14d   %16.4g  %14d\n",
		  isnapshot,igroup,BaseNode[igroup].haloid,BaseNode[igroup].ParentLevel,BaseNode[igroup].Nsub,
							  BaseNode[igroup].Mtot,-1,-1,-1,-1.0,-1);
	}
      }
      fclose(fd);
    }
  }
  finish_myprogressbar(&interrupted);
}


void output_mergers_for_haloid(struct node_data * tree[],int64 *Ngroups,const int64 haloid,const float mineta)
{
  FILE *fp=NULL,*fp1=NULL;
  char outfname[MAXLEN];
  int isnapshot;
  int64 igroup=0;

  struct node_data *sibling;
  struct node_data *bigchild;
  struct node_data *node=NULL;

  for(isnapshot=PARAMS.MAX_SNAPSHOT_NUM;isnapshot>=PARAMS.MIN_SNAPSHOT_NUM;isnapshot--) {
    if (Ngroups[isnapshot] > 0) {
      igroup = 0;
      node = tree[isnapshot];
      while(node[igroup].haloid != haloid  && igroup < Ngroups[isnapshot]) 
	igroup++;
	  
      
      if (igroup < Ngroups[isnapshot] && node[igroup].haloid == haloid)
	break;
    } else {
      /* Assumes that once Ngroups is zero for some snapshot, it is zero for all previous
	 snapshots. If your simulation does not have this 'usual' behaviour and you really
	 want the history for the disappearing halos, comment this else section. */
      igroup = Ngroups[isnapshot] + 1; /* to satisfy the following if condition*/
      break;
    }
  }


  if(igroup >=Ngroups[isnapshot]) {
    fprintf(stderr,"Could not locate halo with id = %"STR_FMT" ...\n",haloid);
    fprintf(stderr,"returning..\n");
  } else {
    my_snprintf(outfname,MAXLEN,"%s/%s_%"STR_FMT"%s",PARAMS.OUTPUT_DIR,"mergers_for_haloid",haloid,".txt");
    fp = my_fopen(outfname,"w");
    my_snprintf(outfname,MAXLEN,"%s/%s_%"STR_FMT"%s",PARAMS.OUTPUT_DIR,"growth_history_for_haloid",haloid,".txt");
    fp1 = my_fopen(outfname,"w");
    fprintf(fp1,"###############################################################################\n");
    fprintf(fp1,"#   snapshot         groupnum          haloid          mtot         Fofhaloid #\n");
    fprintf(fp1,"###############################################################################\n");
    
    bigchild = &node[igroup];
    while(bigchild->BigChild !=NULL) {
      bigchild = bigchild->BigChild;
      if(bigchild->Parent->Nchild > 1) {
	sibling = bigchild -> Sibling;
	if (sibling->Mtot/bigchild->Mtot > mineta)
	  fprintf(fp," At z = %7.2f  primary mass = %12.4g other mass = %12.4g  ratio = %10.4f\n",bigchild->z, 
		  bigchild->Mtot,sibling->Mtot,(sibling->Mtot/bigchild->Mtot));
      }
      
      if(bigchild->haloid ==0)
	fprintf(stderr,"z = %7.2f  Nmergers = %4"STR_FMT"  TotNmergers = %6"STR_FMT"  NDisruptions = %4"STR_FMT"  TotNdisrupt = %6"STR_FMT"  NDissolution = %4"STR_FMT
		"  TotNDissolutions = %6"STR_FMT"\n",bigchild->z,bigchild->Nmergers,bigchild->TotNmergers,bigchild->NDisruptions,
		bigchild->TotNDisruptions,bigchild->NDissolutions,bigchild->TotNDissolutions);
      
      
      fprintf(fp1,"%10hd    %14"STR_FMT"   %14"STR_FMT"    %12.6g    %14"STR_FMT" \n",bigchild->snapshot,bigchild->nodeloc,bigchild->haloid,bigchild->Mtot,bigchild->FofHalo->haloid);
    }
    fclose(fp);
    fclose(fp1);
  }
}


void output_mergers(struct node_data * tree[],int64 *Ngroups)
{
  int printedbig = 0;
  FILE *Fof=NULL;
  FILE *Fofprop=NULL;

  FILE *Isolated=NULL;
  char outfname[MAXLEN];

  struct node_data *thischild;
  struct node_data *bigchild;
  struct node_data *node;
  struct node_data *thisfofhalo;
  struct node_data *bigchildfofhalo;
  const float MinMergerRatio = 0.1;
  int interrupted;
  
  my_snprintf(outfname,MAXLEN,"%s/%s",PARAMS.OUTPUT_DIR,"mergers.txt");
  Fof = my_fopen(outfname,"w");

  my_snprintf(outfname,MAXLEN,"%s/%s",PARAMS.OUTPUT_DIR,"isolated.txt");
  Isolated = my_fopen(outfname,"w");

  fprintf(Fof,"###########################################################################################################################################################################################\n");
  fprintf(Fof,"# Snap        ParentHaloid     ParentMtot    ParentHLevel      BigChSnap     OtherChSnap        BigChId        OtherChId        BigChHLevel     OtherChHlevel   BigChMtot     OtherChMtot  \n");
  fprintf(Fof,"#  i              l                d             i                i                i               l                l               i                 i             d              d       \n");
  fprintf(Fof,"###########################################################################################################################################################################################\n");


  my_snprintf(outfname,MAXLEN,"%s/%s",PARAMS.OUTPUT_DIR,"haloprops.txt");
  Fofprop = my_fopen(outfname,"w");
  fprintf(Fofprop,"#################################################################################################################\n");
  fprintf(Fofprop,"# HaloSnap      ParentId        HaloId      HaloLoc       HaloMass         xcen            ycen            zcen  \n");
  fprintf(Fofprop,"#    i             l               l           l             d               f              f                f   \n");
  fprintf(Fofprop,"#################################################################################################################\n");

  fprintf(Isolated,"#########################################################\n");
  fprintf(Isolated,"# Snap          Haloid            GroupNum         Mtot  \n");
  fprintf(Isolated,"#   i              l                 l               d   \n");
  fprintf(Isolated,"#########################################################\n");

  fprintf(stderr,"\n\nwriting mergers.txt files\n");
  init_my_progressbar(PARAMS.MAX_SNAPSHOT_NUM-PARAMS.MIN_SNAPSHOT_NUM+1,&interrupted);

  for(int isnapshot=PARAMS.MIN_SNAPSHOT_NUM;isnapshot<=PARAMS.MAX_SNAPSHOT_NUM;isnapshot++) {
    my_progressbar(isnapshot-PARAMS.MIN_SNAPSHOT_NUM,&interrupted);
    if (Ngroups[isnapshot] > 0) {
      /* 		  fprintf(stderr,"\n\n Now writing out all mergers at snapshot %d Ngroups = %"STR_FMT" \n",isnapshot,Ngroups[isnapshot]); */
      node = tree[isnapshot];
      for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup++) {
	if(node[igroup].Nchild > 1) {
	  /* 			  fprintf(stderr,"for igroup = %d nchild = %d\n",igroup,node[igroup].Nchild); */
	  bigchild   = node[igroup].BigChild;
	  thischild  = bigchild;
	  if(bigchild->Sibling !=NULL) {
	    printedbig=0;
	    do
	      {
		thischild = thischild->Sibling;
		
		/* 
		   There may be a situation where the big child is not from the previous snapshot
		   but all the rest of the children are. There needs to be a check for that.
		*/
		
		
		fprintf(Fof,"%4d   %14"STR_FMT"   %14.5g   %10d   %8d      %8d    %14"STR_FMT"   %14"STR_FMT"   %15d   %15d  %14.5g   %14.5g     \n",
			isnapshot,node[igroup].haloid,node[igroup].Mtot,node[igroup].ParentLevel,bigchild->snapshot,thischild->snapshot,bigchild->haloid,
			thischild->haloid,bigchild->ParentLevel,thischild->ParentLevel,bigchild->Mtot,thischild->Mtot);
		
		if(printedbig ==0) {
		  fprintf(Fofprop,"%4d  %14"STR_FMT"   %14"STR_FMT"   %8"STR_FMT"   %14.5lf   %14.5lf  %14.5lf  %14.5lf  \n",
			  bigchild->snapshot,node[igroup].haloid,bigchild->haloid,bigchild->nodeloc,bigchild->Mtot,bigchild->xcen,
			  bigchild->ycen,bigchild->zcen);
		  printedbig = 1;
		}
		
		fprintf(Fofprop,   "%4d  %14"STR_FMT"   %14"STR_FMT"   %8"STR_FMT"   %14.4f  %14.4f  %14.4f  %14.4f  \n",
			thischild->snapshot,node[igroup].haloid,thischild->haloid,thischild->nodeloc,thischild->Mtot,thischild->xcen,
			thischild->ycen,thischild->zcen);
		
	      }while(thischild->Sibling !=NULL);
	  } else {
	    fprintf(stderr,"\nError: The parent claims to have more than one child but I can only find one\n");
	    fprintf(stderr,"Snapshot = %d  group num = %"STR_FMT"  big child = %"STR_FMT" \n",isnapshot,igroup,bigchild->nodeloc);
	    fprintf(stderr,"Aborting...\n");
	    exit(EXIT_FAILURE);
	  }
	  
	} else {
	  if(node[igroup].Nchild == 1) {
	    /* 
	       this is where I write out the haloes that are not doing anything ...as in possible isolated halos. 
	       Not really required since a combination of parents_xxx.txt, mergers.txt and table_xxx.txt will tell
	       you which halos did not have a merger. 
	    */
	    
	    fprintf(Isolated,"%4d   %14"STR_FMT"    %14"STR_FMT"   %14.4g \n",isnapshot,node[igroup].haloid,node[igroup].nodeloc,node[igroup].Mtot);
	    
	  }
	}
      }
    }
  }
  
  fclose(Fof);
  fclose(Fofprop);
  fclose(Isolated);

  finish_myprogressbar(&interrupted);

  /*Here on I just repeat the previous steps -- and produce mergers bigger than some ratio*/


  my_snprintf(outfname,MAXLEN,"%s/%s",PARAMS.OUTPUT_DIR,"bigcompletedmergers.txt");
/*   fprintf(stderr,"About to  open file %s\n",outfname); */
  Fof = my_fopen(outfname,"w");

  fprintf(Fof,"######################################################################################################################################\n");
  fprintf(Fof,"# Snap        ParentHaloid     ParentMtot    BigChSnap   OtherChSnap        BigChId        OtherChId        BigChMtot     OtherChMtot \n");
  fprintf(Fof,"#   i              l               d             i           i                 l              l                d              d       \n");
  fprintf(Fof,"######################################################################################################################################\n");

  fprintf(stderr,"\n\nwriting to '%s'\n",outfname);
  /*   for(isnapshot=PARAMS.MAX_SNAPSHOT_NUM-1;isnapshot >= PARAMS.MIN_SNAPSHOT_NUM;isnapshot--) */
  init_my_progressbar(PARAMS.MAX_SNAPSHOT_NUM-PARAMS.MIN_SNAPSHOT_NUM+1,&interrupted);
  for(int isnapshot=PARAMS.MIN_SNAPSHOT_NUM;isnapshot<=PARAMS.MAX_SNAPSHOT_NUM;isnapshot++) {
    my_progressbar(isnapshot-PARAMS.MIN_SNAPSHOT_NUM,&interrupted);
    if (Ngroups[isnapshot] > 0) {
      /* 		  fprintf(stderr,"\n Now writing out all mergers bigger than %7.2f at snapshot %d \n",MinMergerRatio,isnapshot); */
      node = tree[isnapshot];
      for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup++) {
	if(node[igroup].Nchild > 1) {
	  /* 			  fprintf(stderr,"for igroup = %d nchild = %d\n",igroup,node[igroup].Nchild); */
	  bigchild   = node[igroup].BigChild;
	  thischild  = bigchild;
	  if(bigchild->Sibling !=NULL) {
	    printedbig=0;
	    do
	      {
		thischild = thischild->Sibling;
		if(thischild->Mtot/bigchild->Mtot >= MinMergerRatio) {
		  /* 
		     There may be a situation where the big child is not from the previous snapshot
		     but all the rest of the children are. There needs to be a check for that.
		  */
		  
		  fprintf(Fof,"%4d   %14"STR_FMT"   %14.4f   %8d      %8d    %14"STR_FMT"   %14"STR_FMT"    %14.4f   %14.4f     \n",
			  isnapshot,node[igroup].haloid,node[igroup].Mtot,bigchild->snapshot,thischild->snapshot,bigchild->haloid,
			  thischild->haloid,bigchild->Mtot,thischild->Mtot);
		}
	      }while(thischild->Sibling !=NULL);
	  } else {
	    fprintf(stderr,"\nError: The parent claims to have more than one child but I can only find one\n");
	    fprintf(stderr,"Snapshot = %d  group num = %"STR_FMT"  big child = %"STR_FMT" \n",isnapshot,igroup,bigchild->nodeloc);
	    fprintf(stderr,"Aborting...\n");
	    exit(EXIT_FAILURE);
	  }
	}
      }
    }
  }

  fclose(Fof);
  finish_myprogressbar(&interrupted);
  /* fprintf(stderr,"..done\n"); */

  /*Here on I just repeat the previous steps -- and produce all the mergers that are not complete */

  my_snprintf(outfname,MAXLEN,"%s/%s",PARAMS.OUTPUT_DIR,"mergers_started.txt");
/*   fprintf(stderr,"About to  open file %s\n",outfname); */
  Fof = my_fopen(outfname,"w");

  fprintf(Fof,"##################################################################################################################\n");
  fprintf(Fof,"# Snap        FofHaloid         FofMtot             ChId           ChMtot       PrevFofHaloId     PrevFofHaloSnap \n");
  fprintf(Fof,"#  i              l                d                 l                d              l                  i         \n");
  fprintf(Fof,"##################################################################################################################\n");

  fprintf(stderr,"\n\nwriting to file `%s'\n",outfname);
  init_my_progressbar(PARAMS.MAX_SNAPSHOT_NUM-PARAMS.MIN_SNAPSHOT_NUM+1,&interrupted);
  /*   for(isnapshot=PARAMS.MAX_SNAPSHOT_NUM-1;isnapshot >= PARAMS.MIN_SNAPSHOT_NUM;isnapshot--) */
  for(int isnapshot=PARAMS.MIN_SNAPSHOT_NUM;isnapshot<=PARAMS.MAX_SNAPSHOT_NUM;isnapshot++) {
    my_progressbar(isnapshot-PARAMS.MIN_SNAPSHOT_NUM,&interrupted);
    /*       fprintf(stderr,"\n\n Now writing out all mergers that started at snapshot %d Ngroups = %"STR_FMT" \n",isnapshot,Ngroups[isnapshot]); */
    if (Ngroups[isnapshot] > 0)	{
      node = tree[isnapshot];
      for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup++) {
	if(node[igroup].BigChild !=NULL) {
	  thisfofhalo     = node[igroup].FofHalo;
	  bigchild        = node[igroup].BigChild;
	  bigchildfofhalo = bigchild->FofHalo;
	  if( (thisfofhalo->haloid != bigchildfofhalo->haloid) && (thisfofhalo->haloid != node[igroup].haloid )) {
	    if (node[igroup].Mtot/thisfofhalo->Mtot >= MinMergerRatio) {
	      /* If I can't assign these mergers in make_fof_tree -> assign them here.*/
	      /* 					  node[igroup].Nmergers++; */
	      /* 					  node[igroup].TotNmergers++; */
	      
	      fprintf(Fof,"%4d   %14"STR_FMT"   %14.4f   %14"STR_FMT"    %14.4f    %12"STR_FMT"  %12d\n",
		      isnapshot,thisfofhalo->haloid,thisfofhalo->Mtot,
		      node[igroup].haloid,node[igroup].Mtot,bigchildfofhalo->haloid,bigchildfofhalo->snapshot);
	      
	    }
	  }
	}
      }
    }
  }
  
  fclose(Fof);
  /* fprintf(stderr,"..done\n");   */
  finish_myprogressbar(&interrupted);
  
  fprintf(stderr,"\n\nFinished writing out the mergers and the halo properties\n");
  
}

void cumulative_merger_history(struct node_data * tree[],int64 *Ngroups)
{
  FILE *fp=NULL;
  char outfname[MAXLEN];
  int64 totnmergers;
  int64 totnflybys;
  struct node_data *thisnode;
  struct node_data *BaseNode;
  int interrupted;
  
  fprintf(stderr,"\n\nwriting cumulative_merger.txt files\n");
  /* PRINTSTEP = (int)floor(0.1*(PARAMS.MAX_SNAPSHOT_NUM-PARAMS.MIN_SNAPSHOT_NUM-1)); */
  /* SMALLPRINTSTEP = ceil(0.01*(PARAMS.MAX_SNAPSHOT_NUM-PARAMS.MIN_SNAPSHOT_NUM-1)) > 1 ? ceil(0.01*(PARAMS.MAX_SNAPSHOT_NUM-PARAMS.MIN_SNAPSHOT_NUM-1)):1; */
  init_my_progressbar(PARAMS.MAX_SNAPSHOT_NUM-PARAMS.MIN_SNAPSHOT_NUM+1,&interrupted);
  for(int isnapshot=PARAMS.MAX_SNAPSHOT_NUM;isnapshot>=PARAMS.MIN_SNAPSHOT_NUM;isnapshot--) {
    my_progressbar(PARAMS.MAX_SNAPSHOT_NUM-isnapshot,&interrupted);
    if (Ngroups[isnapshot] > 0) {
      my_snprintf(outfname,MAXLEN,"%s/%s%03d%s",PARAMS.OUTPUT_DIR,"cumulative_mergers_",isnapshot,".txt");

      /* 	  fprintf(stderr,"About to  open file %s\n",outfname); */
      fp = my_fopen(outfname,"w");
      fprintf(fp,"#########################################################################################\n");
      fprintf(fp,"#   Snapshot       HaloId         Mtot      CumulNmergers      CumulNflybys      FormationZ  \n");
      fprintf(fp,"#      i              l             d           l                 l              f       \n");
      fprintf(fp,"#########################################################################################\n");
      
      BaseNode = tree[isnapshot];
      for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup++) {
	thisnode = &BaseNode[igroup];
	if (thisnode->isFof == 1) {
	  totnmergers = 0;
	  totnflybys  = 0;
	  /* 			  totnmergers += (thisnode->Nmergers + thisnode->NDisruptions + thisnode->NDissolutions); */
	  /* 			  while( (thisnode = walk_tree(thisnode)) != NULL) */
	  
	  while(thisnode != NULL) {
	    /* 				  if(thisnode->VisitedForCumulativeMerger == 1 ) */
	    /* 					{ */
	    /* 					  fprintf(stderr,"this should not have happened..i have already visited this node\n"); */
	    /* 					  fprintf(stderr,"thisnode->haloid = %"STR_FMT" thisnode->snapshot = %d  groupnum = %"STR_FMT" \n",thisnode->haloid,thisnode->snapshot, */
	    /* 							  thisnode->nodeloc); */
	    /* 					  exit(EXIT_FAILURE); */
	    /* 					} */
	    /* 				  thisnode->VisitedForCumulativeMerger = 1; */
	    totnmergers += (thisnode->Nmergers + thisnode->NDisruptions + thisnode->NDissolutions);
	    totnflybys  += (thisnode->NFlybys);
	    thisnode = partial_walk_tree(thisnode,&BaseNode[igroup]);
	    /* 				  fprintf(stderr,"thisnode->haloid = %"STR_FMT" thisnode->snapshot = %d\n",thisnode->haloid,thisnode->snapshot); */
	  }
	  thisnode = &BaseNode[igroup];
	  thisnode->CumulativeNmergers = totnmergers;
	  thisnode->CumulativeNFlybys  = totnflybys;
	  fprintf(fp,"%10d  %12"STR_FMT"   %12.4g   %10"STR_FMT"      %10"STR_FMT "%12.4f \n", isnapshot, thisnode->haloid,thisnode->Mtot,totnmergers,
		  totnflybys,thisnode->FormationRedshift);
	  /* 			  fprintf(stderr,"haloid = %"STR_FMT"   totnmergers = %"STR_FMT" \n",thisnode->haloid,totnmergers); */
	}
      }
      fclose(fp);
    }
  }
  finish_myprogressbar(&interrupted);
  /* fprintf(stderr,"..done\n"); */
}




void output_interesting_halos(struct node_data * tree[],int64 *Ngroups)
{
  
  FILE *fp=NULL;
  char outfname[MAXLEN];
  int isnapshot;
  int64 igroup;

  struct node_data *thisnode;
  struct node_data *node;
  struct node_data *container;

  my_snprintf(outfname,MAXLEN,"%s/%s",PARAMS.OUTPUT_DIR,"interesting_halos.txt");
  fp = my_fopen(outfname,"w");
  fprintf(fp,"#Sub-sub halos where the subhalo disappears but the sub-sub doesn't \n\n");
  fprintf(fp,"##############################################################################################################\n");
  fprintf(fp,"#   Snapshot       HaloId         Mtot      ParentLevel     ContainerHaloId     NodeFormZ      ContainerFormZ \n");
  fprintf(fp,"#      i              l            d            i                  l                f                 f       \n");
  fprintf(fp,"##############################################################################################################\n");

  for(isnapshot=PARAMS.MAX_SNAPSHOT_NUM;isnapshot>=PARAMS.MIN_SNAPSHOT_NUM;isnapshot--) {
    if (Ngroups[isnapshot] > 0) { 
      node = tree[isnapshot];
      for(igroup=0;igroup<Ngroups[isnapshot];igroup++) {
	thisnode = &node[igroup];
	if (thisnode->ParentLevel > 2) {
	  container = thisnode->ContainerHalo;
	  if(container->Parent != NULL && thisnode->Parent != NULL) {
	    if(container->Parent->haloid != container->haloid && thisnode->Parent->haloid == thisnode->haloid)
	      fprintf(fp,"%10d  %12"STR_FMT"   %12.4g    %10d    %12"STR_FMT"    %12.4f   %12.4f\n", 
		      isnapshot, thisnode->haloid,thisnode->Mtot,thisnode->ParentLevel,container->haloid,thisnode->FormationRedshift,
		      container->FormationRedshift);
	  }
	}
      }
    }
  }
  
  fclose(fp);

}

void output_suspicious_halos(struct node_data * tree[],int64 *Ngroups)
{
  FILE *fp=NULL;
  FILE *fp1=NULL;
  FILE *fp2=NULL;

  char outfname[MAXLEN];
  int isnapshot;
  int64 ifof;
  
  double Mthresh = 1.0; /* 1d10 Msun */
  float Mratio  = 0.3;

  struct node_data *thisnode,*bigsub;
  struct node_data *node;

  my_snprintf(outfname,MAXLEN,"%s/%s",PARAMS.OUTPUT_DIR,"suspicious_halos_mratio.txt");
  fp = my_fopen(outfname,"w");
  fprintf(fp,"######################################################\n");
  fprintf(fp,"# Mthresh = %12.4g  mass ratio = %f \n",Mthresh,Mratio);
  fprintf(fp,"#######################################################################################################\n");
  fprintf(fp,"# Snapshot       FofHaloId       FofMass          BigSubId      BigSubMass     FofFormZ    BigSubFormZ \n");
  fprintf(fp,"#    i               l              d                l              d              f            f      \n");
  fprintf(fp,"#######################################################################################################\n");

  for(isnapshot=PARAMS.MAX_SNAPSHOT_NUM;isnapshot>=PARAMS.MIN_SNAPSHOT_NUM;isnapshot--) {
    if (Ngroups[isnapshot] > 0) {
      node = tree[isnapshot];
      ifof = 0;
      while(ifof < Ngroups[isnapshot] && (&node[ifof]) != NULL ) {
	thisnode = &node[ifof];
	if(thisnode->Nsub > 1) {
	  bigsub = &node[ifof+1];
	  if( (thisnode->Mtot > Mthresh) &&  bigsub->Mtot/thisnode->Mtot > Mratio )
	    fprintf(fp,"%10d  %12"STR_FMT"   %12.4g     %12"STR_FMT"   %12.4g    %12.4f   %12.4f\n", 
		    isnapshot, thisnode->haloid,thisnode->Mtot,bigsub->haloid,bigsub->Mtot,thisnode->FormationRedshift,
		    bigsub->FormationRedshift);
	}
	ifof += thisnode->Nsub;
      }
    }
  }
  
  fclose(fp);

  
  my_snprintf(outfname,MAXLEN,"%s/%s",PARAMS.OUTPUT_DIR,"suspicious_halos_massgain.txt");
  fp = my_fopen(outfname,"w");
  fprintf(fp,"######################################################\n");
  fprintf(fp,"# Mthresh = %12.4g  mass gain factor  = %f \n",Mthresh,1.0/Mratio);
  fprintf(fp,"######################################################################################################\n");
  fprintf(fp,"# Snapshot    GroupNum      HaloId         Mtot     NextSnapshot      NextGroupNum      MtotNextSnap \n");
  fprintf(fp,"#    i           l             l              d            i                l                  d      \n");
  fprintf(fp,"######################################################################################################\n");


  my_snprintf(outfname,MAXLEN,"%s/%s",PARAMS.OUTPUT_DIR,"suspicious_halos_massloss.txt");
  fp1 = my_fopen(outfname,"w");
  fprintf(fp1,"######################################################\n");
  fprintf(fp1,"# Mthresh = %12.4g  mass loss factor  = %f \n",Mthresh,1.0/Mratio);
  fprintf(fp1,"######################################################################################################\n");
  fprintf(fp1,"# Snapshot    GroupNum       HaloId         Mtot     NextSnapshot     NextGroupNum      MtotNextSnap  \n");
  fprintf(fp1,"#    i           l             l              d            i                l                  d      \n");
  fprintf(fp1,"######################################################################################################\n");

  ifof = 0; /* now ifof will serve as the max. haloid that has been processed */
  for(isnapshot=PARAMS.MAX_SNAPSHOT_NUM;isnapshot>=PARAMS.MIN_SNAPSHOT_NUM;isnapshot--) {
    if(Ngroups[isnapshot] > 0) {
      node = tree[isnapshot];
      for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup++) {
	thisnode = &node[igroup];
	if(thisnode->haloid >= ifof) {
	  while(thisnode != NULL && thisnode->BigChild != NULL) {
	    /* mass gain and the smaller of the two needs to satisfy the mass threshold */
	    if(thisnode->BigChild->Mtot/thisnode->Mtot < Mratio  && thisnode->BigChild->Mtot > Mthresh)
	      fprintf(fp,"%6d   %12"STR_FMT" %12"STR_FMT"  %12.4g   %10d  %12"STR_FMT"  %13.4g \n",
		      thisnode->BigChild->snapshot,thisnode->BigChild->nodeloc,thisnode->BigChild->haloid,
		      thisnode->BigChild->Mtot,thisnode->snapshot,thisnode->nodeloc,thisnode->Mtot);
	    
	    /* mass loss and the smaller of the two needs to satisfy the mass threshold */
	    if(thisnode->BigChild->Mtot/thisnode->Mtot > (1.0/Mratio)  && thisnode->Mtot > Mthresh)
	      fprintf(fp1,"%6d   %12"STR_FMT" %12"STR_FMT"  %12.4g   %10d   %12"STR_FMT" %13.4g \n",
		      thisnode->BigChild->snapshot,thisnode->BigChild->nodeloc,thisnode->BigChild->haloid,
		      thisnode->BigChild->Mtot,thisnode->snapshot,thisnode->nodeloc,thisnode->Mtot);
	    
	    thisnode = thisnode->BigChild;
	  }	
	  
	  ifof++;
	}
      }
    }
  }
  
  fclose(fp);
  fclose(fp1);



  my_snprintf(outfname,MAXLEN,"%s/%s",PARAMS.OUTPUT_DIR,"suspicious_halos_appear.txt");
  fp = my_fopen(outfname,"w");
  fprintf(fp,"######################################################\n");
  fprintf(fp,"# Mthresh = %12.4g  \n",Mthresh);
  fprintf(fp,"###############################################################################################################\n");
  fprintf(fp,"# Snapshot          Groupnum           Haloid         Mtot             Fofhaloid       FofMass         FofNsub \n");
  fprintf(fp,"#    i                 l                  l            d                   l              d               l    \n");
  fprintf(fp,"###############################################################################################################\n");


  my_snprintf(outfname,MAXLEN,"%s/%s",PARAMS.OUTPUT_DIR,"suspicious_halos_dissolve.txt");
  fp1 = my_fopen(outfname,"w");
  fprintf(fp1,"######################################################\n");
  fprintf(fp1,"# Mthresh = %12.4g  \n",Mthresh);
  fprintf(fp1,"#################################################################################################################################\n");
  fprintf(fp1,"# Snapshot          Groupnum          Haloid          Mtot       ParentLevel       ContainerHaloId    Parenthaloid    ParentMass \n");
  fprintf(fp1,"#    i                 l                 l             d           i                    l                   l              d     \n");
  fprintf(fp1,"#################################################################################################################################\n");


  my_snprintf(outfname,MAXLEN,"%s/%s",PARAMS.OUTPUT_DIR,"suspicious_halos_disappear.txt");
  fp2 = my_fopen(outfname,"w");
  fprintf(fp2,"######################################################\n");
  fprintf(fp2,"# Mthresh = %12.4g  \n",Mthresh);
  fprintf(fp2,"############################################################\n");
  fprintf(fp2,"# Snapshot          Groupnum          Haloid          Mtot  \n");
  fprintf(fp2,"#   i                 l                 l              d    \n");
  fprintf(fp2,"############################################################\n");

  for(isnapshot=PARAMS.MAX_SNAPSHOT_NUM;isnapshot>=PARAMS.MIN_SNAPSHOT_NUM;isnapshot--) {
    if (Ngroups[isnapshot] > 0) {
      node = tree[isnapshot];
      for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup++) {
	thisnode = &node[igroup];
	if(thisnode->Mtot > Mthresh && thisnode->BigChild==NULL && thisnode->isFof==0) {
	  fprintf(fp,"%10d   %14"STR_FMT"   %14"STR_FMT"   %12.6g      %14"STR_FMT"   %12.4g    %12"STR_FMT" \n",
		  thisnode->snapshot,thisnode->nodeloc,thisnode->haloid,thisnode->Mtot,thisnode->FofHalo->haloid,
		  thisnode->FofHalo->Mtot,thisnode->FofHalo->Nsub);
	  
	}
			  
	if(thisnode->Parent !=NULL && thisnode->Parent->haloid != thisnode->haloid && thisnode->Mtot > Mthresh) {
	  fprintf(fp1,"%10d   %14"STR_FMT"   %14"STR_FMT"   %12.4g    %12hd    %14"STR_FMT"    %12"STR_FMT"  %12.6g   \n",
		  thisnode->snapshot,thisnode->nodeloc,thisnode->haloid,thisnode->Mtot,thisnode->ParentLevel,
		  thisnode->ContainerHalo->haloid, thisnode->Parent->haloid,thisnode->Parent->Mtot);
	}
	
	
	if(thisnode->Parent==NULL && thisnode->Mtot > Mthresh) {
	  fprintf(fp2,"%10d   %14"STR_FMT"   %14"STR_FMT"   %12.4g   \n",
		  thisnode->snapshot,thisnode->nodeloc,thisnode->haloid,thisnode->Mtot);
	}
      }
    }
  }
  
  fclose(fp);
  fclose(fp1);
  fclose(fp2);
}

void print_flyby_future_header(FILE *fp)
{

  fprintf(fp,"#####################################################################################################################################################################################################\n");
  fprintf(fp,"#    snapshot       redshift              haloid       groupnum      mtot         boundfofmtot       subhaloid        subgroupnum       submtot       subboundfofmtot         tag          comesback \n");
  fprintf(fp,"#       i               f                   l             l           d                 d                l                 l               d                  d                 i             i      \n");
  fprintf(fp,"#####################################################################################################################################################################################################\n");
}


void output_flyby_futures(struct node_data * tree[],int64 *Ngroups)
{
  FILE *fp=NULL,*fp1=NULL;
  const double Mmin = 100.0; /*10^12 Msun/h lower limit for MW*/
  const double Mmax = 200.0; /*2x10^12 Msun/h -> upper limit for MW*/
  const float  MinFormZ = 6.0;
  char outfname[MAXLEN];
  char command[MAXLEN];

/*   const int NUM_HEADER_LINES=4; */
/*   int numheaderlines_printed = 0; */
/*   char infname[MAXLEN]; */
/*   char buffer[MAXLEN]; */

  struct node_data *thisnode=NULL,*mwnode=NULL,*tmpnode=NULL;
  struct node_data *BaseNode=NULL;
  int flag=0,comesback=0;
  int tag,nflybys=0,ncomebacks=0;
  float destructionz,destructionradius;

  my_snprintf(command,MAXLEN,"rm -f %s/MW_flybys_*.txt",PARAMS.OUTPUT_DIR);
  system(command);

  my_snprintf(outfname,MAXLEN,"%s/flybys_thatcomeback.txt",PARAMS.OUTPUT_DIR);
  fp1 = my_fopen(outfname,"w");
  fprintf(fp1,"# Min. halo mass = %lf max. halo mass = %lf min. formation z = %f \n",
		  Mmin,Mmax,MinFormZ);
  print_flyby_future_header(fp1);

  for(int isnapshot=PARAMS.MAX_SNAPSHOT_NUM;isnapshot>PARAMS.MAX_SNAPSHOT_NUM-2;isnapshot--) {
    if (Ngroups[isnapshot] > 0) {
      BaseNode = tree[isnapshot];
      for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup++) {
	mwnode = &BaseNode[igroup];
	if (mwnode->isFof == 1 && mwnode->Mtot >= Mmin && mwnode->Mtot <= Mmax && mwnode->FormationRedshift > MinFormZ) {
	  /* All right, so we have a MW type halo at z ~ 0 */
	  flag=0;
	  while(mwnode !=NULL) {
	    if(mwnode->isFof==0) {
	      flag=1; /* So the halo was a subhalo at some point in the past..reject */
	      break;
	    }
	    mwnode=mwnode->BigChild;
	  }
		  
	  /*Check if this halo has already been processed at some later snapshot */
	  if(flag!=1) {
	    mwnode = &BaseNode[igroup];
	    my_snprintf(outfname,MAXLEN,"%s/%s_%"STR_FMT"%s",PARAMS.OUTPUT_DIR,"MW_flybys",mwnode->haloid,".txt");
	    fp = fopen(outfname,"r");/* so as to not overwrite halos written at the higher snapshot */	
	    if (fp !=NULL) {
	      fclose(fp);
	      flag=1;
	    }
	  }
	  
	  if(flag !=1) {
	    fp = my_fopen(outfname,"w");
	    print_flyby_future_header(fp);
	    while(mwnode != NULL ) {
	      if(mwnode->nodeloc < (Ngroups[mwnode->snapshot]-1))
		thisnode = mwnode + 1; /* goodness, pointer increment..*/
	      else
		break;/* to protect when primary fof does not have any subs and is the last group in the snapshot */
			  
	      while(thisnode != NULL && thisnode->isFof == 0 && thisnode->nodeloc < (Ngroups[thisnode->snapshot]-1) && 
		    thisnode->FofHalo == mwnode ) {
		/* 							  tag = get_tag(thisnode,mwnode,&destructionz,&destructionradius);/\*disregard the last two things*\/ */
		//MS 7th Dec, 2011 - updated to use the new tagging system that accounts for 
		//subhalos missing snapshots
		double destruction_subtmot=-1.0;
		double destruction_parentmtot=-1.0;
		int destruction_snapshot=-1;
		double destruction_mstar=-1.0;
		
		tag = get_tag(thisnode,mwnode,&destructionz,&destructionradius,&destruction_snapshot,&destruction_subtmot,&destruction_parentmtot,&destruction_mstar);
		if (tag == 3 || tag == 5) { /* flyby tags -> should make it into a function that returns the tags from a name */
		    nflybys++;
		    comesback = 0;
		    
		    /* Now lets see all the descendants of thisnode (while conserving haloid) and see if it comes back */
		    tmpnode = thisnode;
		    while(tmpnode != NULL && tmpnode->isFof==0)
		      tmpnode = tmpnode->Parent;/*walk to the future when thisnode becomes a FOF*/
		    
		    while(tmpnode != NULL && tmpnode->haloid == thisnode->haloid) {
		      if(tmpnode->isFof==0 && tmpnode->FofHalo->haloid==mwnode->haloid) {
			comesback++;
			break;
		      }
		      tmpnode = tmpnode->Parent;
		    }
		    
		    /* What happens if I take out the haloid requirement -> thisnode dissolves in something else 
		       which then has a merger with mwnode (NOTE that mwnode will be around till z~0)*/
		    tmpnode = thisnode;
		    while(tmpnode->isFof==0)
		      tmpnode = tmpnode->Parent;/*walk to the future when thisnode becomes a FOF*/
		    
		    while(tmpnode != NULL) {
		      if(tmpnode->isFof==0 && tmpnode->FofHalo->haloid==mwnode->haloid) {
			comesback +=2;
			break;
		      }
		      tmpnode = tmpnode->Parent;
		    }
		    fprintf(fp,"%10d    %14.4f   %14"STR_FMT" %14"STR_FMT" %16.4lf  %16.4lf %14"STR_FMT"  %14"STR_FMT"  %16.4lf   %16.4lf   %12d  %12d \n",
			    mwnode->snapshot,mwnode->z,mwnode->haloid,mwnode->nodeloc,mwnode->Mtot,mwnode->BoundFofMtot,
			    thisnode->haloid,thisnode->nodeloc,thisnode->Mtot,thisnode->BoundFofMtot,
			    tag,comesback);
		    
		    fprintf(fp1,"%10d    %14.4f   %14"STR_FMT" %14"STR_FMT" %16.4lf  %16.4lf %14"STR_FMT"  %14"STR_FMT"  %16.4lf   %16.4lf   %12d  %12d \n",
			    mwnode->snapshot,mwnode->z,mwnode->haloid,mwnode->nodeloc,mwnode->Mtot,mwnode->BoundFofMtot,
			    thisnode->haloid,thisnode->nodeloc,thisnode->Mtot,thisnode->BoundFofMtot,
			    tag,comesback);
		    
		    if(comesback > 0)
		      ncomebacks++;
		    
		}
		thisnode++; /* really?!! */
	      }
	      mwnode = mwnode->BigChild;
	    }
	    fclose(fp);
	  }
	}
      }
    }
  }
  
  fclose(fp1);
  
/*   nwritten = snprintf(command,MAXLEN,"ls -1 %s/MW_flybys*.txt",PARAMS.OUTPUT_DIR); */
/*   check_string_copy(nwritten,MAXLEN);*/
/*   fp = popen(command,"r"); */
/*   if (fp == NULL) { */
/*     fprintf(stderr,"Error: Failed to run command `%s'\n",command ); */
/*     exit(EXIT_FAILURE); */
/*   } */

/*   numheaderlines_printed = 0; */
/*   while (fgets(infname,sizeof(infname)-1, fp) != NULL)  */
/* 	{ */
/* 	  chomp(infname); */
/* 	  fp1 = my_fopen(infname,"r"); */

/* 	  while(fgets(buffer,MAXLEN,fp1) != NULL) */
/* 		{ */
/* 		  if(buffer[0] != '#') */
/* 			fprintf(fp2,"%s",buffer); */
/* 		  else */
/* 			{ */
/* 			  /\* print the header once  *\/ */
/* 			  numheaderlines_printed++; */
/* 			  if(numheaderlines_printed < NUM_HEADER_LINES) */
/* 				fprintf(fp2,"%s",buffer); */
/* 			} */
/* 		} */
/* 	  fclose(fp1); */
/* 	} */


/*   fclose(fp2); /\* close the concatenated file handler*\/ */
/*   pclose(fp);/\* close the process output handler *\/ */

  fprintf(stderr,"\n\n nflybys = %d ncomebacks = %d   \n\n",nflybys,ncomebacks);
  
}

void output_milkyway_halos(struct node_data * tree[],int64 *Ngroups)
{
  FILE *fp=NULL;
  const double Mmin = 100.0; /*10^12 Msun/h lower limit for MW*/
  const double Mmax = 200.0; /*2x10^12 Msun/h -> upper limit for MW*/
  const float  MinFormZ = 6.0;
  char outfname[MAXLEN];
  char command[MAXLEN];
  struct node_data *thisnode=NULL,*mwnode=NULL,*newnode=NULL;
  struct node_data *node=NULL;
  int mwcount=0,fof_mwcount=0;
  int flag=0;
  int64 parentid;
  int parentsnap;

  my_snprintf(command,MAXLEN,"rm -f %s/MW_*.txt",PARAMS.OUTPUT_DIR);
  system(command);

  for(int isnapshot=PARAMS.MAX_SNAPSHOT_NUM;isnapshot>PARAMS.MAX_SNAPSHOT_NUM-1;isnapshot--) {
    if (Ngroups[isnapshot] > 0) {
      node = tree[isnapshot];
      for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup++) {
	thisnode = &node[igroup];
	if (thisnode->ParentLevel == 1 && thisnode->Mtot >= Mmin && thisnode->Mtot <= Mmax && thisnode->FormationRedshift > MinFormZ) {
	  /* All right, so we have a MW type halo at z ~ 0 */
	  flag=0;
	  fof_mwcount++;
	  while(thisnode !=NULL) {
	    if(thisnode->isFof==0) {
	      flag=1;
	      break;
	    }
	    thisnode=thisnode->BigChild;
	  }
	  
	  thisnode = &node[igroup];
	  my_snprintf(outfname,MAXLEN,"%s/%s_%"STR_FMT"%s",PARAMS.OUTPUT_DIR,"MW",thisnode->haloid,".txt");
	  fp = fopen(outfname,"r");/* so as to not overwrite halos written at the higher snapshot */
	  if (fp !=NULL) {
	    fclose(fp);
	    flag=1;
	  }
	  
	  if(flag !=1) {
	    fp = my_fopen(outfname,"w");
	    fprintf(fp,"# MW type halo. Haloid = %"STR_FMT" final snapshot = %d  mtot = %g  formationz = %12.5g\n",
		    thisnode->haloid,isnapshot,thisnode->Mtot,thisnode->FormationRedshift);
	    fprintf(fp,"#############################################################################################################################################\n");
	    fprintf(fp,"# Snapshot       GroupNum            Haloid            Mtot       ParentLevel       Fofhaloid      ParentID      ParentSnap    CumulNmergers \n");
	    fprintf(fp,"#     i              l                  l                d            i                 l              l             i               l       \n");
	    fprintf(fp,"#############################################################################################################################################\n");
	    
	    while( thisnode != NULL ) {
	      parentid = -1;
	      parentsnap = -1;
	      if(thisnode->Parent != NULL) {
		parentid = thisnode->Parent->haloid;
		parentsnap = thisnode->Parent->snapshot;
	      }
	      fprintf(fp," %6d    %14"STR_FMT"    %14"STR_FMT"   %14.4g    %10d     %14"STR_FMT"     %10"STR_FMT"   %10d   %14u  \n",
		      thisnode->snapshot,thisnode->nodeloc,thisnode->haloid,thisnode->Mtot,thisnode->ParentLevel,thisnode->FofHalo->haloid,
		      parentid,parentsnap,thisnode->CumulativeNmergers);
	      thisnode = walk_tree(thisnode);
	    }
	    fclose(fp);
	    thisnode = &node[igroup];
	    mwnode   = &node[igroup];
	    my_snprintf(outfname,MAXLEN,"%s/%s_%"STR_FMT"%s",PARAMS.OUTPUT_DIR,"MW_allsubhalos",thisnode->haloid,".txt");
	    fp = my_fopen(outfname,"w");
	    fprintf(fp,"# MW type halo. Haloid = %"STR_FMT" final snapshot = %d  mtot = %g  formationz = %12.5g\n",
		    thisnode->haloid,isnapshot,thisnode->Mtot,thisnode->FormationRedshift);
	    
	    fprintf(fp,"####################################################################################################################################################\n");
	    fprintf(fp,"# Snapshot       GroupNum            Haloid            Mtot       ParentLevel       FormationZ     CumulNmergers        ParentID         ParentSnap \n");
	    fprintf(fp,"#     i              l                  l                d            i                 f               l                   l                i      \n");					  
	    fprintf(fp,"####################################################################################################################################################\n");
	    
	    while(mwnode !=NULL) {
	      newnode = tree[mwnode->snapshot];
	      for(int64 isub=mwnode->nodeloc;isub < (mwnode->nodeloc+mwnode->Nsub);isub++) {
		thisnode = &newnode[isub];
		parentid = -1;
		parentsnap = -1;
		if(thisnode->Parent != NULL) {
		  parentid = thisnode->Parent->haloid;
		  parentsnap = thisnode->Parent->snapshot;
		}
		
		fprintf(fp," %6d    %14"STR_FMT"    %14"STR_FMT"   %14.4g    %10d   %14.4g   %14u    %14"STR_FMT"     %10d \n",
			thisnode->snapshot,thisnode->nodeloc,thisnode->haloid,thisnode->Mtot,thisnode->ParentLevel,
			/* 									  thisnode->FormationRedshift,thisnode->DestructionRedshift,thisnode->CumulativeNmergers); */
			thisnode->FormationRedshift,thisnode->CumulativeNmergers,parentid,parentsnap);
		if(thisnode->FofHalo != mwnode) {
		  fprintf(stderr,"This should not have happened: Inside milky way type halos routine \n");
		  fprintf(stderr,"The subhalo claims to have some other halo as the FOF ...exiting\n");
		  exit(EXIT_FAILURE);
		}
	      }
	      
	      mwnode = mwnode->BigChild;
	    }
	    fclose(fp);
	    mwcount++;
	  }
	}
      }
    }
  }
  fprintf(stderr,"Number of FOF MW type halos found at z = 0.0 is : %d  all MW (may have been a subhalo at some point) = %d \n",mwcount,fof_mwcount);
}




void output_gill_data(struct node_data * tree[],int64 *Ngroups)
{
  FILE *fp=NULL;
  const double Mmin = 1e4; /*10^14 clusters*/
  char outfname[MAXLEN];
  char message[MAXLEN];
  struct node_data *clusternode=NULL,*thisnode=NULL,*subnode=NULL;
  struct node_data *node=NULL;

  int64 *cluster_haloids=NULL;
  int64 Nclusters=0;

  int64 *subhaloids=NULL;
  int64 Nsubhaloids=0;
  struct node_data **subhalos;
  int *alltags=NULL;


  int snapshot;
  float dmin=PARAMS.BOXSIZE,tmp_dmin;
  int snapshot_for_dmin;
  float rvir1_at_min=0.0,rvir2_at_min=0.0;
  float destructionz,destructionradius;

  //MS - added on 06/22/12
  float dmax=0.0,tmp_dmax=0.0;
  int snapshot_for_dmax;
  float rvir1_at_max=0.0,rvir2_at_max=0.0;

  float Max_Search_Radius=0.0;
  const float Max_Search_Radius_Factor=10.0;
  int tag;

  snapshot = PARAMS.MAX_SNAPSHOT_NUM;
  node = tree[snapshot];
  for(int64 igroup=0;igroup<Ngroups[snapshot];igroup++) {
    clusternode = &node[igroup];
    if(clusternode->isFof==1 && clusternode->BoundFofMtot >= Mmin && halo_is_always_fof(clusternode)==1) {
      //found a new cluster -- now let's locate all the subhalos
      //that are now present + all the halos that have interacted with
      //this cluster (in the past) and are still surviving
      
      Nclusters++;
      snprintf(message,MAXLEN,"%s","cluster_haloids in function output_gill_data");
      cluster_haloids = my_realloc(cluster_haloids,sizeof(int64),Nclusters,message);
      cluster_haloids[Nclusters-1] = clusternode->haloid;
      
      snprintf(outfname,MAXLEN,"%s/gill_data_%"STR_FMT".txt",PARAMS.OUTPUT_DIR,clusternode->haloid);
      fp = my_fopen(outfname,"w");
      
      subhaloids = NULL;
      subhalos   = NULL;
      alltags    = NULL;
      
      Max_Search_Radius = Max_Search_Radius_Factor*clusternode->Rvir;
      for(int64 jgroup=0;jgroup<Ngroups[snapshot];jgroup++) {
	if(jgroup != igroup) {
	  thisnode = &node[jgroup];
	  tmp_dmin = get_separation_between_centres(thisnode,clusternode);
	  //check if the halo is now within the search radius
	  if(tmp_dmin < Max_Search_Radius) {
	    //now check if the halo was ever a subhalo of "clusternode"
	    while(thisnode != NULL) {
	      if(thisnode->FofHalo->haloid == clusternode->haloid)
		break;
	      
	      thisnode = thisnode->BigChild;
	    }
	    
	    if(thisnode != NULL && thisnode->FofHalo->haloid == clusternode->haloid)  {
	      tag = -1;
	      Nsubhaloids++;
	      snprintf(message,MAXLEN,"%s","subhaloids in function output_gill_data");
	      subhaloids = my_realloc(subhaloids,sizeof(int64),Nsubhaloids,message);
	      subhaloids[Nsubhaloids-1] = thisnode->haloid;
	      
	      snprintf(message,MAXLEN,"%s","subhalos in function output_gill_data");
	      subhalos = my_realloc(subhalos,sizeof(struct node_data *),Nsubhaloids,message);
	      subhalos[Nsubhaloids-1] = &node[jgroup];//store the pointer at z=0
	      
	      
	      //now let's go back to when the subhalo first appeared inside this cluster
	      while(thisnode->BigChild != NULL && thisnode->BigChild->FofHalo->haloid==clusternode->haloid)
		thisnode=thisnode->BigChild;
	      
	      if(thisnode != NULL && thisnode->FofHalo->haloid==clusternode->haloid) {
		double destruction_subtmot=-1.0;
		double destruction_parentmtot=-1.0;
		int destruction_snapshot=-1;
		double destruction_mstar=-1.0;
		
		tag = get_tag(thisnode,thisnode->FofHalo,&destructionz,&destructionradius,&destruction_snapshot,&destruction_subtmot,&destruction_parentmtot,&destruction_mstar);
	      }
	      
	      snprintf(message,MAXLEN,"%s","alltags in function output_gill_data");
	      alltags = my_realloc(alltags,sizeof(int),Nsubhaloids,message);
	      alltags[Nsubhaloids-1] = tag;//store the pointer at z=0
	      
	    }
	  }
	}
      }
      //Now all the subhalos that have interacted with this cluster have been identified. 
      fprintf(stderr,"For cluster # %"STR_FMT" there are %"STR_FMT" halos that have interacted with it and survive to z=0\n",Nclusters,Nsubhaloids);
      snprintf(outfname,MAXLEN,"%s/gill_data_%"STR_FMT".txt",PARAMS.OUTPUT_DIR,clusternode->haloid);
      fp = my_fopen(outfname,"w");
      fprintf(fp,"##########################################################################################################################################################################################################################################\n");
      fprintf(fp,"#  ClusterId        SubhaloId         Rsep       Rvir1         Rvir2        MinRsep     Rvir1_at_min    Rvir2_at_min     Snapshot_at_min    InfallMass     FinalMass     tag      dmax   Rvir1_at_max    Rvir2_at_max   Snapshot_at_max   \n");
      fprintf(fp,"#     l                 l               f          f             f             f             f              f                  i                f              f          i         f         f             f                i\n");
      fprintf(fp,"##########################################################################################################################################################################################################################################\n");
      
      for(int64 isub=0;isub<Nsubhaloids;isub++) {
	subnode = subhalos[isub];
	dmin = PARAMS.BOXSIZE;
	
	rvir2_at_min = subnode->Rvir;
	rvir1_at_min = subnode->FofHalo->Rvir;
	snapshot_for_dmin = subnode->snapshot;
	
	while(subnode != NULL && subnode->FofHalo->haloid==clusternode->haloid) {
	  tmp_dmin = get_separation_between_centres(subnode,subnode->FofHalo);
	  dmin = tmp_dmin < dmin ? tmp_dmin:dmin;
	  if(tmp_dmin < dmin) {
	    dmin = tmp_dmin;
	    snapshot_for_dmin = subnode->snapshot;
	    rvir2_at_min = subnode->Rvir;
	    rvir1_at_min = subnode->FofHalo->Rvir;
	  }
	  subnode = subnode->BigChild;
	}
	
	//MS - Added 06/22/12
	//repeat for finding the max separation. but now the sub doesn't have to contained inside the halo. 
	subnode = subhalos[isub];
	while(subnode != NULL) {
	  if(subnode->FofHalo->haloid == clusternode->haloid)
	    break;
	  
	  subnode = subnode->BigChild;
	}
	
	//so now subnode is at the last time it was a subhalo..now let's walk into the future and find 
	//the maximum separation between the sub and the cluster node.
	snapshot_for_dmax = subnode->snapshot;			  
	dmax = 0.0;
	//new scope
	{
	  struct node_data *tmp_fof=subnode->FofHalo;
	  if (tmp_fof->haloid != clusternode->haloid)
	    {
	      fprintf(stderr,"ERROR: (in output_gill_data) Expected to have cluster with haloid = %"STR_FMT" instead I have haloid = %"STR_FMT"..exiting\n",clusternode->haloid,tmp_fof->haloid);
	      exit(EXIT_FAILURE);
	    }
	  //make sure that the subhalo is still itself. (It should be -- since the sub survives to the last snapshot -- where the haloids are assigned.)
	  while(subnode!=NULL && subnode->haloid == subhaloids[isub] && tmp_fof!=NULL && tmp_fof->haloid == clusternode->haloid) {
	    tmp_dmax = get_separation_between_centres(subnode,tmp_fof);
	    if(tmp_dmax > dmax) {
	      dmax = tmp_dmax;
	      snapshot_for_dmax = subnode->snapshot;
	      rvir2_at_max = subnode->Rvir;
	      rvir1_at_max = subnode->FofHalo->Rvir;
	    }
	    subnode = subnode->Parent;
	    tmp_fof = tmp_fof->Parent;
	  }
	}
	
	subnode = subhalos[isub];
	fprintf(fp,"%14"STR_FMT" %14"STR_FMT"  %14.6f     %14.6f  %14.6f  %14.6f   %14.6f   %14.6f  %10d   %14.6f  %14.6f %10d  %14.6f  %14.6f  %14.6f  %10d\n",
		clusternode->haloid,subnode->haloid,get_separation_between_centres(subnode,clusternode),
		subnode->Rvir,clusternode->Rvir, dmin,rvir1_at_min,rvir2_at_min,snapshot_for_dmin,subnode->InfallMass,
		subnode->Mtot,alltags[isub],dmax,rvir1_at_max,rvir2_at_max,snapshot_for_dmax);
      }//tmp_fof will not be visible outside of this closing brace.
      
      
      //we are done outputting data for this cluster
      free(subhaloids);
      free(subhalos);
      free(alltags);
      fclose(fp);
    }
  }
  
  free(cluster_haloids);
}



struct node_data * partial_walk_tree(struct node_data *this,struct node_data *start)
{
  struct node_data *tmp=NULL;
  tmp = this;
  if( tmp->BigChild != NULL)
    return tmp->BigChild;
  else
    {
      if(tmp->Sibling !=NULL)
		return tmp->Sibling;
      else
		{
		  while( (tmp->Sibling==NULL) && (tmp->Parent !=NULL) && (tmp->snapshot < start->snapshot))
			tmp = tmp->Parent;
	  
		  if(tmp->Parent != NULL && tmp->Parent->snapshot >=start->snapshot)
			tmp=NULL;
		  else
			{
			  if(tmp->Sibling !=NULL)
				tmp = tmp->Sibling;
			  else
				tmp = NULL;
			}
		}
    }
  return tmp;
  
}



struct node_data * walk_tree(struct node_data *start)
{
  struct node_data *tmp;
  tmp = start;
  
  if( tmp->BigChild != NULL)
    return tmp->BigChild;
  else
    {
      if(tmp->Sibling !=NULL)
	return tmp->Sibling;
      else
	{
	  while( (tmp->Sibling==NULL) && (tmp->Parent !=NULL))
	    tmp = tmp->Parent;
	  
	  if(tmp->Sibling !=NULL)
	    tmp = tmp->Sibling;
	  else
	    tmp = NULL;
	}
    }
  return tmp;
  
}
