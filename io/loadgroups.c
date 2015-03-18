#include "defs.h"
#include "loadgroups.h"

//To do the re-ordering
#include "sglib.h"

#ifdef BGC2
#include "bgc2.h"
#endif

#include "utils.h"
#include "read_param.h"

#define    FOFMINLEN        (32)
#define    SUBHALOMINLEN    (20)


void reorder_groups_on_array(const int64 Ngroups,struct group_data *group)
{
  // The last two fields, viz., parentgroupforparticle and parentsnapshotforparticle,  are also particle level information but they do not
  // contain any info at this point. But I will swap them as well, since I might (inadvertently) call this function from elsewhere. 

  int64 Nallocated = Ngroups;//some "reasonable" estimate of the max. number of particles in a halo. Tune to suit your needs
  float *sqr_radius=NULL;

	for(int64 igroup=0;igroup<Ngroups;igroup++) {
  
    if(group[igroup].N > Nallocated) {
      free(sqr_radius);
      sqr_radius = NULL;
      Nallocated = group[igroup].N;
    }

    if(sqr_radius == NULL) {
      sqr_radius = my_malloc(sizeof(*sqr_radius), Nallocated);
    }

    for(int64 j=0;j<group[igroup].N;j++) {
      const float dx = periodic(group[igroup].xcen - group[igroup].x[j]);
      const float dy = periodic(group[igroup].ycen - group[igroup].y[j]);
      const float dz = periodic(group[igroup].zcen - group[igroup].z[j]);
      sqr_radius[j] = dx*dx + dy*dy + dz*dz;
    }
    struct group_data *thisgroup = &(group[igroup]);
#ifdef SUSSING_TREES

#define MULTIPLE_ARRAY_EXCHANGER(vartype,a,i,j) { SGLIB_ARRAY_ELEMENTS_EXCHANGER(float,sqr_radius,i,j); \
    SGLIB_ARRAY_ELEMENTS_EXCHANGER(id64, thisgroup->id,i,j);		\
    SGLIB_ARRAY_ELEMENTS_EXCHANGER(int, thisgroup->type,i,j);		\
    SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->x,i,j);		\
    SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->y,i,j);		\
    SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->z,i,j);		\
    SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->ParticleEnergy,i,j);	\
    SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->vx,i,j);		\
    SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->vy,i,j);		\
    SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->vz,i,j);		\
    SGLIB_ARRAY_ELEMENTS_EXCHANGER(int64, thisgroup->parentgroupforparticle,i,j); \
    SGLIB_ARRAY_ELEMENTS_EXCHANGER(int  , thisgroup->parentsnapshotforparticle,i,j); }

#else
#define MULTIPLE_ARRAY_EXCHANGER(vartype,a,i,j) { \
      SGLIB_ARRAY_ELEMENTS_EXCHANGER(float,sqr_radius,i,j);			\
      SGLIB_ARRAY_ELEMENTS_EXCHANGER(id64, thisgroup->id,i,j);		\
      SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->x,i,j);		\
      SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->y,i,j);		\
      SGLIB_ARRAY_ELEMENTS_EXCHANGER(float, thisgroup->z,i,j);		\
      SGLIB_ARRAY_ELEMENTS_EXCHANGER(int, thisgroup->type,i,j);		\
      SGLIB_ARRAY_ELEMENTS_EXCHANGER(int64, thisgroup->parentgroupforparticle,i,j); \
      SGLIB_ARRAY_ELEMENTS_EXCHANGER(int, thisgroup->parentsnapshotforparticle,i,j); }
#endif

    /* SGLIB_ARRAY_QUICK_SORT(float, sqr_radius, thisgroup->N, SGLIB_NUMERIC_COMPARATOR , MULTIPLE_ARRAY_EXCHANGER); */
		/*if the sort is taking too long, replace with HEAP_SORT -> the worst case performance is better*/
		SGLIB_ARRAY_HEAP_SORT(float, sqr_radius, thisgroup->N, SGLIB_NUMERIC_COMPARATOR , MULTIPLE_ARRAY_EXCHANGER);
#undef MULTIPLE_ARRAY_EXCHANGER

  }

	if(sqr_radius != NULL) {
		free(sqr_radius);
	}
}

static inline void mark_particle_for_deletion(struct group_data *this_group, const int64 partloc)
{
	assert(this_group != NULL && "The group pointer can not be NULL while removing a particle");
	id64 *id = this_group->id;

	id[partloc] = -1;
}

static inline void swap_groups(struct group_data *groups, const int64 first, const int64 second)
{
	struct group_data tmp = groups[first];
	groups[first] = groups[second];
	groups[second] = tmp;
}


void remove_particles_from_all_groups(struct group_data *groups, const int Ngroups)
{
	assert(groups != NULL && "The group pointer can not be NULL");
	for(int64 igroup=0;igroup<Ngroups;igroup++) {
		struct group_data *this_group = &(groups[igroup]);
		id64 *id = this_group->id;
		int *type = this_group->type;
		float *x = this_group->x;
		float *y = this_group->y;
		float *z = this_group->z;
		
#ifdef SUSSING_TREES
		float *ParticleEnergy = this_group->ParticleEnergy;
		float *vx = this_group->vx;
		float *vy = this_group->vy;
		float *vz = this_group->vz;
#endif

		int64 *parentgroupforparticle  = this_group->parentgroupforparticle;
		int *parentsnapshotforparticle = this_group->parentsnapshotforparticle;

		int64 num_moved=0;
		const int64 orig_N = this_group->N;
		int64 j=0;
		while(j < this_group->N) {
			if(id[j] != -1) {
				j++;
				continue;
			}
			const int64 source=this_group->N-1;
			id[j] = id[source];
			type[j] = type[source];
			x[j] = x[source];
			y[j] = y[source];
			z[j] = z[source];

#ifdef SUSSING_TREES
			ParticleEnergy[j] = ParticleEnergy[source];
			vx[j] = vx[source];
			vy[j] = vy[source];
			vz[j] = vz[source];
#endif

			parentgroupforparticle[j]    = parentgroupforparticle[source];
			parentsnapshotforparticle[j] = parentsnapshotforparticle[source];

			num_moved++;
			this_group->N--;
		}
		if(num_moved > 0) {
			assert(this_group->N+num_moved == orig_N && "Sum of particles removed and particles left should equal original number of particles") ;
			/* fprintf(stderr,"igroup = %"STR_FMT" number of duplicates removed = %"STR_FMT" original Npart = %"STR_FMT" current N = %"STR_FMT"\n",igroup,num_moved,orig_N,this_group->N); */
		}
	}

	//Now check that the removal worked
	for(int64 igroup=0;igroup<Ngroups;igroup++) {
		struct group_data *this_group = &(groups[igroup]);
		id64 *id = this_group->id;
		for(int64 j=0;j<this_group->N;j++) {
			if(! (*id != -1)) {
				fprintf(stderr,"ERROR: about to crash\n");
				fprintf(stderr,"igroup = %"STR_FMT" j = %"STR_FMT" this_group->N = %"STR_FMT"\n",igroup,j,this_group->N);
			}
			assert(*id != -1 && "id's are correct");
			id++;
		}
	}
}

/* This function should be called after the halos have been
	 sorted such that fofs and subs are contiguous. 
 */
void remove_duplicate_particles(const int64 Ngroups, struct group_data *groups)
{
	int64 FofIndex=0;//actual index and not a fof haloid since that can be different
	const int mem_increase_fac=1.1;
	int64 max_npart = 10000000;
	long *partids          = my_malloc(sizeof(*partids),  max_npart);
	short *halolevel       = my_malloc(sizeof(*halolevel),max_npart);
	int64 *haloindex       = my_malloc(sizeof(*haloindex),   max_npart);
	int64 *partloc_in_halo = my_malloc(sizeof(*partloc_in_halo), max_npart);
	while(FofIndex < Ngroups) {
		int64 location=0;
		/* First take all of the particle ids, the groups they are located in and the total number of particles in the groups */
		for(int64 igroup=FofIndex;igroup< (FofIndex+groups[FofIndex].Nsub);igroup++) {
			if( !(groups[igroup].FOFHalo == FofIndex) ) {
				fprintf(stderr,"ERROR: ABout to crash. igroup = %"STR_FMT " FofIndex = %"STR_FMT" group->FOFHalo = %"STR_FMT"\n",igroup,FofIndex,groups[igroup].FOFHalo);
			}
			assert(groups[igroup].FOFHalo == FofIndex && "The FOFHalo indices should have been setup correctly");

			for(int64 j=0;j<groups[igroup].N;j++) {
				const long this_id = groups[igroup].id[j];
				if(location == max_npart) {
					//reallocate memory
					max_npart *= mem_increase_fac;
					while(max_npart <=location)
						max_npart+=1000;

					partids   = my_realloc(partids,sizeof(*partids),max_npart,"partids");
					halolevel = my_realloc(halolevel,sizeof(*halolevel),max_npart,"halolevel");
					haloindex = my_realloc(haloindex,sizeof(*haloindex),max_npart,"haloindex");
				} 

				/* if(this_id == 2047766) { */
				/* 	fprintf(stderr,"this_id = %ld igroup = %"STR_FMT" j = %"STR_FMT" halolevel = %d group->N = %"STR_FMT" containerindex = %"STR_FMT" groupnum = %"STR_FMT"\n", */
				/* 					this_id,igroup,j,groups[igroup].ParentLevel,groups[igroup].N,groups[igroup].ContainerIndex,groups[igroup].groupnum); */
				/* } */
				
				partids[location]   = this_id;
				//Make sure that Parentlevel is set up before 
				assert(groups[igroup].ParentLevel >= 1 && "ParentLevels need to have been initialized already");
				halolevel[location] = groups[igroup].ParentLevel;
				haloindex[location] = igroup;
				partloc_in_halo[location] = j;
				location++;
			}
		}
		/* fprintf(stderr,"haloindex[%"STR_FMT"-1] = %"STR_FMT"\n",location,haloindex[location-1]); */
		const int64 Npart = location;
		
#define MULTIPLE_ARRAY_EXCHANGER(vartype,a,i,j) { \
			SGLIB_ARRAY_ELEMENTS_EXCHANGER(short,halolevel,i,j);						\
			SGLIB_ARRAY_ELEMENTS_EXCHANGER(long, partids,i,j);							\
			SGLIB_ARRAY_ELEMENTS_EXCHANGER(int64,haloindex,i,j);						\
			SGLIB_ARRAY_ELEMENTS_EXCHANGER(int64, partloc_in_halo,i,j);		}
		
		SGLIB_ARRAY_QUICK_SORT(id64, partids, Npart, SGLIB_NUMERIC_COMPARATOR , MULTIPLE_ARRAY_EXCHANGER);
#undef MULTIPLE_ARRAY_EXCHANGER
		
		int64 ipart=0;
		while(ipart < Npart) {
			const long this_partid = partids[ipart];
			int64 nmatches=0;
			int64 partindex=ipart+1;
			short maxparentlevel=halolevel[ipart];
			/* int64 haloindex_maxparentlevel = haloindex[ipart]; */
			while(partindex < Npart && partids[partindex] == this_partid) {
				if(halolevel[partindex] > maxparentlevel) {
					maxparentlevel = halolevel[partindex];
					/* haloindex_maxparentlevel = haloindex[partindex]; */
				}
				nmatches++;
				partindex++;
			}

			//There are no duplicates -> update ipart and continue with the loop
			if(nmatches==0) {
				ipart++;
				continue;
			}

			/* So we have duplicate particles -> Need to remove this particles from at least one halo */
			/* fprintf(stderr,"ipart = %"STR_FMT" nmatches = %"STR_FMT" this_id = %"STR_ID"\n",ipart,nmatches,this_partid); */
			for(int64 i=ipart;i<=(ipart+nmatches);i++) {
				/* fprintf(stderr,"partids[%"STR_FMT"] = %"STR_ID" this_partid = %"STR_ID"\n",i,partids[i],this_partid); */
				assert(partids[i] == this_partid && "Particle ids must be identical");
				//Check that the same particle does not belong to two different subhalos at the same parentlevel
				//MS 3rd Mar, 2015 - this condition does not hold. Two different subhalos may contain the same particle. 
				/* if(halolevel[i] == maxparentlevel) { */
				/* 	if( !(haloindex_maxparentlevel == haloindex[i])) { */
				/* 		fprintf(stderr,"ERROR: About to crash: partid = %"STR_ID" haloindex[%"STR_FMT"] = %"STR_FMT"  haloindex_maxparentlevel = %"STR_FMT"  halolevel[%"STR_FMT"] = %d maxparentlevel = %d\n", */
				/* 						partids[i], i,haloindex[i],haloindex_maxparentlevel,i,halolevel[i],maxparentlevel); */
				/* 	} */
				/* 	assert(haloindex_maxparentlevel == haloindex[i] && "The same particle can not be in two different halos at the same parent level"); */
				/* } */

				/* Do we need to remove this particle from its group */
				if(halolevel[i] < maxparentlevel) {
					//Remove this particle from its group
					const int64 this_haloindex = haloindex[i];
					assert(this_haloindex >=0 && this_haloindex < (FofIndex + groups[FofIndex].Nsub));
					struct group_data *this_group = &(groups[this_haloindex]);
					const int64 partloc = partloc_in_halo[i];
					assert(partloc >=0 && partloc < this_group->N);
					/* fprintf(stderr,"i = %"STR_FMT" before removing haloindex = %"STR_FMT" group->N = %"STR_FMT"\n",i,this_haloindex,this_group->N); */
					mark_particle_for_deletion(this_group,partloc);
					/* fprintf(stderr,"i = %"STR_FMT " removed particle %"STR_FMT" from haloindex = %"STR_FMT" group->N = %"STR_FMT"\n",i,partloc,this_haloindex,this_group->N); */
				}
			}
			ipart = ipart + nmatches + 1;
		}
		
		assert(groups[FofIndex].Nsub >=1 && "Number of subhalos must be at least 1");
		FofIndex += groups[FofIndex].Nsub;
	}
	assert(FofIndex == Ngroups && "The number of halos processed should be *exactly* equal to the number of halos present");

	free(partids);
	free(halolevel);
	free(haloindex);
	free(partloc_in_halo);

	//The duplciate particles have already been marked -- now remove them
	remove_particles_from_all_groups(groups,Ngroups);
}


#if defined(SUSSING_TREES) || defined(AHF_INPUT)
//Recursively assign the parentlevels
void assign_parentlevel(struct group_data *this_group, const int64 Ngroups, struct group_data *groups)
{
	const int64 hosthaloindex = this_group->ContainerIndex;
	assert(hosthaloindex >= 0 && hosthaloindex < Ngroups && "Host halo index must be within bounds");
	struct group_data *containerhalo = &(groups[hosthaloindex]);

	if(containerhalo->ParentLevel == -1) {
		if(! (containerhalo->ContainerIndex != hosthaloindex)) {
			fprintf(stderr,"ERROR: About to crash -- hosthalo  = %"STR_FMT " host->parentlevel = %d this->Nsub = %"STR_FMT" parentlevel = %d\n",hosthaloindex,containerhalo->ParentLevel,this_group->Nsub,this_group->ParentLevel);
		}
		assert(containerhalo->ContainerIndex != hosthaloindex && "Bug in fixing container indices. If ParentLevel is not set, then halo can not contain itself - infinite recursion will occur");
		assign_parentlevel(containerhalo, Ngroups, groups);
	}

	
	assert(containerhalo->ParentLevel >=1 && "ParentLevels must be at least 1 (for FOFs) and larger for sub/sub-subs etc");
	this_group->ParentLevel = containerhalo->ParentLevel + 1;
}	
#endif


#ifndef SUSSING_TREES
#ifndef ASCII_DATA
#ifndef BGC2
int64 returnNhalo(const char* buf)
{
  FILE* fsub=NULL;
  int64 Nsub;
  size_t one=1;
  fsub = my_fopen(buf,"r");
  fread(&Nsub,sizeof(int64),one,fsub);
  fclose(fsub);
  return Nsub;

}
#endif
#endif
#endif

#ifdef BGC2
int64 returnNhalo(const char *buf, const int fof_only)
{
  
  
}
#endif

#ifdef AHF_INPUT
int64 returnNhalo_AHF(const char* buf,const int fof_only)
{
  long Ngroups=0;
  Ngroups = getnumlines(buf,'#');
  if(fof_only != 1) {
    return Ngroups;
  } else {
    FILE* fp=NULL;
    char buffer[MAXLINESIZE];
    int64 FOFHalo;
    int64 Nhalos=0;
    fp=my_fopen(buf,"r");
    for(int64 i=0;i<Ngroups;i++) {
      if(fgets(buffer,MAXLINESIZE,fp) != NULL) {

				if(buffer[0] != '#') {
					sscanf(buffer,"%*d %"STR_FMT,&FOFHalo);
					//default AHF outputs host halos as -1
					//SUSSING had this set to 0 (this is the only
					//difference between returnNhalo_AHF and returnNhalo_SUSSING
					if(FOFHalo == -1) {
						Nhalos++;
					}
				}
      } else {
				fprintf(stderr,"ERROR: Could not parse the halo file `%s' for %ld values\n",buf,Ngroups);
				exit(EXIT_FAILURE);
      }
    }
    fclose(fp);
    return Nhalos;
  }
  return 0;
}

void loadgroups(int num, struct group_data *group)
{
  
  char catalogue_fname[MAXLEN];
  char partids_fname[MAXLEN];
	char substruct_fname[MAXLEN];
  FILE *fp=NULL;
  
  long Ngroups;
  long haloid;
  char buffer[MAXLINESIZE];
  my_snprintf(catalogue_fname,MAXLEN,"%s/%s_%03d.z%5.3f.AHF_halos"       , PARAMS.GROUP_DIR, PARAMS.GROUP_BASE,num,REDSHIFT[num]);
  my_snprintf(partids_fname  ,MAXLEN,  "%s/%s_%03d.z%5.3f.AHF_particles"   , PARAMS.GROUP_DIR, PARAMS.GROUP_BASE,num,REDSHIFT[num]);
	my_snprintf(substruct_fname,MAXLEN,  "%s/%s_%03d.z%5.3f.AHF_substructure", PARAMS.GROUP_DIR, PARAMS.GROUP_BASE,num,REDSHIFT[num]);

  fp = my_fopen(partids_fname,"r");
  fscanf(fp,"%ld",&Ngroups);
  fclose(fp);

  fp = my_fopen(catalogue_fname,"r");
  int64 i=0;
  while(i<Ngroups) {
    if(fgets(buffer,MAXLINESIZE,fp) != NULL) {
      if(buffer[0] != '#') {
				/* group[i].Mtot = 0.0; */
				group[i].groupnum = i;
				group[i].NParents = 0;
				group[i].Switched = 0;
				
				group[i].N_per_wedge = 0;
				/* initialise the parent finding variables*/
				group[i].ParentId =-1;
				group[i].ParentSnapshot = -1;
				group[i].Ncommon = 0;
				group[i].Rank = 0.0;
				group[i].snapshot = (short) num;
				
				group[i].ParentLevel = -1;
				const int nread  = sscanf(buffer,"%ld %"STR_FMT" %"STR_FMT" %lf %"STR_FMT" %f %f %f %f %f %f" ,
																	&group[i].haloID,
																	&group[i].FOFHalo,
																	&group[i].Nsub,
																	&group[i].Mtot,
																	&group[i].N,
																	&group[i].xcen,
																	&group[i].ycen,
																	&group[i].zcen,
																	&group[i].vxcen,
																	&group[i].vycen,
																	&group[i].vzcen);
				
				/* fprintf(stderr,"group[%"STR_FMT"].haloid = %ld \n",i,group[i].haloID); */
				assert(nread == 11 && "Could not read 11 fields as expected");
				
				if(group[i].FOFHalo == -1) {
					group[i].isFof = 1;
					group[i].FOFHalo = i;
					group[i].ParentLevel = 1;
					group[i].ContainerIndex = i;      
				} else {
					group[i].isFof = 0;
					group[i].ContainerIndex = group[i].FOFHalo;      
				}
				i++;
      }
    }
  }
  fclose(fp);

	/* We have read in the entire halo catalog - however, the subhalos and the hosts are not contiguous yet. First
		 make sure that the halos appear FOF, sub, sub.., next FOF, sub...etc

		 Easiest way to do this is to sort on FofHalo field of the group struct using sglib

		 -MS 23rd Feb, 2015.
		 
	 */

 
#error Not implemented. 
	


	
	
  
  fp = my_fopen(partids_fname,"r");
  fgets(buffer,MAXLINESIZE,fp);
  sscanf(buffer,"%ld",&Ngroups);
  fprintf(stderr,"partids_fname = `%s' Ngroups = %ld \n",partids_fname,Ngroups);
  i=0;
  while(i<Ngroups){
    if(fgets(buffer,MAXLINESIZE,fp) != NULL) {
      if(buffer[0] != '#') {
				/* fprintf(stderr,"group[%"STR_FMT"].haloid = %ld buffer = `%s' ",i,group[i].haloID,buffer); */
				int nread = sscanf(buffer, "%ld %ld",&(group[i].N),&haloid);
				assert(nread == 2);
				assert(haloid == group[i].haloID);
				
				group[i].id             = my_malloc(sizeof(*group[i].id) ,group[i].N);
				group[i].type           = my_malloc(sizeof(*group[i].type) ,group[i].N);
				group[i].ParticleEnergy = my_malloc(sizeof(*group[i].ParticleEnergy),group[i].N);
				group[i].x              = my_malloc(sizeof(*group[i].x),group[i].N);
				group[i].y              = my_malloc(sizeof(*group[i].y),group[i].N);
				group[i].z              = my_malloc(sizeof(*group[i].z),group[i].N);
				group[i].vx             = my_malloc(sizeof(*group[i].vx),group[i].N);
				group[i].vy             = my_malloc(sizeof(*group[i].vy),group[i].N);
				group[i].vz             = my_malloc(sizeof(*group[i].vz),group[i].N);
				
				group[i].parentgroupforparticle    = my_malloc(sizeof(*group[i].parentgroupforparticle),group[i].N);
				group[i].parentsnapshotforparticle = my_malloc(sizeof(*group[i].parentsnapshotforparticle),group[i].N);
				
				for(int64 j=0; j<group[i].N; j++) {
					if(fgets(buffer,MAXLINESIZE,fp) != NULL ) {
						group[i].type[j] = 1;
						sscanf(buffer, "%ld %f %f %f %f %f %f %f",
									 &(group[i].id[j]),
									 &(group[i].ParticleEnergy[j]),       // [Msun/h (km/sec)^2]
									 &(group[i].x[j]),                    // [kpc/h]
									 &(group[i].y[j]),                    // [kpc/h]
									 &(group[i].z[j]),                    // [kpc/h]
									 &(group[i].vx[j]),                   // [km/sec]
									 &(group[i].vy[j]),                   // [km/sec]
									 &(group[i].vz[j])                    // [km/sec]
							);
						group[i].parentgroupforparticle[j]    = -1;
						group[i].parentsnapshotforparticle[j] = -1;
					} else {
						fprintf(stderr,"ERROR: Could not parse file for individual particles \n");
						exit(EXIT_FAILURE);
					}
				}
				i++;
      }
    } else {
      fprintf(stderr,"ERROR: Could not read halo membership info\n");
      exit(EXIT_FAILURE);
    }
  }
  fclose(fp);
  fprintf(stderr,"done\n");
}
#endif //AHF_INPUT




#ifdef SUSSING_TREES
int64 returnNhalo_SUSSING(const char* buf,const int fof_only)
{
  long Ngroups=0;
  Ngroups = getnumlines(buf,'#');
  if(fof_only != 1) {
    return Ngroups;
  } else {
    FILE* fp=NULL;
    char buffer[MAXLINESIZE];
    int64 FOFHalo;
    int64 Nhalos=0;
    fp=my_fopen(buf,"r");
    for(int64 i=0;i<Ngroups;i++) {
      if(fgets(buffer,MAXLINESIZE,fp) != NULL) {
				if(buffer[0] != '#') {
					sscanf(buffer,"%*d %"STR_FMT,&FOFHalo);
					if(FOFHalo == 0) {
						Nhalos++;
					}
				}
				
      } else {
				fprintf(stderr,"ERROR: Could not parse the halo file `%s' for %ld values\n",buf,Ngroups);
				exit(EXIT_FAILURE);
      }
    }
    fclose(fp);
    return Nhalos;
  }

  return 0;
}

void loadgroups(int num,struct group_data *group)
{
  //Loads in the ID data only for the SUSSING mergertrees  
  
  char catalogue_fname[MAXLEN];
  char partids_fname[MAXLEN];
	/* char substruct_fname[MAXLEN]; */
  FILE *fp=NULL;
  
/* #ifndef ASCII_DATA */
/*   my_snprintf(catalogue_fname,MAXLEN,  "%s/%s_%03d/subhalo_tab_%03d.0", PARAMS.GROUP_DIR, PARAMS.GROUP_BASE,num,num); */
/*   my_snprintf(partids_fname,MAXLEN,  "%s/%s_%03d/subhalo_ids_%03d.0", PARAMS.GROUP_DIR, PARAMS.GROUP_BASE,num,num); */
/* #else */
  long Ngroups;
  char buffer[MAXLINESIZE];
	my_snprintf(catalogue_fname,MAXLEN,"%s/%s%05d.z%5.3f.AHF_halos", PARAMS.GROUP_DIR, PARAMS.GROUP_BASE,num,REDSHIFT[num]);
  my_snprintf(partids_fname  ,MAXLEN,"%s/%s%05d.z%5.3f.AHF_particles"   , PARAMS.GROUP_DIR, PARAMS.GROUP_BASE,num,REDSHIFT[num]);
	/* my_snprintf(substruct_fname,MAXLEN,"%s/%s%05d.z%5.3f.AHF_substructure"   , PARAMS.GROUP_DIR, PARAMS.GROUP_BASE,num,REDSHIFT[num]); */
	
  fp = my_fopen(partids_fname,"r");
  fscanf(fp,"%ld",&Ngroups);
  fclose(fp);

  fp = my_fopen(catalogue_fname,"r");
  int64 i=0;
  while(i<Ngroups) {
    if(fgets(buffer,MAXLINESIZE,fp) != NULL) {
      if(buffer[0] != '#') {
				/* group[i].Mtot = 0.0; */
				group[i].groupnum = i;
				group[i].NParents = 0;
				group[i].Switched = 0;
				
				group[i].N_per_wedge = 0;
				/* initialise the parent finding variables*/
				group[i].ParentId =-1;
				group[i].ParentSnapshot = -1;
				group[i].Ncommon = 0;
				group[i].Rank = 0.0;
				group[i].snapshot = (short) num;
				
				group[i].ParentLevel = -1;
				int nread;

				nread  = sscanf(buffer,"%ld %"STR_FMT" %"STR_FMT" %lf %"STR_FMT" %f %f %f %f %f %f" ,
												&group[i].haloID,
												&group[i].FOFHalo,
												&group[i].Nsub,
												&group[i].Mtot,
												&group[i].N,
												&group[i].xcen,
												&group[i].ycen,
												&group[i].zcen,
												&group[i].vxcen,
												&group[i].vycen,
												&group[i].vzcen
					);
				
				/* fprintf(stderr,"group[%"STR_FMT"].haloid = %ld \n",i,group[i].haloID); */
				assert(nread == 11 && "Could not read 11 fields as expected");


				/* Increment the number of subhalos since we need to count the halo
					 itself as one subhalo. This way, even if the halo has no other subhalos,
					 the halo itself is counted. There are numerous loops scattered around
					 in the code that depend on Nsub>=1. 
				 */

				if(group[i].FOFHalo == 0) {
					group[i].isFof = 1;
					group[i].FOFHalo = i;//The Fofhalo is itself
					group[i].ParentLevel = 1;
					group[i].ContainerIndex = i;
					group[i].Nsub++;
				} else {
					group[i].isFof = 0;
					//This was the SUSSING mergertree way of storing the FOF host information
					long hosthaloindex = group[i].FOFHalo - ((long) 1e12* (long) (num) );
					hosthaloindex--;//0-based indexing.
					group[i].ContainerIndex = hosthaloindex;

					if(group[i].Nsub > 0) {
						group[i].Nsub++;
					}
					
					/* This might not be true -- since the host halo might actually not be the FOF we are after
						 and it might be true that the actual FOF halo has not been read-in yet ( I don't think
						 that's the case but I will assume that to be safe).
					 */
					group[i].FOFHalo = -1;
				}
				i++;
      }
    }
  }
  fclose(fp);


  fp = my_fopen(partids_fname,"r");
  fgets(buffer,MAXLINESIZE,fp);
  sscanf(buffer,"%ld",&Ngroups);
  fprintf(stderr,"Reading particles from file  = `%s' Ngroups = %ld ...",partids_fname,Ngroups);
  i=0;
	while(i<Ngroups){
    if(fgets(buffer,MAXLINESIZE,fp) != NULL) {
      if(buffer[0] != '#') {
				int nread;
				/* fprintf(stderr,"group[%"STR_FMT"].haloid = %ld buffer = `%s' ",i,group[i].haloID,buffer); */
				long haloid;
				long numpart_in_halo;
				nread = sscanf(buffer, "%ld %ld",&numpart_in_halo,&haloid);
				assert(nread == 2 && "Parsed both numpart and haloid");

				/* fprintf(stderr,"haloid = %ld \n",haloid); */
				
				//Since I have re-ordered the groups before, this check will not be valid any longer
				/* assert(haloid == group[i].haloID); */

				/* Sussing trees convention for creating the haloid
				 */
				long haloindex = haloid - ((long) 1e12* (long) (num));
				haloindex--;

				/* haloid now contains the index of the group in question. Check that that is true */
				assert(haloindex >=0 && haloindex < Ngroups && "Checking that halo index is within bounds");
				if( ! (haloindex >= 0 && haloindex < Ngroups && haloid == group[haloindex].haloID && numpart_in_halo == group[haloindex].N) ) {
					fprintf(stderr,"ERROR: About to crash. Here's some info about the group/particle that caused the crash\n");
					fprintf(stderr," i = %"STR_FMT" haloid = %ld haloindex = %ld Ngroups = %"STR_FMT" numpart_in_halo = %"STR_FMT"\n",i,haloid,haloindex,Ngroups,numpart_in_halo);
					fprintf(stderr,"group[haloindex].haloid = %"STR_FMT" group.N = %"STR_FMT"\n",group[haloindex].haloID, group[haloindex].N);
				}

				assert(haloindex >= 0 && haloindex < Ngroups && "haloindex must be valid");
				assert(haloindex >= 0 && haloindex < Ngroups && "haloindex must be valid");
				assert(haloid == group[haloindex].haloID && "Checking that haloid matches between particles and halo catalogs");
				assert(numpart_in_halo == group[haloindex].N && "Checking that the number of particles match between particle and halo catalogs");

				/* Now allocate memory and prepare to parse the particles file for the particle positions in the group */
				group[haloindex].id             = my_malloc(sizeof(*group[haloindex].id) ,group[haloindex].N);
				group[haloindex].type           = my_malloc(sizeof(*group[haloindex].type) ,group[haloindex].N);
				group[haloindex].ParticleEnergy = my_malloc(sizeof(*group[haloindex].ParticleEnergy),group[haloindex].N);
				group[haloindex].x              = my_malloc(sizeof(*group[haloindex].x),group[haloindex].N);
				group[haloindex].y              = my_malloc(sizeof(*group[haloindex].y),group[haloindex].N);
				group[haloindex].z              = my_malloc(sizeof(*group[haloindex].z),group[haloindex].N);
				group[haloindex].vx             = my_malloc(sizeof(*group[haloindex].vx),group[haloindex].N);
				group[haloindex].vy             = my_malloc(sizeof(*group[haloindex].vy),group[haloindex].N);
				group[haloindex].vz             = my_malloc(sizeof(*group[haloindex].vz),group[haloindex].N);
				
				group[haloindex].parentgroupforparticle    = my_malloc(sizeof(*group[haloindex].parentgroupforparticle),group[haloindex].N);
				group[haloindex].parentsnapshotforparticle = my_malloc(sizeof(*group[haloindex].parentsnapshotforparticle),group[haloindex].N);
				
				for(int64 j=0; j<group[haloindex].N; j++) {
					if(fgets(buffer,MAXLINESIZE,fp) != NULL ) {
						group[haloindex].type[j] = 1;//This is incorrect because particles might of various types
						nread = sscanf(buffer, "%ld %f %f %f %f %f %f %f",
													 &(group[haloindex].id[j]),
													 &(group[haloindex].ParticleEnergy[j]),       // [Msun/h (km/sec)^2]
													 &(group[haloindex].x[j]),                    // [kpc/h]
													 &(group[haloindex].y[j]),                    // [kpc/h]
													 &(group[haloindex].z[j]),                    // [kpc/h]
													 &(group[haloindex].vx[j]),                   // [km/sec]
													 &(group[haloindex].vy[j]),                   // [km/sec]
													 &(group[haloindex].vz[j])                    // [km/sec]
							);
						assert(nread == 8 && "Parsed 8 fields for the particle properties");
						group[haloindex].parentgroupforparticle[j]    = -1;
						group[haloindex].parentsnapshotforparticle[j] = -1;
					} else {
						fprintf(stderr,"ERROR: Could not parse file for individual particles \n");
						exit(EXIT_FAILURE);
					}
				}
				i++;
      }
    } else {
      fprintf(stderr,"ERROR: Could not read halo membership info\n");
      exit(EXIT_FAILURE);
    }
  }
  fclose(fp);
	fprintf(stderr,"..done\n");

	/* Now assign the parentlevel information */
	for(int64 igroup=0;igroup<Ngroups;igroup++) {
		if(group[igroup].ParentLevel == 1) continue;//FOF halos already have parentlevels assigned -- no need to process them
		/* fprintf(stderr,"Assigning parentlevel for igroup = %"STR_FMT"\n",igroup); */
		struct group_data *this_group = &(group[igroup]);
		assign_parentlevel(this_group, Ngroups, group);
	}

	
/* fprintf(stderr,"BEFORE sorting \n"); */
	/* for(int64 igroup=0;igroup<Ngroups;igroup++) { */
	/* 	fprintf(stderr,"igroup %07"STR_FMT" group->haloID = %10"STR_FMT" group->N = %10"STR_FMT" group->FOFHalo = %10"STR_FMT" group->ContainerIndex = %10"STR_FMT" group->parentlevel = %d group->Nsub=%"STR_FMT"\n" , */
	/* 					igroup,group[igroup].haloID, group[igroup].N,group[igroup].FOFHalo,group[igroup].ContainerIndex,group[igroup].ParentLevel,group[igroup].Nsub); */
	/* } */

	
	//I am creating a linear array to help sort the group struct in the next step
	long *allfofhaloids    = my_malloc(sizeof(*allfofhaloids), Ngroups);
	int64 *containerindices = my_malloc(sizeof(*containerindices), Ngroups);
	int64 *currentindices   = my_malloc(sizeof(*currentindices), Ngroups);
	for(int64 igroup=0;igroup<Ngroups;igroup++) {
		int64 fofindex=group[igroup].ContainerIndex;
		/* fprintf(stderr,"igroup = %"STR_FMT" fofindex = %"STR_FMT" group[igroup].ContainerIndex = %"STR_FMT"\n",igroup,fofindex,group[igroup].ContainerIndex); */
		while(group[fofindex].isFof != 1) {
			/* fprintf(stderr,"[LOOP] igroup = %"STR_FMT" fofindex = %"STR_FMT" group[igroup].ContainerIndex = %"STR_FMT" group[igroup].FOFHalo = %"STR_FMT"\n",igroup,fofindex,group[igroup].ContainerIndex,group[igroup].FOFHalo); */
			fofindex = group[fofindex].ContainerIndex;
			/* fprintf(stderr,"[LOOP] igroup = %"STR_FMT" fofindex = %"STR_FMT" group[igroup].ContainerIndex = %"STR_FMT" group[igroup].FOFHalo = %"STR_FMT"\n",igroup,fofindex,group[igroup].ContainerIndex,group[igroup].FOFHalo); */
		}
		assert(fofindex >= 0 && fofindex < Ngroups && "Fofindex must be valid within [0,Ngroups)");
		allfofhaloids[igroup] = fofindex;
		containerindices[igroup]=group[igroup].ContainerIndex;
		currentindices[igroup] = igroup;
	}

	fprintf(stderr,"Sorting FOFs+Subhalos into contiguous order...");
	/* Now rearrange the groups so that the Fofs and their associated subhalos are contiguous */
#define MULTIPLE_ARRAY_EXCHANGER(vartype,a,i,j) {												\
		SGLIB_ARRAY_ELEMENTS_EXCHANGER(long,allfofhaloids,i,j);							\
		SGLIB_ARRAY_ELEMENTS_EXCHANGER(int64,containerindices,i,j);					\
		SGLIB_ARRAY_ELEMENTS_EXCHANGER(int64,currentindices,i,j);						\
		SGLIB_ARRAY_ELEMENTS_EXCHANGER(struct group_data, group,i,j);				}
	SGLIB_ARRAY_QUICK_SORT(long, allfofhaloids, Ngroups, SGLIB_NUMERIC_COMPARATOR , MULTIPLE_ARRAY_EXCHANGER);
#undef MULTIPLE_ARRAY_EXCHANGER

	/* One more re-ordering to make sure that the FOFs come first and then the subhalos rather than some scrambled order (but all hosthaloindices must be contiguous because of the sorting step */
	{
		int64 start=0;
		while(start < Ngroups) {
			int64 end = start+1;
			int64 fofloc = start;//assume that the fof is the first -- will get updated in the next while loop if that is not true
			const int64 start_fof = allfofhaloids[start];
			/* fprintf(stderr,"start = %"STR_FMT" start_fof = %" STR_FMT" \n",start,start_fof); */
			while(end < Ngroups && allfofhaloids[end] == start_fof) {
				if(group[end].isFof == 1) fofloc = end;
				end++;
			}
			assert(fofloc >=0 && fofloc < Ngroups && "FOF groups must be valid index");
			assert(group[fofloc].isFof == 1 && "Fof location must actually be a fof");
			if(fofloc != start) {
				/* fprintf(stderr,"Before swap: NOT FOF  group[%"STR_FMT"].isFof == %d Nsub = %"STR_FMT". FOF group[%"STR_FMT"].isFof = %d Nsub = %"STR_FMT" end-start = %"STR_FMT"\n",start,group[start].isFof,group[start].Nsub,fofloc,group[fofloc].isFof,group[fofloc].Nsub,end-start); */
				swap_groups(group, start, fofloc);
				int64 tmp = containerindices[start];
				containerindices[start]=containerindices[fofloc];
				containerindices[fofloc]=tmp;

				tmp = currentindices[start];
				currentindices[start]=currentindices[fofloc];
				currentindices[fofloc]=tmp;
				/* fprintf(stderr,"After swap: FOF  group[%"STR_FMT"].isFof == %d Nsub=%"STR_FMT". NOT FOF group[%"STR_FMT"].isFof = %d Nsub = %"STR_FMT" end-start = %"STR_FMT"\n",start,group[start].isFof,group[start].Nsub,fofloc,group[fofloc].isFof,group[fofloc].Nsub,end-start); */
			}
			assert(group[start].isFof == 1 && "FOF must be at the correct index");
			group[start].Nsub = end-start;//Update the number of subhalos to contain all the sub-subs etc
			/* fprintf(stderr,"start = %"PRId64" end = %"PRId64"  Nsub = %"PRId64" start+Nsub = %"PRId64"\n",start,end,group[start].Nsub,start+group[start].Nsub); */
			/* for(i=start;i<end;i++) { */
			/* 	fprintf(stderr,"i=%"PRId64" isfof=%d Nsub=%"PRId64" Parentlevel=%d\n",i,group[i].isFof,group[i].Nsub,group[i].ParentLevel); */
			/* } */
			start = end;
		}
	}
	free(allfofhaloids);
	fprintf(stderr,"...done\n");

	
	fprintf(stderr,"Fixing the FOFHalo and ContainerIndices...");
	//Now fix the fofindex and containerindices 
	for(int64 fofindex=0;fofindex<Ngroups;fofindex+=group[fofindex].Nsub) {
		if( ! (group[fofindex].isFof == 1) ) {
			fprintf(stderr,"ERROR: About to crash - Fof halo doesn't have fof flag set. fofindex = %"STR_FMT" Nsub = %"STR_FMT" haloid = %"STR_FMT" Parentlevel = %d ContainerIndex = %"STR_FMT"\n",
							fofindex,group[fofindex].Nsub,group[fofindex].haloID,group[fofindex].ParentLevel,group[fofindex].ContainerIndex);
		}
		assert(group[fofindex].isFof == 1 && "fofindex must point to a FOF group");
		assert(group[fofindex].Nsub >= 1 && "Number of subhalos must be at least 1");
		group[fofindex].ContainerIndex = fofindex;
		for(int64 igroup=fofindex;igroup< (fofindex+group[fofindex].Nsub);igroup++) {
			group[igroup].FOFHalo = fofindex;
			group[igroup].nodeloc = igroup;
			const int64 orig_containerindex = containerindices[igroup];
			/* const int64 curr_containerindex = currentindices[orig_containerindex]; */
			int64 jgroup;
			for(jgroup=fofindex; jgroup < (fofindex+group[fofindex].Nsub);jgroup++) {
				/* fprintf(stderr,"fofindex = %"PRId64" Nsub = %"PRId64" igroup = %"PRId64" jgroup = %"PRId64" orig_containerindex = %"PRId64" currentindices[jgroup] = %"PRId64" group[jgroup].ParentLev = %d group[igroup].ParentLev = %d\n", */
				/* 				fofindex,group[fofindex].Nsub,igroup,jgroup,orig_containerindex,currentindices[jgroup],group[jgroup].ParentLevel,group[igroup].ParentLevel); */
				if(currentindices[jgroup] == orig_containerindex) {
					break;
				}
			}
			if(igroup == fofindex) {
				if( ! (jgroup == fofindex) ) {
					fprintf(stderr,"ERROR: About to crash. igroup = %"STR_FMT" jgroup = %"STR_FMT" fofindex = %"STR_FMT" groupnum = %"STR_FMT" isfof = %d\n",
									igroup,jgroup,fofindex,group[igroup].groupnum,group[igroup].isFof);
				}
				assert(jgroup == fofindex && "Bug in fixing container indices. FOF must refer to itself as the container");
			}
			assert(currentindices[jgroup] == orig_containerindex && "ERROR: The impossible has happened. Could not find container index for halo");
			group[igroup].ContainerIndex = jgroup;
		}
	}
	free(containerindices);
	free(currentindices);
	fprintf(stderr,"...done\n");
	
	/* fprintf(stderr,"AFTER: sorting\n"); */
	/* for(int64 igroup=0;igroup<Ngroups;igroup++) { */
	/* 	fprintf(stderr,"igroup %07"STR_FMT" group->haloID = %10"STR_FMT" group->N = %10"STR_FMT" group->FOFHalo = %10"STR_FMT" group->ContainerIndex = %10"STR_FMT" group->parentlevel = %d group->Nsub=%"STR_FMT"\n" , */
	/* 					igroup,group[igroup].haloID, group[igroup].N,group[igroup].FOFHalo,group[igroup].ContainerIndex,group[igroup].ParentLevel,group[igroup].Nsub); */
	/* } */


	fprintf(stderr,"Removing duplicate particles...");
	/* All the particles have been read-in. Now remove duplicates */
	remove_duplicate_particles(Ngroups, group);
	fprintf(stderr,"....done\n");

	time_t t_sectionstart=time(NULL);
	fprintf(stderr,"Reordering the particles inside the groups ....");
	/* And re-order particles based on distance from center */
	reorder_groups_on_array(Ngroups,group);
	fprintf(stderr,"....done\n");
	time_t t_sectionend=time(NULL);
	print_time(t_sectionstart,t_sectionend,"Re-ordering particles inside groups");
}  
#endif //SUSSING


#ifdef SUBFIND
void loadgroups(int num,struct group_data* group)
{
  int64 Ngroups = 0;
  int64 NumPart = 0;

  FILE *fpos=NULL;
  FILE *ftype=NULL;
  FILE *fid=NULL;
  FILE *fcat=NULL;
  size_t one=1;

  int64 i,j;
  float rmax,temp_rmax;

  int64 *GroupLen=NULL;
  int64 *GroupOffset=NULL;
  float *xpos=NULL;
  float *ypos=NULL;
  float *zpos=NULL;
  char *type=NULL;
  id64 *id=NULL;
  float tmp_pos[3];

  char catalogue_fname[MAXLEN];
  char particles_fname[MAXLEN];
  char parttypes_fname[MAXLEN];
  char partids_fname[MAXLEN];
  int64 FOF_Parent;

  float *GroupCM[3];
#ifdef GET_GROUPVEL
  float *GroupCMV[3];
#endif

#ifndef FOF_ONLY
  FILE* fprop=NULL;
  FILE* fsubcat=NULL;
  FILE* fsubprop=NULL;
  int64 Nsub = 0;
  char subhalos_fname[MAXLEN];
  char subprop_fname[MAXLEN];
  char fofprop_fname[MAXLEN];
  int64 *GroupSubs=NULL;
  float *GroupMtot=NULL;
  float *GroupMgas=NULL;
  int64 *SubLen=NULL;
  int64 *SubOffset=NULL;
  int64 *SubParent=NULL;
  float *SubMtot=NULL;
  float *SubMgas=NULL;
  int64 sumgroups;
  int64 *tempindex=NULL;
  float *SubCM[3];
#ifdef GET_GROUPVEL
  float *SubCMV[3];
#endif 


#else

  //FOF_ONLY
  double mtot,mgas;
#endif

  my_snprintf(particles_fname,MAXLEN,  "%s/%s_%03d.pos", PARAMS.GROUP_DIR, PARAMS.GROUP_BASE,num);
  my_snprintf(parttypes_fname,MAXLEN,  "%s/%s_%03d.types", PARAMS.GROUP_DIR,PARAMS.GROUP_BASE,num);
  my_snprintf(partids_fname,  MAXLEN,  "%s/%s_%03d.ids", PARAMS.GROUP_DIR,PARAMS.GROUP_BASE,num);
  my_snprintf(catalogue_fname,MAXLEN,  "%s/%s_%03d.fofcat", PARAMS.GROUP_DIR,PARAMS.GROUP_BASE,num);
#ifndef FOF_ONLY
  my_snprintf(fofprop_fname,  MAXLEN,  "%s/%s_%03d.fofprop", PARAMS.GROUP_DIR,PARAMS.GROUP_BASE,num);
  my_snprintf(subhalos_fname, MAXLEN,  "%s/%s_%03d.subcat", PARAMS.GROUP_DIR,PARAMS.GROUP_BASE,num);
  my_snprintf(subprop_fname,  MAXLEN,  "%s/%s_%03d.subprop", PARAMS.GROUP_DIR,PARAMS.GROUP_BASE,num);
#endif
  
  
  fcat = my_fopen(catalogue_fname, "r");
  my_fread(&Ngroups, sizeof(int64), one, fcat);

  if (Ngroups > 0) 
	{

		fpos  = my_fopen(particles_fname,"r");
		ftype = my_fopen(parttypes_fname, "r");
		fid   = my_fopen(partids_fname, "r");

#ifndef FOF_ONLY

		fprop    = my_fopen(fofprop_fname, "r");
		fsubcat  = my_fopen(subhalos_fname, "r");
		fsubprop = my_fopen(subprop_fname, "r");
#endif

		GroupLen    = my_malloc(sizeof(*GroupLen),Ngroups);
		GroupOffset = my_malloc(sizeof(*GroupOffset),Ngroups);

		my_fread(GroupLen   ,sizeof(int64),(size_t) Ngroups,fcat);
		my_fread(GroupOffset,sizeof(int64),(size_t) Ngroups,fcat);
	  
#ifndef FOF_ONLY
	  
		GroupSubs   = my_malloc(sizeof(*GroupSubs),Ngroups);//not recasting to (int64 *) -> should be automatic
		GroupMtot   = my_malloc(sizeof(*GroupMtot),Ngroups);
		GroupMgas   = my_malloc(sizeof(*GroupMgas),Ngroups);

		for(i=0;i<Ngroups;i++)
			my_fread(&GroupSubs[i],sizeof(int64),one,fcat);

		/* Read in from fof_prop file*/

		my_fread(&Ngroups,sizeof(int64),one,fprop);
		for(i=0;i<3;i++)
		{
			GroupCM[i]  = my_malloc(sizeof(*GroupCM[i]),Ngroups);
#ifdef GET_GROUPVEL
			GroupCMV[i] = my_malloc(sizeof(*GroupCMV[i]),Ngroups);
#endif

		}

		for(i=0;i<Ngroups;i++)
			for(j=0;j<3;j++)
				my_fread(&GroupCM[j][i],sizeof(float),one,fprop); /* Warning: This was written directly as a block of 3. I am reading it back in as 3 sets of 1. */

#ifdef GET_GROUPVEL
		for(i=0;i<Ngroups;i++)
			for(j=0;j<3;j++)
				my_fread(&GroupCMV[j][i],sizeof(float),one,fprop); /* Warning: This was written directly as a block of 3. I am reading it back in as 3 sets of 1. */
#endif


		for(i=0;i<Ngroups;i++)
			my_fread(&GroupMtot[i],sizeof(float),one,fprop);

		for(i=0;i<Ngroups;i++)
			my_fread(&GroupMgas[i],sizeof(float),one,fprop);

		fclose(fprop);

		/* Read in from subhalo catalogue*/
	  
		my_fread(&Nsub,sizeof(int64),one,fsubcat);
		SubLen    = my_malloc(sizeof(*SubLen),Nsub);
		SubOffset = my_malloc(sizeof(*SubOffset),Nsub);
		SubParent = my_malloc(sizeof(*SubParent),Nsub);

		for(i=0;i<Nsub;i++)
			my_fread(&SubLen[i],sizeof(int64),one,fsubcat);

		for(i=0;i<Nsub;i++)
			my_fread(&SubOffset[i],sizeof(int64),one,fsubcat);

		for(i=0;i<Nsub;i++)
			my_fread(&SubParent[i],sizeof(int64),one,fsubcat);

		fclose(fsubcat);


		/* read in from the subprop file*/
		my_fread(&Nsub,sizeof(int64),one,fsubprop);
		SubMtot = my_malloc(sizeof(*SubMtot),Nsub);
		SubMgas = my_malloc(sizeof(*SubMgas),Nsub);
		for(i=0;i<3;i++)
		{
			SubCM[i]  = my_malloc(sizeof(*SubCM[i]),Nsub);
#ifdef GET_GROUPVEL
			SubCMV[i] = my_malloc(sizeof(*SubCMV[i]),Nsub);
#endif
		}


		for(i=0;i<Nsub;i++)
			for(j=0;j<3;j++)
				my_fread(&SubCM[j][i],sizeof(float),one,fsubprop);

#ifdef GET_GROUPVEL
		for(i=0;i<Nsub;i++)
			for(j=0;j<3;j++)
				my_fread(&SubCMV[j][i],sizeof(float),one,fsubprop);
#endif


		for(i=0;i<Nsub;i++)
			my_fread(&SubMtot[i],sizeof(float),one,fsubprop);
	  
		for(i=0;i<Nsub;i++)
			my_fread(&SubMgas[i],sizeof(float),one,fsubprop);
	  
		fclose(fsubprop);

#endif	  

		/* Read in the positions  */
		my_fread(&NumPart,sizeof(int64),one,fpos);
		/* 	  fprintf(stderr,"There are " STR_FMT " particles in the groups file `%s'\n",NumPart,particles_fname); */
		my_fread(&NumPart,sizeof(int64),one,ftype);
		my_fread(&NumPart,sizeof(int64),one,fid);

	  
		xpos = my_malloc(sizeof(*xpos),NumPart);
		ypos = my_malloc(sizeof(*ypos),NumPart);
		zpos = my_malloc(sizeof(*zpos),NumPart);
		type = my_malloc(sizeof(*type),NumPart);
		id   = my_malloc(sizeof(*id),NumPart);

		for(i=0;i<NumPart;i++)
		{
			/* read positions */
			my_fread(&tmp_pos[0],sizeof(float),(size_t) 3,fpos);

			xpos[i] = tmp_pos[0];
			ypos[i] = tmp_pos[1];
			zpos[i] = tmp_pos[2];

			/* read types */
			my_fread(&type[i],sizeof(char),one,ftype);

			/* read ids */
			my_fread(&id[i],sizeof(id64),one,fid);
		}

		fclose(fcat);
		fclose(fpos);
		fclose(ftype);
		fclose(fid);
	  

		/* So everything has been read and all the files have been closed
			 Now, lets insert the subhalos into the proper places in the structure 
			 If we are reading FOF only halos then pretty much nothing needs to be done.
	  
			 Keep in mind that the subfind does return the FOF halos, all the particles 
			 present in the subhalos will add up to be <= the original FOF halo, i.e., the
			 position array will probably contain elements that do not belong to any subhalo.
	  
		*/


#ifndef FOF_ONLY

		/* output came from subfind */

		/* Decide which are the parent halos.
			 Notice that the loop below goes to Ngroups and NOT Nsub.
		*/


		tempindex = my_malloc(sizeof(*tempindex),Ngroups);
	  
		for(i=0;i<Nsub;i++)
		{
			group[i].isFof=0;
			group[i].Nsub = 0;
		}
	  
		sumgroups=0;
		for(i=0;i<Ngroups;i++)
		{
			tempindex[i] = sumgroups;
			sumgroups += GroupSubs[i];
		  
		}
	  
		for(i=0;i<Ngroups;i++) {
			if(tempindex[i] >= Nsub ) {
				fprintf(stderr,"WARNING: Possible case of FOF group falling below resolution limit \n");
				fprintf(stderr,"WARNING: i = %"STR_FMT" tempindex = %"STR_FMT" is greater than Nsub = %"STR_FMT" with Ngroups = %"STR_FMT"  \n",i,tempindex[i],Nsub,Ngroups);
			} else {
				group[tempindex[i]].isFof = 1;
				group[tempindex[i]].Nsub = GroupSubs[i];
			}
		}
		free(tempindex);

		FOF_Parent=0;
		for(i=0;i<Nsub;i++) {
			/* only the fof's have a definite
				 parentlevel at this point. all
				 others (subhalos) need to be assigned
				 one at a later stage.
			*/
	
			group[i].ParentLevel = -1;
			if(group[i].isFof==1) {
				group[i].ParentLevel = 1;
				FOF_Parent = i;
			}
	
			group[i].N = SubLen[i];
      group[i].nodeloc = i;
			group[i].snapshot = num;
			group[i].redshift = REDSHIFT[num];
			
			group[i].Mtot = SubMtot[i];
			group[i].Mgas = SubMgas[i];
			group[i].xcen  = xpos[SubOffset[i]];
			group[i].ycen  = ypos[SubOffset[i]];
			group[i].zcen  = zpos[SubOffset[i]];
	
#ifdef GET_GROUPVEL
			group[i].vxcen = SubCMV[0][i];
			group[i].vycen = SubCMV[1][i];
			group[i].vzcen = SubCMV[2][i];
#endif
			group[i].groupnum = i;
			group[i].N_per_wedge = 0;
			/* initialise the parent finding variables*/
			group[i].ParentId =-1;
			group[i].NParents = 0;
			group[i].Switched = 0;
			group[i].ParentSnapshot = -1;
			group[i].Ncommon = 0;
			group[i].Rank = 0.0;
			group[i].NpartinParent = 0;
			group[i].snapshot = (short) num;
			group[i].FOFHalo = FOF_Parent;
			group[i].ContainerIndex = FOF_Parent;
	
	
			group[i].x    =  my_malloc(sizeof(*group[i].x),SubLen[i]);
			group[i].y    =  my_malloc(sizeof(*group[i].y),SubLen[i]);
			group[i].z    =  my_malloc(sizeof(*group[i].z),SubLen[i]);
			group[i].id   =  my_malloc(sizeof(*group[i].id),SubLen[i]);
			group[i].type =  my_malloc(sizeof(*group[i].type),SubLen[i]);
	
			group[i].parentgroupforparticle    = my_malloc(sizeof(*group[i].parentgroupforparticle),SubLen[i]);
			group[i].parentsnapshotforparticle = my_malloc(sizeof(*group[i].parentsnapshotforparticle),SubLen[i]);
	
			rmax = 0.0;
	
			for(j=0;j<SubLen[i];j++) {
				group[i].x[j] = xpos[SubOffset[i]+j];
				group[i].y[j] = ypos[SubOffset[i]+j];
				group[i].z[j] = zpos[SubOffset[i]+j];
				group[i].id[j] = id[SubOffset[i]+j];
				group[i].type[j] = (int) type[SubOffset[i]+j];
	  
				/* 			  fprintf(stderr,"group[%d].id[%d]=%ld\n",i,j,group[i].id[j]); */
	  
				temp_rmax = sqrt(SQR_PERIODIC(group[i].x[j],group[i].xcen) + SQR_PERIODIC(group[i].y[j],group[i].ycen) + SQR_PERIODIC(group[i].z[j],group[i].zcen) );
				if( temp_rmax > rmax)
					rmax = temp_rmax;
	  
				/* initialise the parent finding variables*/
				group[i].parentgroupforparticle[j]    = -1;
				group[i].parentsnapshotforparticle[j] = -1;
	  
			}
			group[i].Rmax = rmax;
		}
      
#else
		/* output came from fof groups */

		for(i=0;i<Ngroups;i++) {
			group[i].N = GroupLen[i];
			group[i].NpartinParent = 0;
			group[i].nodeloc = i;
			group[i].snapshot = num;
			group[i].redshift = REDSHIFT[num];
			
			group[i].x    = my_malloc(sizeof(*group[i].x),GroupLen[i]);
			group[i].y    = my_malloc(sizeof(*group[i].y),GroupLen[i]);
			group[i].z    = my_malloc(sizeof(*group[i].z),GroupLen[i]);
			group[i].id   = my_malloc(sizeof(*group[i].id),GroupLen[i]);
			group[i].type = my_malloc(sizeof(*group[i].type),GroupLen[i]);
	
			group[i].parentgroupforparticle    = my_malloc(sizeof(*group[i].parentgroupforparticle),GroupLen[i]);
			group[i].parentsnapshotforparticle = my_malloc(sizeof(*group[i].parentsnapshotforparticle),GroupLen[i]);
	
			/* Fake value since all the FOF are considered parents. Almost */ 
			group[i].groupnum = i;
			group[i].NParents = 0;
			group[i].Switched = 0;
			group[i].isFof = 1;
			group[i].Nsub = 0;
			group[i].FOFHalo = i;
			group[i].ParentLevel = 1;
			group[i].ContainerIndex = i;
	
			group[i].N_per_wedge = 0;
			/* initialise the parent finding variables*/
			group[i].ParentId =-1;
			group[i].ParentSnapshot = -1;
			group[i].Ncommon = 0;
			group[i].Rank = 0.0;
			group[i].snapshot = (short) num;
	
			/*WARNING: Place-holder for future computation of actual cm in FOF_ONLY case*/
			group[i].xcen  = xpos[GroupOffset[i]];
			group[i].ycen  = ypos[GroupOffset[i]];
			group[i].zcen  = zpos[GroupOffset[i]];
	
#ifdef GET_GROUPVEL
			group[i].vxcen = GroupCMV[0][i];
			group[i].vycen = GroupCMV[1][i];
			group[i].vzcen = GroupCMV[2][i];
#endif 
			mtot=0.0;
			mgas=0.0;
	
			/* This rmax should not be taken too seriously for FOF groups. The center needs to be re-calculated */
			rmax = 0.0;
	
			for(j=0;j<GroupLen[i];j++) {
				group[i].x[j]    = xpos[GroupOffset[i]+j]; 
				group[i].y[j]    = ypos[GroupOffset[i]+j]; 
				group[i].z[j]    = zpos[GroupOffset[i]+j]; 
				group[i].id[j]   = id[GroupOffset[i]+j]; 
				group[i].type[j] = (int) type[GroupOffset[i]+j]; 
	  
				/* 			  mtot+=PARAMS.MASSARR[group[i].type[j]]; */
				/* 			  if (group[i].type[j] == 0) */
				/* 				mgas+=PARAMS.MASSARR[0]; */
	  
				temp_rmax = sqrt(SQR_PERIODIC(group[i].x[j],group[i].xcen) + SQR_PERIODIC(group[i].y[j],group[i].ycen) + SQR_PERIODIC(group[i].z[j],group[i].zcen) );
				if( temp_rmax > rmax)
					rmax = temp_rmax;
	  
				group[i].parentgroupforparticle[j]    = -1;
				group[i].parentsnapshotforparticle[j] = -1;
	  
			}
			/* 		  group[i].Mtot = mtot; */
			/* 		  group[i].Mgas = mgas; */
			group[i].Mtot = GroupMtot[i];
			group[i].Mgas = GroupMgas[i];
			group[i].Rmax = rmax;
		}
#endif	  


		/*free up memory*/
		free(xpos);
		free(ypos);
		free(zpos);
		free(type);
		free(id);
		free(GroupLen);
		free(GroupOffset);
	  
#ifdef GET_GROUPVEL
		for(i=0;i<3;i++) {
			free(GroupCM[i]);
			free(GroupCMV[i]);
		}
#endif
	  
#ifndef FOF_ONLY
			  
		free(SubMtot);
		free(SubMgas);
		free(SubLen);
		free(SubOffset);
		free(SubParent);
		free(GroupSubs);
		free(GroupMtot);
		free(GroupMgas);
	  
#ifdef GET_GROUPVEL
		for(i=0;i<3;i++) {
			free(SubCM[i]);
			free(SubCMV[i]);
		}
#endif


#endif
	  
	} else {
		fclose(fcat);
		fprintf(stderr,"No groups found in file..freeing the group pointer \n");
		my_free((void **) &(group));
	}
}



#endif



#undef    FOFMINLEN
#undef    SUBHALOMINLEN


