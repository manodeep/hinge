#include "defs.h"
#include "read_param.h"
#include "utils.h"
#include "missinghalos.h"
#include "loadgroups.h"
#include "sglib.h"

int64 binary_search(id64 key,id64 *arr,int64 size);
void print_header_ghost_halos(FILE *fp);
int comparator_for_bsearch_id64(const void *key, const void *member);

int comparator_for_bsearch_id64(const void *key, const void *member)
{
  const id64 *keyptr = key;
  const id64 *memptr = member;
  const id64 keyval = *keyptr;
  const id64 memval = *memptr;

  return (keyval < memval) ? -1:( (keyval == memval) ?  0:1 );
}


int64 binary_search(id64 key,id64 *arr,int64 size)
{
  int64 left, right, middle;
  left = 0;
  right = size-1;
  
  while(left <= right)
	{
	  middle = (int64) (left + right)/2;
	  
	  if(key == arr[middle])
		return(middle);
	  
	  if(key > arr[middle])
		left = middle+1;
	  else
		right = middle-1;
		
	}

  return(-1);
}


void print_header_ghost_halos(FILE *fp)
{

  fprintf(fp,"####################################################################################################################################\n");
  fprintf(fp,"#     snapshot         haloid          xcen            ycen              zcen              vxcen          vycen             vzcen   \n");
  fprintf(fp,"#       i                 l              f               f                f                  f              f                 f     \n");
  fprintf(fp,"####################################################################################################################################\n");
}



void output_missing_halo_centres(struct node_data *tree[],int64 *Ngroups)
{
  struct node_data *BaseNode=NULL,*thisnode=NULL,*tmpnode=NULL;
  struct group_data *group0=NULL;

  char fname[MAXLEN];
  char snapshot_name[MAXLEN];

  id64 *TrackIds=NULL;
  int64 *TrackNums=NULL; //Just a fake groupnumber
  int64 *TrackHaloIds=NULL;
  int64 *Track_NumPart_in_Group=NULL;
  short *EarliestSnapshot_to_Track=NULL;
  float *SubCM[3],*SubCMV[3];
  float *zero_part_pos[3];

  int64 *Nfound_in_Group=NULL;


  short EarliestSnapshot = PARAMS.MAX_SNAPSHOT_NUM;
  int64 Nids=0;
  int64 Nmissing=0;
  int64 npart=0;
  int64 offset=0;
  struct io_header header1;
  float tmp_pos[3],tmp_vel[3];
  
  id64 tmp_id;

  int64 location,nmisses;
  int64_t TOTNPART=NUMPART;
  int64_t check_totnpart=0;
  FILE *fdpos=NULL,*fdvel=NULL,*fdids=NULL;
  int32_t dummy,ifile,nfiles;
  short prevsnapshot;
  int64 index,totnleft;

  FILE *fp=NULL;
  int64 AllNmissing[PARAMS.MAX_SNAPSHOT_NUM+1];
  for(short i=PARAMS.MIN_SNAPSHOT_NUM;i<=PARAMS.MAX_SNAPSHOT_NUM;i++) {
		AllNmissing[i] = 0;
	}

  my_snprintf(fname,MAXLEN,"rm -f %s/ghost_halos_???.txt",PARAMS.OUTPUT_DIR);
  fprintf(stderr,"executing system command `%s' ",fname);
  system(fname);//remove all the pre-existing ghost_halos_???.txt files.
  fprintf(stderr,"...done\n");

  my_snprintf(snapshot_name,MAXLEN,"%s/%s_%03d",PARAMS.SNAPSHOT_DIR,PARAMS.SNAPSHOT_BASE,PARAMS.MAX_SNAPSHOT_NUM);
  header1 = get_gadget_header(snapshot_name);
  nfiles = header1.num_files;

  check_totnpart = 0;
  for(int i = 0; i < 6; i++)
    {
      check_totnpart += header1.npartTotal[i];
      check_totnpart += (((int64_t) header1.npartTotalHighWord[i]) << 32);
    }
  
  if(check_totnpart != TOTNPART)
    {
      fprintf(stderr,"ERROR: COSMO struct is not initialized correctly\n");
      fprintf(stderr,"Got TOTNPART  = %"PRId64" from COSMO struct but snapshot header shows %"PRId64" particles\n",TOTNPART,check_totnpart);
      fprintf(stderr,"exiting...\n");
      exit(EXIT_FAILURE);
    }


  for(short isnapshot=PARAMS.MAX_SNAPSHOT_NUM;isnapshot>PARAMS.MIN_SNAPSHOT_NUM;isnapshot--)
	{
	  prevsnapshot = isnapshot-1;
	  if(prevsnapshot == PARAMS.MIN_SNAPSHOT_NUM) break;

	  if(Ngroups[isnapshot] > 0)
		{
		  Nids=0;
		  offset = Nids;
		  Nmissing=0;
		  EarliestSnapshot = PARAMS.MAX_SNAPSHOT_NUM;
		  BaseNode = tree[isnapshot];
		  group0 = allocate_group(Ngroups[isnapshot]);
		  loadgroups(isnapshot,group0);
		  
		  for(int64 igroup=0;igroup<Ngroups[isnapshot];igroup++)
			{
			  thisnode=&BaseNode[igroup];
			  tmpnode = thisnode->BigChild;
			  
			  if(tmpnode != NULL && (thisnode->snapshot-tmpnode->snapshot) > 1)
				{
				  Nmissing++;
				  
/* 				  Nids +=  (group0[igroup].N > MAXPARTLOC) ? MAXPARTLOC:group0[igroup].N; */

				  if(PARAMS.MAX_RANK_LOC <= 0)
					{
					  Nids += group0[igroup].N;
					}
				  else
					{
					  Nids += (PARAMS.MAX_RANK_LOC < group0[igroup].N) ? PARAMS.MAX_RANK_LOC:group0[igroup].N; 
					}

				  //tmpnode exists => I have to track the group to tmpnode->snapshot + 1
				  if( (tmpnode->snapshot+1) < EarliestSnapshot) 
					EarliestSnapshot = tmpnode->snapshot+1;

				  TrackIds=my_realloc(TrackIds,sizeof(id64),Nids,"TrackIds");
				  TrackNums=my_realloc(TrackNums,sizeof(int64),Nids,"TrackNums");
				  Track_NumPart_in_Group = my_realloc(Track_NumPart_in_Group,sizeof(int64),Nids,"Track_NumPart_in_Group");
				  EarliestSnapshot_to_Track=my_realloc(EarliestSnapshot_to_Track,sizeof(short),Nids,"EarliestSnapshot_to_Track");

				  for(short missing_snapshot=isnapshot-1;missing_snapshot>=tmpnode->snapshot+1;missing_snapshot--) {
						AllNmissing[missing_snapshot]++;
					}
				  
				  TrackHaloIds    = my_realloc(TrackHaloIds,sizeof(int64),Nmissing,"TrackHaloIds");
				  TrackHaloIds[Nmissing-1] = thisnode->haloid;

				  for(int64 k=offset;k<Nids;k++)
					{
					  TrackIds[k] = group0[igroup].id[k-offset];
					  TrackNums[k] = Nmissing-1;
					  Track_NumPart_in_Group[k] = k-offset;
					  EarliestSnapshot_to_Track[k] = tmpnode->snapshot+1;
					}


				  offset = Nids;
				}//tmpnode != NULL
			}//igroup loop
		  
		  free_group(group0,Ngroups[isnapshot]);

		  
		  //new scope
		  {
			int64 bad=0;
			for(int64 k=0;k<Nids;k++)
			  {
				if(TrackIds[k] <= 0 || TrackIds[k] > TOTNPART)
				  bad++;
				
/* 				fprintf(stderr,"TrackIds[%"STR_FMT"] = %"STR_ID_FMT" TrackNums[%"STR_FMT"] = %"STR_FMT"\n", */
/* 						k,TrackIds[k],k,TrackNums[k]); */
			  }
			if(bad > 0)
			  {
				fprintf(stderr,"ERROR: number of particleids less than 0 or greater than TOTNPART (%"PRId64") = %"STR_FMT" ..exiting\n",TOTNPART,bad);
				exit(EXIT_FAILURE);
			  }
		  }
		  
		  //Now all the particles that need to be tracked have been fit into TrackIds and TrackNums
		  //Use the previous snapshot 
		  for(int i=0;i<3;i++)
			{
			  SubCM[i]  = my_calloc(sizeof(*SubCM[i]),Nmissing);
			  SubCMV[i] = my_calloc(sizeof(*SubCMV[i]),Nmissing);
			  zero_part_pos[i] = my_calloc(sizeof(*zero_part_pos[i]),Nmissing);
			}
		  Nfound_in_Group = my_calloc(sizeof(*Nfound_in_Group),Nmissing);


		  //sort the missing particles such that particles that need to be tracked the longest appear at the front

#define MULTIPLE_ARRAY_EXCHANGER(type,a,i,j) { SGLIB_ARRAY_ELEMENTS_EXCHANGER(short,EarliestSnapshot_to_Track,i,j); SGLIB_ARRAY_ELEMENTS_EXCHANGER(int64,TrackNums,i,j); \
			SGLIB_ARRAY_ELEMENTS_EXCHANGER(id64,TrackIds,i,j);SGLIB_ARRAY_ELEMENTS_EXCHANGER(int64,Track_NumPart_in_Group,i,j); }

		  
		  //There probably are particles that need to be tracked more than one snapshot. Track them through
		  //successive snapshots until they are done. 
		  short snapshot=prevsnapshot;
		  while(snapshot >= EarliestSnapshot && snapshot > PARAMS.MIN_SNAPSHOT_NUM && Nids > 0)
			{
			  //sort the array by earliest snapshot
			  SGLIB_ARRAY_QUICK_SORT(short, EarliestSnapshot_to_Track, Nids, SGLIB_NUMERIC_COMPARATOR, MULTIPLE_ARRAY_EXCHANGER);

			  npart=0;
			  //count how many particles we need to track at this `snapshot' 
			  for(int64 k=0;k<Nids;k++)
				if(EarliestSnapshot_to_Track[k] <= snapshot)
				  npart++;
				else
				  break;
			  
			  //new scope -- make sure that the sorting in earliest snapshot works ..(as I think it should -- earliest snapshots in the beg. and later ones
			  // at the end. so the later snapshots/particles willl get done earlier and then I can merrily realloc and discard those last parts
			  {
				int64 bad=0;
				for(int64 k=npart;k<Nids;k++)
				  if(EarliestSnapshot_to_Track[k] <= snapshot)
					{
					  fprintf(stderr,"This does not work. I have earliest snapshot to track[%"STR_FMT"] = %d\n",k,EarliestSnapshot_to_Track[k]);
					  bad++;
					}
				
				if(bad > 0) 
				  {
					fprintf(stderr,"nbad = %"STR_FMT" exiting\n",bad);
					exit(EXIT_FAILURE);
				  }
			  }//done with the limited scope
			  
			  //So the previous check worked - reset Nids to contain the actual number of particles that need to be
			  //checked at this snapshot (including earlier snapshots - those will come in the next iteration of the snapshot while loop)
			  Nids = npart;
			  TrackIds=my_realloc(TrackIds,sizeof(id64),Nids,"TrackIds");
			  TrackNums=my_realloc(TrackNums,sizeof(int64),Nids,"TrackNums");
			  Track_NumPart_in_Group = my_realloc(Track_NumPart_in_Group,sizeof(int64),Nids,"Track_NumPart_in_Group");
			  EarliestSnapshot_to_Track=my_realloc(EarliestSnapshot_to_Track,sizeof(short),Nids,"EarliestSnapshot_to_Track");

			  //Now sort all the arrays back by particle id (TrackIds)
			  SGLIB_ARRAY_QUICK_SORT(id64, TrackIds, Nids, SGLIB_NUMERIC_COMPARATOR, MULTIPLE_ARRAY_EXCHANGER);


/* 			  fprintf(stderr,"Working on snapshot = %4d Number of missing particles = %12"STR_FMT" earliest snapshot to track to = %4d\n", */
/* 					  snapshot,Nids,EarliestSnapshot); */
			  
			  //Now lets initialise all the variables
			  for(int64 j=0;j<Nmissing;j++)
				{
				  Nfound_in_Group[j]  = 0;
				  for(int i=0;i<3;i++)
					{
					  SubCM[i][j]=0.0;
					  SubCMV[i][j]=0.0;
					  zero_part_pos[i][j] = 0.0;
					}
				  
				}
			  
			  ifile=0;
			  totnleft = Nids;
			  while(ifile < nfiles && totnleft > 0)
				{
				  if (nfiles == 1) {
						my_snprintf(snapshot_name,MAXLEN,"%s/%s_%03d",PARAMS.SNAPSHOT_DIR,PARAMS.SNAPSHOT_BASE,snapshot);
					} else {
						my_snprintf(snapshot_name,MAXLEN,"%s/%s_%03d.%d",PARAMS.SNAPSHOT_DIR,PARAMS.SNAPSHOT_BASE,snapshot,ifile);
					}

				  fprintf(stderr,"snapshotfile = `%s' Nmissing = %"STR_FMT" Nids = %"STR_FMT"  totnleft = %"STR_FMT" Earliestsnap = %hd \n",
						  snapshot_name,Nmissing,Nids,totnleft,EarliestSnapshot);


				  fdpos=my_fopen(snapshot_name,"r");
				  my_fread(&dummy, sizeof(dummy), 1, fdpos);
				  my_fread(&header1, sizeof(header1), 1, fdpos);
				  my_fread(&dummy, sizeof(dummy), 1, fdpos);
				  my_fread(&dummy, sizeof(dummy), 1, fdpos);
				  
				  fdvel=my_fopen(snapshot_name,"r");
				  my_fread(&dummy, sizeof(dummy), 1, fdvel);
				  my_fread(&header1, sizeof(header1), 1, fdvel);
				  my_fread(&dummy, sizeof(dummy), 1, fdvel);
				  my_fread(&dummy, sizeof(dummy), 1, fdvel);
				  for(int k=0; k<6; k++) /* skip pos data */
					fseek(fdvel, header1.npart[k]*sizeof(float)*3, SEEK_CUR);
				  
				  my_fread(&dummy, sizeof(dummy), 1, fdvel);
				  my_fread(&dummy, sizeof(dummy), 1, fdvel);
				  
				  fdids=my_fopen(snapshot_name,"r");
				  my_fread(&dummy, sizeof(dummy), 1, fdids);
				  my_fread(&header1, sizeof(header1), 1, fdids);
				  my_fread(&dummy, sizeof(dummy), 1, fdids);
				  my_fread(&dummy, sizeof(dummy), 1, fdids);
				  for(int k=0; k<6; k++) /* skip pos data */
					fseek(fdids, header1.npart[k]*sizeof(float)*3, SEEK_CUR);
				  my_fread(&dummy, sizeof(dummy), 1, fdids);
				  my_fread(&dummy, sizeof(dummy), 1, fdids);
				  for(int k=0; k<6; k++) /* skip vel data */
					fseek(fdids, header1.npart[k]*sizeof(float)*3, SEEK_CUR);
				  my_fread(&dummy, sizeof(dummy), 1, fdids);
				  my_fread(&dummy, sizeof(dummy), 1, fdids);

				  nmisses=0;
				  index=0;
				  for(int k=0; k<6; k++)
					{ 
					  for(int n=0; n<header1.npart[k]; n++)
						{ 
						  my_fread(&tmp_pos[0],sizeof(float),3,fdpos);
						  my_fread(&tmp_vel[0],sizeof(float),3,fdvel);
						  my_fread(&tmp_id,sizeof(id64),1,fdids);
						  if(tmp_id <= 0 || tmp_id > TOTNPART)
							{
							  fprintf(stderr,"ERROR: I have id = %"STR_ID_FMT" but it should be between 0 and %"PRId64" \n",tmp_id,TOTNPART);
							  exit(EXIT_FAILURE);
							}
						  
						  index++;
						  location = binary_search(tmp_id,TrackIds,Nids);

						  if(location != -1)
							{
							  if(location >= Nids || location < 0)
								{
								  fprintf(stderr,"ERROR: location is beyond array range location = %"STR_FMT" tmp_id = %"STR_ID_FMT" index = %"STR_FMT"\n",location,tmp_id,index);
								  exit(EXIT_FAILURE);
								}
							  if(TrackIds[location] != tmp_id)
								{
								  fprintf(stderr,"ERROR: Did not get the correct location.. TrackIds[%"STR_FMT"] = %"STR_ID_FMT"  tmp_id = %"STR_ID_FMT"\n",
										  location,TrackIds[location],tmp_id);
								  exit(EXIT_FAILURE);
								}
							  else
								{
								  totnleft--;
								  Nfound_in_Group[TrackNums[location]]++;
/* 								  fprintf(stderr,"totnleft = %"STR_FMT" TrackNums[%"STR_FMT"] = %"STR_FMT" Nfound = %"STR_FMT"\n", */
/* 										  totnleft,location,TrackNums[location],Nfound_in_Group[TrackNums[location]]); */
								  if(Nfound_in_Group[TrackNums[location]] == 1)
									{
									  for(int ii=0;ii<3;ii++)
										zero_part_pos[ii][TrackNums[location]] = tmp_pos[ii];
										  
									}
								
								  for(int ii=0;ii<3;ii++)
									{
									  SubCM[ii][TrackNums[location]] += periodic(tmp_pos[ii] - zero_part_pos[ii][TrackNums[location]]);
									  SubCMV[ii][TrackNums[location]] += tmp_vel[ii];
									}
								
								  
								}	//so I really found that particle TrackIds[location] == tmp_id					  
							}// location != -1
						  else
							nmisses++;

						}//n=0, n<header1.npart[k]
					}//k=0,k<6


				  fclose(fdpos);
				  fclose(fdvel);
				  fclose(fdids);
				  ifile++;

				  fprintf(stderr,"after reading the snapshot: totnleft = %"STR_FMT"  index = %"STR_FMT" nmisses=%"STR_FMT" index-nmisses=%"STR_FMT"\n",totnleft,index,nmisses,index-nmisses);
				}//while ifile < nfiles
			  
			  for(int j=0;j<3;j++)
				{
				  for(int64 i=0;i<Nmissing;i++)
					{
					  if(Nfound_in_Group[i] > 0)
						{
						  SubCM[j][i] /= (double) Nfound_in_Group[i];
						  SubCM[j][i] += zero_part_pos[j][i];
						  SubCM[j][i] = periodic_wrap(SubCM[j][i]);
						  
						  SubCMV[j][i] /= (double)Nfound_in_Group[i];
						}
					}
				}


			  //output the data
			  my_snprintf(fname,MAXLEN,"%s/ghost_halos_%03d.txt",PARAMS.OUTPUT_DIR,snapshot);
			  fp = my_fopen_carefully(fname,&print_header_ghost_halos);

			  for(int64 i=0;i<Nmissing;i++)
				if(Nfound_in_Group[i] > 0)
				  fprintf(fp,"%10hd     %14"STR_FMT"   %14.4f   %14.4f   %14.4f   %14.4f   %14.4f   %14.4f\n",snapshot,TrackHaloIds[i],SubCM[0][i],SubCM[1][i],SubCM[2][i],
						  SubCMV[0][i],SubCMV[1][i],SubCMV[2][i]);

/* 				else */
/* 				  fprintf(stderr,"Nfound_in_group[%"STR_FMT"] = %"STR_FMT"\n",i,Nfound_in_Group[i]); */

			  
			  fclose(fp);
			  snapshot--;
			}//snapshot > EarliestSnapshot
		  

		  my_free((void **) &TrackIds);
		  my_free((void **) &TrackNums);
		  my_free((void **) &Track_NumPart_in_Group);
		  my_free((void **) &EarliestSnapshot_to_Track);

		  for(int i=0;i<3;i++)
			{
			  my_free((void **) &(SubCM[i]));
			  my_free((void **) &(SubCMV[i]));
			  my_free((void **) &(zero_part_pos[i]));
			}

		  my_free((void **) &Nfound_in_Group);
		  my_free((void **) &TrackHaloIds);

		}//Ngroups > 0
	}//isnapshot loop
  
  for(short isnapshot=PARAMS.MIN_SNAPSHOT_NUM;isnapshot<=PARAMS.MAX_SNAPSHOT_NUM;isnapshot++)
	fprintf(stderr,"AllNmissing[%d]  = %"STR_FMT"\n",isnapshot,AllNmissing[isnapshot]);

  
}
