#include "../header/hpfem.h"
#include "../header/exvar.h"
#include "../header/repartition_BSFC.h"
//#define DEBUG

/* decent values for this partitioning scheme, the user can
   set them in the as parameters to tune for better performance */
#define BINS_PER_PROC 50 /* minimum amount of coarse bins on each processor */
#define SUBBINS_PER_BIN 50 /* amount of subbins a bin is divided into */
#define MAX_REFINEMENT_LEVEL 10 /* amount of refinement of the bins */

// stuff for the load-balancing weights used
#define NON_EMPTY_CELL (double) 1.4
#define EMPTY_CELL (double) 1


void BSFC_create_refinement_info(int* number_of_cuts, 
				 float* global_actual_work_allocated,
				 float total_weight, 
				 float* work_percent_array,
				 unstructured_communication verts_in_cuts_info,
				 float**, int myid, int numprocs);

void BSFC_create_bins(int num_local_objects,
		      BSFC_VERTEX_PTR sfc_vert_ptr, 
		      int* amount_of_bits_used, int size_of_unsigned,
		      float* global_actual_work_allocated, 
		      float *work_percent_array, float* total_weight_ptr,
		      int* balanced_flag, unstructured_communication* verts_in_cuts_info,
		      int* number_of_cuts,  
		      int bins_per_proc,
		      int myid, int numprocs);

void BSFC_update_element_proc(int myid, int numprocs, 
			      HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
			      BSFC_VERTEX_PTR sfc_vert_ptr);

/* Space filling curve (BSFC) partioning routine */

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void repartition(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, int time_step)
{
  int ierr, i, j, k;                     /* local variables */
  int num_local_objects;              /* the number of objects this processor owns */
  BSFC_VERTEX_PTR sfc_vert_ptr;        /* array that stores the sfc objects */
  float* objs_wgt = NULL;             /* array of objects weights */
  float *global_actual_work_allocated = NULL; /* cumulative actual work allocated */
  float total_weight = 0;   /* sum of the work i.e. the amount of work for all objects */
  float *work_percent_array = NULL;   /* the cumulative percent of work each
					 processor should ideally get */
  int balanced_flag;                  /* flag to indicate if all partitions
					 are balanced */
  float* wgts_in_cut_ptr = NULL;      /* array of weights for sfc objects in 
					 a cut bin */
  int number_of_cuts = 0; /* maximum amount of cuts in a coarse bin on this processor */
  int amount_of_bits_used = 0;        /* amount of bits used in calculating the
					 bin an sfc object belongs to */
  int refinement_level_counter = 0;   /* counter to keep track of how many 
					 levels of bin refinement have been performed */
  int myid, numprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  int size_of_unsigned = sizeof(unsigned);
  //printf("proc %d has entered the repartitioning scheme...\n",myid);
  
  /* get application data (number of objects, ids, weights, and coords */
  int elm_counter = 0;  //used to keep count how many active elements are on this processor
  int no_of_buckets = HT_Elem_Ptr->get_no_of_buckets();
  HashEntryPtr entryp;
  Element* EmTemp;
  for(i=0;i<no_of_buckets;i++)
    {
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
	{
	  EmTemp = (Element*)(entryp->value);
#ifdef DEBUG	   
	    printf("lbweight %e",EmTemp->get_lb_weight());
	    assert(EmTemp->get_lb_weight()< 1000000000000000);
#endif	    
	  if(!EmTemp->get_refined_flag()) {
	    EmTemp->put_lb_weight((double) EmTemp->get_no_of_dof());

	    EmTemp->put_myprocess(-1);
	    
	    total_weight += EmTemp->get_lb_weight();
	    EmTemp->put_new_old(BSFC_NEW);
	    elm_counter++;
	  }
	  entryp = entryp->next;
	}
    }
  //printf("there are %d active elements on proc %d\n",elm_counter,myid);
  //elements that share constrained nodes are grouped into 1 objects (since they cannot be separated onto different processors)  
  num_local_objects = 0; 
  for(i=0;i<no_of_buckets;i++)
    {
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
	{
	  EmTemp = (Element*)(entryp->value);
#ifdef DEBUG	   
	    printf("lbweight %e",EmTemp->get_lb_weight());
	    assert(EmTemp->get_lb_weight()< 1000000000000000);
#endif

	  if(!EmTemp->get_refined_flag() && EmTemp->get_new_old() == BSFC_NEW) {
	    //check for constrained nodes on the vertex nodes
	    j = 4;
	    while(j<8)  {
	      Node* ndtemp = (Node*) HT_Node_Ptr->lookup((EmTemp->getNode()+j*KEYLENGTH));
	      if(ndtemp->getinfo() == S_S_CON) {
		BSFC_combine_elements(j-4, EmTemp, HT_Elem_Ptr, HT_Node_Ptr, -1);
		j = 10;  //exit out of the loop because we found a constrained node...
	      }
	      j++;
	    }
	    if(j==8)  // no constrained nodes...
	      EmTemp->copy_key_to_lb_key();
		
	    
	    num_local_objects++;
// 	    if(time_step == 300 && myid == 0) {
// 	      printf("%d is %u %u \n",num_local_objects, *(EmTemp->get_lb_key()), *(EmTemp->get_lb_key()+1));
// 	    }
	  }
	  entryp = entryp->next;
	}
    }

  
  sfc_vert_ptr = (BSFC_VERTEX_PTR) malloc(num_local_objects * sizeof(BSFC_VERTEX));
  // fill up the sfc_vert_ptr array which stores all the necessary info about the load-balancing objects
  j = 0; 
  for(i=0;i<no_of_buckets;i++)
    {
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
	{
	  EmTemp = (Element*)(entryp->value);
	  if(!EmTemp->get_refined_flag() && EmTemp->get_new_old() > 0) {
// 	    if(j == 1489)
// 	      j = j;
	    sfc_vert_ptr[j].lb_weight = EmTemp->get_lb_weight();
#ifdef DEBUG
	    assert(EmTemp->get_lb_weight()<1000000000000);
#endif
	    unsigned* elem_key = EmTemp->get_lb_key();
	    for(k=0;k<KEYLENGTH;k++)
	      sfc_vert_ptr[j].sfc_key[k] = elem_key[k];
	 
	    j++;
// 	    if(time_step == 300 && myid == 0) {
// 	      printf("%d is %u %u \n",j, *(EmTemp->get_lb_key()), *(EmTemp->get_lb_key()+1));
// 	    }
	  }
	  entryp = entryp->next;
	}
    }    
  assert(j == num_local_objects);
  global_actual_work_allocated=(float*) malloc(sizeof(float)*(numprocs+1));

  work_percent_array = (float*) malloc(sizeof(float) * numprocs);
  unstructured_communication verts_in_cut_info;
  verts_in_cut_info.used_flag = 0;
  /*create bins, fill global weight vector and perform initial partition of the bins*/
  BSFC_create_bins(num_local_objects, sfc_vert_ptr, 
		   &amount_of_bits_used, size_of_unsigned,
		   global_actual_work_allocated, work_percent_array, 
		   &total_weight, &balanced_flag, &verts_in_cut_info,
		   &number_of_cuts, BINS_PER_PROC, myid, numprocs);
 
  if(balanced_flag != BSFC_BALANCED) { 
    int* local_balanced_flag_array; /* used to indicate which partitions on this 
				       processor are already balanced - useful
				       for when more than one cut in a bin */
    int* ll_bins_head; /* used to indicate the beginning of the linklist */
    float* work_prev_allocated = NULL; /* stores the weights of all 
					  objects before a cut */
    
    if(verts_in_cut_info.recv_count == 0 || myid == 0) 
      balanced_flag = BSFC_BALANCED;
    BSFC_create_refinement_info(&number_of_cuts, 
				global_actual_work_allocated, 
				total_weight, work_percent_array,
				verts_in_cut_info, &work_prev_allocated,
				myid, numprocs);
    
    ll_bins_head = (int*) malloc(sizeof(int) * (1+number_of_cuts));
    
    if(number_of_cuts == 0)
      balanced_flag = BSFC_BALANCED;
    
    if(ll_bins_head != NULL)
      ll_bins_head[number_of_cuts] = 0;  /* the first link list starts off
					    at array location number_of_cuts-1 ! */
    for(i=0;i<number_of_cuts;i++)
      ll_bins_head[i] = -1;
    
    local_balanced_flag_array = (int*) malloc(sizeof(int) * (1+number_of_cuts));
    
    for(i=0;i<number_of_cuts;i++)
      local_balanced_flag_array[i] = BSFC_BALANCED;
    local_balanced_flag_array[number_of_cuts] = balanced_flag;
    
    /* refine bins until a satisfactory partition tolerance is obtained */
    while(balanced_flag != BSFC_BALANCED &&
	  refinement_level_counter < MAX_REFINEMENT_LEVEL) { 
      BSFC_refine_partition(&balanced_flag, 
			    &amount_of_bits_used, verts_in_cut_info.recv_count,
			    verts_in_cut_info.recv_sfc_vert,  
			    work_percent_array, total_weight,
			    global_actual_work_allocated, 
			    number_of_cuts,
			    ll_bins_head, work_prev_allocated, 
			    SUBBINS_PER_BIN, local_balanced_flag_array, myid, numprocs);
      refinement_level_counter++;
    }
    // printf("proc %d has refined %d levels===============\n",myid,refinement_level_counter);
    free(local_balanced_flag_array);
    free(work_prev_allocated);
    free(ll_bins_head);
  }
  
  /* if the load-balancing objects were moved to different processors,
     we need to move them back now */
 
  if(verts_in_cut_info.used_flag != 0) {  // move back and free up the space...
    int recv_count = 0;
    sfc_vertex* send_sfc_vert = new sfc_vertex[verts_in_cut_info.send_count];  //this array is actually for receiving the info...
    // fill up the send array...
    int* proc_counter = new int[numprocs];
    proc_counter[0] = 0;
    for(i=1;i<numprocs;i++)
      proc_counter[i] = proc_counter[i-1] + verts_in_cut_info.send_procs_ptr[i-1];
    //done filling up the send array
    int tag = 21504;
    
    recv_count = 0;
    MPI_Request* send_request = new MPI_Request[numprocs];
    MPI_Request* recv_request = new MPI_Request[numprocs];
    for(i=0;i<numprocs;i++) {
      if(i!= myid) {
	//  receive necessary info here...
	if(verts_in_cut_info.send_procs_ptr[i] != 0) {
	  j = MPI_Irecv((send_sfc_vert+proc_counter[i]), verts_in_cut_info.send_procs_ptr[i],
			LB_VERT_TYPE, i, tag, MPI_COMM_WORLD, (send_request+i));
	  
	  
	}
	// send necessary info here...
	if(verts_in_cut_info.recv_procs_ptr[i] != 0) {
	  j = MPI_Isend(&(verts_in_cut_info.recv_sfc_vert[recv_count]), verts_in_cut_info.recv_procs_ptr[i],
			LB_VERT_TYPE, i, tag, MPI_COMM_WORLD, (recv_request+i));
	}
      }
      recv_count += verts_in_cut_info.recv_procs_ptr[i];
    }
    // wait until the info is sent and received...
    for(i=0;i<numprocs;i++)
      if(i!= myid)
	{
	  if(verts_in_cut_info.send_procs_ptr[i] != 0) {
	    MPI_Status status;
	    j = MPI_Wait(&(send_request[i]), &status);
	  }
	  if(verts_in_cut_info.recv_procs_ptr[i] != 0) {
	    MPI_Status status;
	    j = MPI_Wait(&(recv_request[i]), &status);
	  }
	}
	  
    recv_count = 0;
    for(i=0;i<myid;i++)
      recv_count += verts_in_cut_info.recv_procs_ptr[i];
    for(i=0;i<num_local_objects;i++)
      if(sfc_vert_ptr[i].cut_bin_flag == BSFC_CUT) {
	if(sfc_vert_ptr[i].sfc_key[0] == (unsigned) 3157082553)
	  j = 0;
	if(sfc_vert_ptr[i].destination_proc != myid) {
	  j = sfc_vert_ptr[i].destination_proc;
	  sfc_vert_ptr[i].destination_proc = 
	    send_sfc_vert[proc_counter[sfc_vert_ptr[i].destination_proc]].destination_proc;
	  proc_counter[j] += 1;
	}
	else { // if i need to send to myself... 
	  sfc_vert_ptr[i].destination_proc = verts_in_cut_info.recv_sfc_vert[recv_count].destination_proc;
	  recv_count++;
	}
      }
    delete []proc_counter;
    delete []send_request;
    delete []recv_request;
    delete []send_sfc_vert;
  }
  free(global_actual_work_allocated);
  free(work_percent_array);  

  BSFC_update_element_proc(myid, numprocs, HT_Elem_Ptr, HT_Node_Ptr, sfc_vert_ptr);
    
  //debug stuff...
/*  unsigned* max = new unsigned[numprocs];
  for(i=0;i<numprocs;i++)
    max[i] = 0;
  unsigned* min = new unsigned[numprocs];
  for(i=0;i<numprocs;i++)
    min[i] = ~0;
  unsigned* ustore = new unsigned[numprocs];
  for(i=0;i<num_local_objects;i++) {
    if(max[sfc_vert_ptr[i].destination_proc] < sfc_vert_ptr[i].sfc_key[0])
      max[sfc_vert_ptr[i].destination_proc] = sfc_vert_ptr[i].sfc_key[0];
    if(min[sfc_vert_ptr[i].destination_proc] > sfc_vert_ptr[i].sfc_key[0])
      min[sfc_vert_ptr[i].destination_proc] = sfc_vert_ptr[i].sfc_key[0];
  }
  i = MPI_Allreduce(min, ustore, numprocs, MPI_UNSIGNED, MPI_MIN, MPI_COMM_WORLD);
  if(myid == 0) {
   // printf(" time step %d ****************** the minimums are: ", time_step);
    for(i=0;i<numprocs;i++)
      // printf("%u ",ustore[i]);
    //printf("\n");
  }
  delete []min;
  i = MPI_Allreduce(max, ustore, numprocs, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
  if(myid == 0) {
    //printf(" time step %d ****************** the maximums are: ", time_step);
    for(i=0;i<numprocs;i++)
     // printf("%u ",ustore[i]);
    //printf("\n");
  }
  delete []max;
  delete []ustore;  */
  

  // more debugging stuff...
/*  double* proc_load = new double[numprocs];
  for(i=0;i<numprocs;i++)
    proc_load[i] = 0;
  //meshplotter(HT_Elem_Ptr, HT_Node_Ptr, time_step+10000);
  for(i=0;i<no_of_buckets;i++)
    {
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
	{
	  EmTemp = (Element*)(entryp->value);
	  if(!EmTemp->get_refined_flag()) {
	    if(EmTemp->get_new_old() > 0)
	      proc_load[EmTemp->get_myprocess()] += EmTemp->get_lb_weight();
	  }
	  entryp = entryp->next;
	}
    }
  double* procload2 = new double[numprocs];
  i = MPI_Allreduce(proc_load, procload2, numprocs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if(myid == 0)
    for(i=0;i<numprocs;i++)
     // printf("proc %d load is %e total weight is %e\n",i, procload2[i], total_weight);
  delete []procload2;
  delete []proc_load;*/ 
  //done debug stuff
  free(sfc_vert_ptr);

  BSFC_update_and_send_elements(myid, numprocs, HT_Elem_Ptr, HT_Node_Ptr, time_step);

  
  return;
}

/* create info before starting the multi-level refinement of the bins 
 */

void BSFC_create_refinement_info(int* number_of_cuts, 
				 float* global_actual_work_allocated,
				 float total_weight, 
				 float* work_percent_array,
				 unstructured_communication verts_in_cut_info,
				 float** work_prev_allocated_ptr,
				 int myid, int numprocs)
{
  float my_work, work2;
  int i = 0, j;
  //printf("the total weight is %e\n",total_weight);
  /* find out how many cuts are in this bin. */
  my_work = global_actual_work_allocated[myid];
  while(my_work > total_weight * work_percent_array[myid-i]) {
    i++;
  }
  *number_of_cuts = i;
  if(verts_in_cut_info.recv_count == 0 || myid ==0) {
    *number_of_cuts = 0;
    return;
  }

  /* create link list for objects in the array.  link list
     is set up so that the objects in the array are
     traversed consecutively */
  for(i=0;i<(verts_in_cut_info.recv_count-1);i++)
    verts_in_cut_info.recv_sfc_vert[i].next_sfc_vert_index = i+1;
  verts_in_cut_info.recv_sfc_vert[verts_in_cut_info.recv_count-1].next_sfc_vert_index = -1;
     
  /* update work previously allocated to include work in all bins
     with higher keys than this bin */
  work2 = 0;
  for(i=0;i<verts_in_cut_info.recv_count;i++) 
    work2 += verts_in_cut_info.recv_sfc_vert[i].lb_weight;
  
  work2 = global_actual_work_allocated[myid] - work2;
  
  // note that work_prev_allocated data will only be correct for the
  // actual processors that are getting examined
  *work_prev_allocated_ptr = 
    (float*) malloc(sizeof(float)  * numprocs);
  for(i=myid-*number_of_cuts+1;i<=myid;i++)
    *(*work_prev_allocated_ptr+i) = work2;
  
  return;
}
