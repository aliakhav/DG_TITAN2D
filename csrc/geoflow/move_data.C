#include "../header/hpfem.h"
#include "../header/exvar.h"
#include "../header/geoflow.h"
//#define DEBUG
//#define PRINT_MOVE

void move_data(int nump, int myid, HashTable* El_Table, HashTable* NodeTable)
{

#ifdef PRINT_MOVE
   MPI_Barrier(MPI_COMM_WORLD);
  if(myid == 0)    printf("========================================================================================\n");
#endif
  int i, j, k, ierr;
  /* assume that no elements share a neighboring element on another processor */
  int* num_send_recv = new int[nump];
  for(i=0;i<nump;i++)
    num_send_recv[i] = 0;

  /* count how many elements we should send and receive from other procs */
  HashEntryPtr* buck = El_Table->getbucketptr();
  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i))
      {
	HashEntryPtr currentPtr = *(buck+i);
	while(currentPtr)
	  {
	    Element* Curr_El=(Element*)(currentPtr->value);
#ifdef DEBUG	   
	    printf("lbweight %e",Curr_El->get_lb_weight());
	    assert(Curr_El->get_lb_weight()< 1000000000000000);
#endif	    
 	    if(!(Curr_El->get_refined_flag()))//if this is a refined element don't involve!!!
	      {
		int* neigh_proc = Curr_El->get_neigh_proc();
		assert(neigh_proc);
		for(j=0;j<8;j++) {
		  assert(neigh_proc[j]>-100);
		  if(neigh_proc[j] >= 0 ) {
		    num_send_recv[neigh_proc[j]] += 1;
		  }
		}
		      
	      }
	    currentPtr=currentPtr->next;      	    
	  }
      }

  num_send_recv[myid] = 0;  // don't need to send info to myself

  ElemPack** send_array = new ElemPack*[nump];
  ElemPack** recv_array = new ElemPack*[nump];
  for(i=0;i<nump;i++)
    if(num_send_recv[i] != 0) {
      send_array[i] = new ElemPack[num_send_recv[i]];
      recv_array[i] = new ElemPack[num_send_recv[i]];      
    }


  int* proc_counter = new int[nump];
  for(i=0;i<nump;i++)
    proc_counter[i] = 0;

  int juuukk = 0;
  /* put (GHOST) elements to be moved in the proper arrays */
  for(i=0;i<El_Table->get_no_of_buckets();i++)
    if(*(buck+i))
      {
	HashEntryPtr currentPtr = *(buck+i);
	while(currentPtr)
	  {
	    Element* Curr_El=(Element*)(currentPtr->value);
#ifdef DEBUG	   
	    printf("lbweight %e",Curr_El->get_lb_weight());
	    assert(Curr_El->get_lb_weight()< 1000000000000000);
#endif	    
 	    if(!(Curr_El->get_refined_flag()))//if this is a refined element don't involve!!!
	      {
		int* neigh_proc = Curr_El->get_neigh_proc();
		for(j=0;j<8;j++) 
		  if(neigh_proc[j] != myid && neigh_proc[j] >= 0) {
		    Pack_element(Curr_El, 
		     (send_array[neigh_proc[j]]+proc_counter[neigh_proc[j]]), 
				 NodeTable, 
				 neigh_proc[j]);
		    (send_array[neigh_proc[j]]+
		     proc_counter[neigh_proc[j]])->refined = GHOST;
		    proc_counter[neigh_proc[j]] += 1;
		  }
		}
	    currentPtr=currentPtr->next;      	    
	  }
      }  

  MPI_Request* request = new MPI_Request[2*nump]; //first nump is for send, 2nd nump is for recv
  
#ifdef PRINT_MOVE
  for(i=0;i<nump;i++)
    printf("proc %d is sending/receiving %d from %d\n", myid, num_send_recv[i], i);
#endif
  int send_tag = 22674;

  for(i=0;i<nump;i++)
   if(num_send_recv[i] != 0) {
     j = MPI_Isend((void*) send_array[i], num_send_recv[i], 
		   ELEMTYPE, i, send_tag+myid, MPI_COMM_WORLD, (request+i));
     j = MPI_Irecv((void*) recv_array[i], num_send_recv[i],  
		   ELEMTYPE, i, send_tag+i, MPI_COMM_WORLD, (request+nump+i));
   }
  
  int add_counter = 0, update_counter = 0;
  MPI_Status status;
  for(i=0;i<nump;i++)
    if(num_send_recv[i] != 0) {
      ierr = MPI_Wait((request+nump+i), &status);

      for(j=0;j<num_send_recv[i];j++) {
	Element* elm = (Element*) (El_Table->lookup((recv_array[i]+j)->key));
	if(elm == NULL) { // this elm doesn't exist on this proc
	  Element* new_elm = new Element();
	  construct_el(new_elm, (recv_array[i]+j), NodeTable, myid);
	  El_Table->add(new_elm->pass_key(), new_elm);
	  add_counter++;
	}
	else { //this elm is on this proc - copy state variables and slopes to the existing elm
	  double* dPtr = elm->get_el_solution();
	  double* d2Ptr = (recv_array[i]+j)->el_solution;
	  for(k=0;k<ELM_DOF;k++){
	    dPtr[k] = d2Ptr[k];
	    //   if (d2Ptr[0]>0)
	    //   printf("el_solution %e",d2Ptr[k]);
	  }
	  
	  dPtr = elm->get_prev_el_solution();
	  d2Ptr = (recv_array[i]+j)->prev_el_solution;
#ifdef DEBUG	   
	    printf("lbweight %e",(recv_array[i]+j)->lb_weight);
	    assert((recv_array[i]+j)->lb_weight< 10000000000000);
#endif

	  for(k=0;k<ELM_DOF;k++){ 
	    dPtr[k] = d2Ptr[k];
	    //   if (d2Ptr[0]>0)
	    //   printf("el_solution %e",d2Ptr[k]);
	  }
	  update_counter++;
	}
      }
    }
#ifdef PRINT_MOVE
  printf("proc %d has added %d  and updated %d ghost elements \n", myid, add_counter, update_counter);
#endif
  for(i=0;i<nump;i++)
    if(num_send_recv[i] != 0) 
      ierr = MPI_Wait((request+i), &status);
  
  delete []request;
  for(i=0;i<nump;i++)
    if(num_send_recv[i] != 0) {
      delete [](send_array[i]);
      delete [](recv_array[i]);
    }
  delete []send_array;
  delete []recv_array;
  delete []num_send_recv;
  delete []proc_counter;
      
  return;
}


/* delete the ghost elements that were put in the element hashtable */
void delete_ghost_elms(HashTable* El_Table, int myid) 
{
  int i;
  int delete_counter = 0;
  HashEntryPtr* buck = El_Table->getbucketptr();
  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i))
      {
	HashEntryPtr currentPtr = *(buck+i);
	while(currentPtr)
	  {
	    Element* Curr_El=(Element*)(currentPtr->value);
	    
 	    if(Curr_El->get_refined_flag() == GHOST) {
	      HashEntryPtr deletePtr = currentPtr;
	      currentPtr = currentPtr->next;
	      El_Table->remove(Curr_El->pass_key());
	      delete Curr_El;
	      delete_counter++;
	    }
	    else
	      currentPtr=currentPtr->next;      	    
	  }
      }
#ifdef PRINT_MOVE
  printf("proc %d has deleted %d ghost elms \n",myid, delete_counter);
#endif
  return;
}

void create_delete_memory() {

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  
  double* junk = new double[100000];
  junk[0] = 0;
  delete []junk;

  return;
}
