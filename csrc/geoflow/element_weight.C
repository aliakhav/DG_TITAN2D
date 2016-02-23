#include "../header/hpfem.h"

double element_weight(HashTable* El_Table, HashTable* NodeTable, int myid, int nump) {
  int i,j,k, counter;
  double tiny = GEOFLOW_TINY;
  int el_counter = 0;
  double evalue = 1;
  double sub_weight[2] = {0,0}; // second number is to keep track of the number of objects
  
  //-------------------go through all the elements of the subdomain and  
  //-------------------find the edge states

  HashEntryPtr* buck = El_Table->getbucketptr();
  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i))
      {
	HashEntryPtr currentPtr = *(buck+i);
	while(currentPtr)
	  {
	    Element* Curr_El=(Element*)(currentPtr->value);
	    int r_flag = Curr_El->get_refined_flag(); 
	    if(r_flag == 0 )  {  //if this is a refined element don't involve!!! 
	      sub_weight[0] += *(Curr_El->get_el_error());
	      sub_weight[1] += 1;
	      
	    }
	    
	    currentPtr=currentPtr->next;      	    
	  }
      }

  double global_weight[2];
  i = MPI_Allreduce(sub_weight, global_weight, 2,
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  global_weight[0] = 1.5*(global_weight[0]-global_weight[1]*WEIGHT_ADJUSTER)/global_weight[1]+GEOFLOW_TINY;

  return global_weight[0];
}
