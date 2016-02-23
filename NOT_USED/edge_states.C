#include "../header/hpfem.h"

void calc_edge_states(HashTable* El_Table, HashTable* NodeTable,
		      int myid)
{
  int i,j,k, counter;
  double tiny = GEOFLOW_TINY;
  int el_counter = 0;
  double evalue = 1;

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
	    if(r_flag == 0 )//if this is a refined element don't involve!!!
	      {
		double pheight = *(Curr_El->get_state_vars());
		Curr_El->calc_edge_states(El_Table, NodeTable, myid);
		double pheight2 = *(Curr_El->get_state_vars());
		if(pheight != pheight2)
		    printf("prolbem of changing height here,,,.....\n");
		el_counter++;
	      }
	    currentPtr=currentPtr->next;      	    
	  }
      }

  return;
}
