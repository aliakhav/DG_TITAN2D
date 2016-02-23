#include "../header/hpfem.h"
#include "../header/geoflow.h"

void slopes(HashTable* El_Table, HashTable* NodeTable, double gamma)
{
  int i;
  //-------------------go through all the elements of the subdomain------------------------
  //-------------------and   --------------------------
  
  HashEntryPtr* buck = El_Table->getbucketptr();
  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i))
      {
	HashEntryPtr currentPtr = *(buck+i);
	while(currentPtr)
	  {
	    Element* Curr_El=(Element*)(currentPtr->value);
 	    if(!(Curr_El->get_refined_flag()))//if this is a refined element don't involve!!!
	      Curr_El->get_slopes(El_Table, NodeTable, gamma);
		
	    
	    currentPtr=currentPtr->next;      	    
	  }
      }

  return;
}
