#include "../header/hpfem.h"

void calc_average(HashTable* El_Table, HashTable* NodeTable)
{
  int i;

  HashEntryPtr* buck = El_Table->getbucketptr();
  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i))
      {
	HashEntryPtr currentPtr = *(buck+i);
	while(currentPtr)
	  {
	    Element* Curr_El=(Element*)(currentPtr->value);
 	    if(Curr_El->get_refined_flag() == 0)
	      Curr_El->calc_average(NodeTable);
	    currentPtr=currentPtr->next;      	    
	  }
      }

  return;
}

