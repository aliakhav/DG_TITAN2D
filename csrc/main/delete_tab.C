
#include "../header/hpfem.h"

void Delete_Table(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr)
{

  int i, j, k;
  HashEntryPtr entryp;

  int elements = HT_Elem_Ptr->get_no_of_buckets();
  int nodes    = HT_Node_Ptr->get_no_of_buckets();   
  for(i=0;i<elements;i++){ 
    entryp = *(HT_Elem_Ptr->getbucketptr() + i); 
    while(entryp)
      { 
	Element* EmTemp = (Element*)(entryp->value);
	delete EmTemp;
	entryp = entryp->next;
      }
  }

  for(i=0;i<nodes;i++){ 
    entryp = *(HT_Node_Ptr->getbucketptr() + i); 
    while(entryp)
      { 
	Node* NdTemp = (Node*)(entryp->value);
	delete NdTemp;
	entryp = entryp->next;
      }
  }

  delete HT_Elem_Ptr;
  delete HT_Node_Ptr;

}
