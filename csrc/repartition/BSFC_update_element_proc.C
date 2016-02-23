#include "../header/hpfem.h"
#include "../header/exvar.h"
#include "../header/repartition_BSFC.h"
 
// this routine updates the processor that an element is assigned to

void BSFC_update_element_proc(int myid, int numprocs, 
			      HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
			      BSFC_VERTEX_PTR sfc_vert_ptr)
{
  int i, j, k, no_of_buckets = HT_Elem_Ptr->get_no_of_buckets();
  HashEntryPtr entryp;
  Element* EmTemp;
			    
  j = 0; 
  for(i=0;i<no_of_buckets;i++)
    {
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
	{
	  EmTemp = (Element*)(entryp->value);
	  if(!EmTemp->get_refined_flag()) {
	    if(*(EmTemp->pass_key()) == (unsigned) 74470341)
	      k = i;
	    if(EmTemp->get_new_old() == 1) {  //only this element
	      if(EmTemp->get_myprocess() != sfc_vert_ptr[j].destination_proc) { // this element will get moved to a new processor
		EmTemp->put_myprocess(sfc_vert_ptr[j].destination_proc);
	      }
	      j++;
	    }
	    else if(EmTemp->get_new_old() > 1) {  //multiple elements connected by a constrained node
	      //check for constrained nodes on the vertex nodes
	      k = 4;
	      while(k<8)  {
		int ll = 1;
		Node* ndtemp = (Node*) HT_Node_Ptr->lookup((EmTemp->getNode()+k*KEYLENGTH));
		if(ndtemp->getinfo() == S_S_CON) {
		  BSFC_combine_elements(k-4, EmTemp, HT_Elem_Ptr, HT_Node_Ptr,
					sfc_vert_ptr[j].destination_proc);
		  k = 8;  //exit out of the loop because we found a constrained node...
		}
		k++;
	      }	      
	      j++;
	    }
	  }
	  entryp = entryp->next;
	}
    }

  return;
}

