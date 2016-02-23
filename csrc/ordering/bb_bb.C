#include "../header/hpfem.h"

void bb_bb(HashTable* ht_node_ptr, HashTable* ht_elem_ptr, int* nsveb,
	   BBLink* BBHead, int* DofCounter)
{

  int        i, j, k;
  Node*      NdTemp;
  int        order;
  unsigned*  keyP; 
  Element*   EmTemp;
  void*      p;
  
  HashEntryPtr  entryp;

  BBLink*   BB_new;

  BBLink*  BBTail = BBHead;
  while(BBTail->next)
      BBTail = BBTail->next;

  BBLink*  BB_old = BBTail;
/*

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if(myid==2)
        {

	   ofstream dout("hash2.dat", ios::out);
           HashEntryPtr  entryp2;
           Element* dent;
           entryp2 = *(ht_elem_ptr->getbucketptr() + 499);
           while(entryp2)
            {
             dent=(Element*)entryp2->value;
             dout<<entryp2->key[0]<<" "<<entryp2->key[1]<<setw(30);
             dout<<*(dent->pass_key())<<" "<<*(dent->pass_key()+1)<<endl<<flush;
             entryp2=entryp2->next;
            }
        }
*/
  
  int no_of_buckets = ht_elem_ptr->get_no_of_buckets();

  for(i=0;i<no_of_buckets;i++)
    {
      entryp = *(ht_elem_ptr->getbucketptr() + i);
      while(entryp)
	{
	  EmTemp = (Element*)(entryp->value);
	  if(!EmTemp->get_refined_flag())
	    {
	      keyP   = EmTemp->pass_key();
	      p      = ht_node_ptr->lookup(keyP);
	      NdTemp = (Node*)p;
	      NdTemp->putinfo(BUBBLE);
	      order = *(EmTemp->get_order()+4);
	      NdTemp->put_order(order);
	      BB_new = new BBLink(keyP);
	      BB_old->next = BB_new;  
	      BB_new->pre  = BB_old;
	      BB_old = BB_new;  
	      *(nsveb+5) += pow((order-1), 2)*EQUATIONS;
	      NdTemp->putdof(*DofCounter, *DofCounter+pow((order-1),2)*EQUATIONS-1);
	      *DofCounter = *DofCounter+pow((order-1),2)*EQUATIONS;
	    }
	  entryp = entryp->next;
	}
    }
}
