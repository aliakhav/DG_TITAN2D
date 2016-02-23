#include"../header/hpfem.h"
#include<time.h>

int howmanyelements(HashTable* HT_Elem_Ptr)
{

  double errorsq;
  double solsq;
  Element* EmTemp;
  HashEntryPtr entryp;
  int myid;
  int counter=0;
  int i;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
 
  int buckets=HT_Elem_Ptr->get_no_of_buckets(); 
  for(i=0; i<buckets; i++)
    {
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
        {
          EmTemp = (Element*)(entryp->value);
          if(!EmTemp->get_refined_flag())
            {
                counter++;
            }
          entryp = entryp->next;
        }
    }
 return counter;
}
