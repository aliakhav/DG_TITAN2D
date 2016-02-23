#include"../header/hpfem.h"
#include<time.h>

double random_error(HashTable* HT_Elem_Ptr)
{

  double errorsq;
  double solsq;
  Element* EmTemp;
  HashEntryPtr entryp;
  int myid;
  int counter=0;
  int i;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
 
  
  for(i=0; i<HT_Elem_Ptr->get_no_of_buckets(); i++)
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




	double TARGET=0;

  for(i=0; i<HT_Elem_Ptr->get_no_of_buckets(); i++)
    {
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
	{
	  EmTemp = (Element*)(entryp->value);
	  if(!EmTemp->get_refined_flag())
	    {
	      //srand(time(0)+myid*1000+i*100);
	      srand(myid*1000+i*100);
		double help=rand()%999;
	      errorsq=pow(help/1000, 2)/counter;
	      TARGET+=sqrt(errorsq)/counter;
	      solsq=0;
	      EmTemp->putel_sq(solsq, errorsq); 
//		cout<<myid<<"  *************the random number: "<<help<<endl<<flush; 
//		cout<<myid<<" element error squareroot: "<<sqrt(errorsq)<<endl<<flush;
	    }
	  entryp = entryp->next;
	}
    }
  return TARGET;
  
}
