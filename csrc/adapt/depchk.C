#include "../header/hpfem.h"



void depchk(Element* EmTemp, HashTable* HT_Elem_Ptr, int* ifg, Element* refined[], int* count)

  /*---
    refined[] stores the address of ready-for-refinement element of the sub-domain
    refined_temp[] stores the address of ready-for-refinement element triggered by one element refinement
    count is counting the number of refinement of the subdomain
    j is counting the number of refinement triggered by one element refinement
    ---------------*/
{

  int i, j, k;
  Element* element;
  Element* Neigh;
  Element* refined_temp[128];
  void* p;
  int myid, numprocs;
  unsigned send_buf[4*KEYLENGTH];
  unsigned recv_buf[4*KEYLENGTH];

  MPI_Status     status;
  MPI_Request    request;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  for(j=1;j<128;j++) refined_temp[j] = NULL;
  refined_temp[0] = EmTemp;
 
  j = 0; k = 0;
  element = EmTemp;//-- EmTemp is the trigger of this round of refinement

  while(element&&(j<10))//--element is temporary varible
    {
      for(i=0;i<4;i++)//-- checking the four neighbors to identify which must be refined
	{
	  int neigh_proc = *(element->get_neigh_proc()+i);
	  
	  if((neigh_proc!=-1)&&(neigh_proc!=-2))//-- if there is a neighbor
	    {
	      
	      Neigh = (Element*)(HT_Elem_Ptr->lookup(element->get_neighbors()+i*KEYLENGTH));
	      
	      if(Neigh != NULL && neigh_proc == myid) //-- if this neighbor is in the same proc as element is
		{
		  if((!Neigh->get_refined_flag())&&element->get_gen()>Neigh->get_gen())//-- if the neighbor is bigger, then it must be refined
		    {
		      int flag = 1; int m = 0;
		      while(refined_temp[m]) 
			{ 
			  if(*refined_temp[m]->pass_key() == *Neigh->pass_key()) 
			    {flag = 0; break;}
			  else m++;
			}
		      
		      if(flag)//-- if this neighbor has not yet been marked
			{
			  j++; 
			  refined_temp[j] = Neigh;
			}
		    }
		}
	      else//-- need neighbor's generation infomation
		{
		  if(element->get_gen()>*(element->get_neigh_gen()+i))//--stop this round of refinement
		    {
		      *ifg = 0;
		      for(int m=0;m<128;m++) refined_temp[m] = NULL;
		      break;
		    }
		}
	    }			    
			    
	}
      if(!*ifg) break;
      k++;
      element = refined_temp[k];//--check next
    }
  
  if(*ifg)
    {
      if(j<10)//-- 10 is the maximum tolerence of related refinement
	{
	  int m = 0;
	  while(refined_temp[m])
	    {
	      int sur = 0; int mi = 0;
	      while(refined[mi])
		{
		  if(*(refined[mi]->pass_key()) == *(refined_temp[m]->pass_key())) //-- KEYLENGTH should be considered
		    {
		      sur = 1;
		      break;
		    }
		  mi++;
		}
	      if(!sur)
		{
		  refined[*count] = refined_temp[m];
		  *count = *count+1;
		}
	      m++;
	    }
	}
      else
	{
	  *ifg = 0;//-- refuse to do the refinement
	}
    }
}



