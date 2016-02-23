#include "../header/hpfem.h"

void rm_surplus_dof(NNLink* NNHead, SSLink* SSHead, HashTable* ht_elem_ptr, 
		    HashTable* ht_node_ptr, int myid, int numprocs, int* nsveb)
{

  Element*       EmTemp;
  Node*          NdTemp;
  unsigned       KeyTemp[KEYLENGTH];
  HashEntryPtr   entryp;
  unsigned*      keyP;
  void*          p;
  int*           assocP;
  int*           dofP;
  int            power;
  int            NN_total;

  NNLink*  NN_old;
  NNLink*  NN_new;
  SSLink*  SS_old;
  SSLink*  SS_new;

  int  Npre;

  int i, j, k;
  int on_off;

  MPI_Status    status;
  MPI_Request   request;

  int  NNtype = 1;
  int  SStype = 2;
  int  VVtype = 3;
  int  EEtype = 4;
  int  BBtype = 5;
  int  nsvebtype = 6;
  int  keytype = 7;
  int  rmtype = 8;

  unsigned k_buf[KEYLENGTH+1];//--the last digit is used as a counter
  int rm = -1;//--indicating rm numbers

  if(myid == 0)
    {
      for(i=1;i<numprocs-1;i++)//-- proc numproc-1 need not to send out message
	{
	  j = 0;//--counter, indicating how many N of proc i be numbered twice 
	  on_off = ON;
	  while(on_off)
	    {
	      MPI_Irecv(&k_buf, KEYLENGTH, MPI_UNSIGNED, i, keytype, MPI_COMM_WORLD, &request);
	      MPI_Wait(&request, &status);
	      if((*k_buf)||(*(k_buf+1)))//--if the communication is not completed 
		{
		  NN_old = NNHead->next;//---scan the NNLink 
		  while(NN_old)
		    {
		      //if(*(NN_old->key) == *k_buf) 
		      if ( compare_key ( NN_old->key, k_buf ) )
			{
			  //j++;
			  j += EQUATIONS;
			  break;
			}
		      else
			{
			  if(!NN_old->next)//--if at the end of the NNLink, insert the new one to the tail
			    {
			      NN_new = new NNLink(k_buf);
			      NN_old->next = NN_new;
			      NN_new->pre  = NN_old;
			      break;
			    }
			  else NN_old = NN_old->next;//--else, compared with the next one 
			}
		    }
		}
	      else 
		{
		  MPI_Isend(&j, 1, MPI_INT, i, rmtype, MPI_COMM_WORLD, &request);
		  MPI_Wait(&request, &status);
		  on_off = OFF;//---communication with proc i completed 
		}
	    }
	}

      NN_old = NNHead->next; i = -1; 
      while(NN_old)
	{ 
	  //i++;                 
	  i += EQUATIONS;
	  p = ht_node_ptr->lookup(NN_old->key); 
	  
	  for(j=0;j<KEYLENGTH;j++) k_buf[j] = *(NN_old->key+j);//| prepare sending buf 
	  k_buf[KEYLENGTH] = i;                                //| last is global ordering 
	  for(j=1;j<numprocs;j++)
	    {//---share NN ordering with other proc 
	      MPI_Isend(&k_buf, KEYLENGTH+1, MPI_UNSIGNED, j, keytype, MPI_COMM_WORLD, &request);
	      MPI_Wait(&request, &status);
	    }
	  if(p) 
	    {
	      NdTemp = (Node*)p;
	      NdTemp->putglnum(i);//--put NN 
	    }
	  NN_old=NN_old->next; 
	  //cout <<"\n"<<*k_buf<<"  "<<i<<"\n"; 
	} 
      NN_total = i;         //--total num of NN 
      k_buf[KEYLENGTH] = 0; //--flag for indicating sending was finished 
      for(j=1;j<numprocs;j++)
	{
	  MPI_Isend(&k_buf, KEYLENGTH+1, MPI_UNSIGNED, j, keytype, MPI_COMM_WORLD, &request); 
	  MPI_Wait(&request, &status);//---NN sending finished 
	}
    }
  
  else //--proc other than 0 
    { 

      NN_old = NNHead->next;//--preparing for sending NN link to proc 0 
      while(NN_old) 
	{ 
	  for(j=0;j<KEYLENGTH;j++) k_buf[j] = *((NN_old)->key + j);//--key is written to buffer 
	  MPI_Isend(&k_buf, KEYLENGTH, MPI_UNSIGNED, 0, keytype, MPI_COMM_WORLD, &request); 
	  MPI_Wait(&request, &status);
	  NN_old = NN_old->next;
	}
      *k_buf=0; *(k_buf+1)=0;//--indicate proc 0 that NN link sending finished 

      if(myid != numprocs -1)
	{
	  MPI_Isend(&k_buf, KEYLENGTH, MPI_UNSIGNED, 0, keytype, MPI_COMM_WORLD, &request);
	  MPI_Wait(&request, &status);
	  for(;;)
	    {
	      MPI_Irecv(&rm, 1, MPI_INT, 0, rmtype, MPI_COMM_WORLD, &request); 
	      MPI_Wait(&request, &status);
	      if(rm !=-1) break;
	    }
	  nsveb[1] =nsveb[1]  - rm;//--update local NN number 
	}

      k_buf[KEYLENGTH] = NONZERO;//---prepare for recv from proc 0 
      while(k_buf[KEYLENGTH])
	{   //--continue to recv info from proc 0
	  MPI_Irecv(&k_buf, KEYLENGTH+1, MPI_UNSIGNED, 0, keytype, MPI_COMM_WORLD, &request);  
	  MPI_Wait(&request, &status);
	  p = ht_node_ptr->lookup(k_buf);
	  if(p)
	    {
	      NdTemp = (Node*)p;
	      if(k_buf[KEYLENGTH]) NdTemp->putglnum(k_buf[KEYLENGTH]);
	    }
	}

    }
  //--------------------supplus numbering elimination fineshed------------
  
}
