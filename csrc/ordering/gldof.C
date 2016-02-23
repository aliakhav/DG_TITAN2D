#include "../header/hpfem.h"
 
void gldof(HashTable* ht_elem_ptr, HashTable* ht_node_ptr, int* nsveb, 
	   int myid, int numprocs, NNLink* NNHead, SSLink* SSHead, 
	   VVLink* VVHead, EELink* EEHead, BBLink* BBHead,
	   int* Nn, int* Ss, int* Vv, int* Ee, int* Bb)
{
  int            i, j, k, mi;
  InterFacePtr   InTemp;
  Element*       EmTemp;
  Node*          NdTemp;
  unsigned       KeyTemp[KEYLENGTH];
  HashEntryPtr   entryp;
  unsigned*      keyP;
  void*          p;
  int*           assocP;
  int*           dofP;
  int            interf;
  int            NN_total;
  int            order;

  int  on_off;
  int  Npre;

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

  NNPtr  NN_new;
  NNPtr  NN_old;
  SSPtr  SS_new;
  SSPtr  SS_old;
  VVPtr  VV_new;
  VVPtr  VV_old;
  EEPtr  EE_new;
  EEPtr  EE_old;
  BBPtr  BB_new;
  BBPtr  BB_old;

  unsigned k_buf[KEYLENGTH+1];//--the last digit is used as a counter
  int rcv_buf[64][6];
  if(numprocs > 64)
    printf("change rcv_buffer size in gldof.C\n");
  assert(numprocs < 65);
  
  int rm = -1;//--indicating rm numbers 

 for(i=0;i<myid;i++)
    {
      MPI_Isend(nsveb, 6, MPI_INT, i, nsvebtype, MPI_COMM_WORLD, &request);
      MPI_Wait(&request, &status); 
    }
  for(i=myid+1;i<numprocs;i++)
    {
      MPI_Isend(nsveb, 6, MPI_INT, i, nsvebtype, MPI_COMM_WORLD, &request);
      MPI_Wait(&request, &status); 
    }
  for(i=0;i<myid;i++)
    {
      MPI_Irecv(&rcv_buf[i][0], 6, MPI_INT, i, nsvebtype, MPI_COMM_WORLD, &request);
      MPI_Wait(&request, &status);
    } 
  for(i=myid+1;i<numprocs;i++)
    {
      MPI_Irecv(&rcv_buf[i][0], 6, MPI_INT, i, nsvebtype, MPI_COMM_WORLD, &request);
      MPI_Wait(&request, &status);
    } 

  //-----NN------------ 
  //--------global ordering has been done in proc 0 
  *Nn = 0;
  for(i=0;i<myid;i++)
    *Nn = *Nn+rcv_buf[i][1];
  for(i=myid+1;i<numprocs;i++)
    *Nn = *Nn+rcv_buf[i][1];
  *Nn = *Nn+nsveb[1];//--total NN->*Nn  

  //------SS--------- 
  *Ss = 0;
  for(i=0;i<myid;i++)
    *Ss = *Ss + rcv_buf[i][2];//----caculate the beginning position for this proc SS. *Ss is the beginning--- 

  SS_old = SSHead->next;
  i = 1;
  while(SS_old)
    {
      p = ht_node_ptr->lookup(SS_old->key);
      NdTemp = (Node*)p;
      order  = NdTemp->get_order();

      for(j=0;j<KEYLENGTH;j++) k_buf[j] = *(SS_old->key+j);
      k_buf[KEYLENGTH] = i+*Ss+*Nn;                          
      for(j=1;j<myid;j++)
	{
	  MPI_Isend(k_buf, KEYLENGTH+1, MPI_UNSIGNED, j, keytype, MPI_COMM_WORLD, &request);
	  MPI_Wait(&request, &status);
	}
      for(j=myid+1;j<numprocs;j++)
	{
	  MPI_Isend(k_buf, KEYLENGTH+1, MPI_UNSIGNED, j, keytype, MPI_COMM_WORLD, &request);
	  MPI_Wait(&request, &status);
	}
      NdTemp->putglnum(i+*Ss+*Nn);
      SS_old = SS_old->next;
      i += (order-1)*EQUATIONS;
    }
  
  k_buf[KEYLENGTH] = 0;
  
  for(j=1;j<myid;j++)
    {
      MPI_Isend(k_buf, KEYLENGTH+1, MPI_UNSIGNED, j, keytype, MPI_COMM_WORLD, &request);
      MPI_Wait(&request, &status);
    }
  for(j=myid+1;j<numprocs;j++)
    {
      MPI_Isend(k_buf, KEYLENGTH+1, MPI_UNSIGNED, j, keytype, MPI_COMM_WORLD, &request);
      MPI_Wait(&request, &status);
    }//-----------------------sending finished------------- 
  
  if(myid != 0)
    {
      for(j=0;j<myid;j++)
	{
	  k_buf[KEYLENGTH] = NONZERO;//---prepare for recv from proc j 
	  while(k_buf[KEYLENGTH])
	    {//--continue to recv info from proc j
	      MPI_Irecv(k_buf, KEYLENGTH+1, MPI_UNSIGNED, j, keytype, MPI_COMM_WORLD, &request); 
	      MPI_Wait(&request, &status);
	      if(k_buf[KEYLENGTH]) 
		{
		  p = ht_node_ptr->lookup(k_buf);
		  //if(!p) {cout<<"Aborting SS "<<myid<<"\n"; exit(0);}
		  //assert(p);
		  if(p)
		    {
		      NdTemp = (Node*)p;
		      NdTemp->putglnum(k_buf[KEYLENGTH]);
		    }
		}
	    }
	}
      for(j=myid+1;j<numprocs;j++)
	{
	  k_buf[KEYLENGTH] = NONZERO;
	  while(k_buf[KEYLENGTH])
	    {
	      MPI_Irecv(k_buf, KEYLENGTH+1, MPI_UNSIGNED, j, keytype, MPI_COMM_WORLD, &request);  
	      MPI_Wait(&request, &status);
	      if(k_buf[KEYLENGTH])
		{
		  p = ht_node_ptr->lookup(k_buf);
		  //if(!p) {cout<<"Aborting SS "<<myid<<"\n"; exit(0);}
		  //assert(p);
		  if(p)
		    {
		      NdTemp = (Node*)p;
		      NdTemp->putglnum(k_buf[KEYLENGTH]);  
		    }
		}
	    }
	}
    }
  *Ss = 0;//--counting
  for(i=0;i<myid;i++)
    *Ss = *Ss+rcv_buf[i][2];
  for(i=myid+1;i<numprocs;i++)
    *Ss = *Ss+rcv_buf[i][2];
  *Ss = *Ss+nsveb[2];//---total SS->*Ss  

  //------VV---------- 
  
  for(i=0;i<myid;i++)
    *Vv = *Vv + rcv_buf[i][3];//----caculate the beginning position for this proc VV--- 
  VV_old = VVHead->next;
  i = 1;
  while(VV_old)
    {
      p = ht_node_ptr->lookup(VV_old->key);
      NdTemp = (Node*)p;
      NdTemp->putglnum(i+*Vv+*Ss+*Nn);
      VV_old = VV_old->next;
      i += EQUATIONS;
    }

  *Vv = 0;
  for(i=0;i<myid;i++)
    *Vv = *Vv+rcv_buf[i][3];
  for(i=myid+1;i<numprocs;i++)
    *Vv = *Vv+rcv_buf[i][3];
  *Vv = *Vv+nsveb[3];//---total VV->Vv 

  //------EE----------
  
  for(i=0;i<myid;i++)
    *Ee = *Ee + rcv_buf[i][4];//----caculate the beginning position for this proc EE--- 
  EE_old = EEHead->next;
  i = 1;
  while(EE_old)
    {
      p = ht_node_ptr->lookup(EE_old->key); 
      //if(!p) cout<<"Failed"<<"\n";
      NdTemp = (Node*)p;
      order  = NdTemp->get_order();
      NdTemp->putglnum(i+*Ee+*Vv+*Ss+*Nn);
      EE_old = EE_old->next;
      //i++;
      i += (order-1)*EQUATIONS;//---modified in 5.23 
    }
  *Ee = 0;
  for(i=0;i<myid;i++)
    *Ee = *Ee+rcv_buf[i][4];
  for(i=myid+1;i<numprocs;i++)
    *Ee = *Ee+rcv_buf[i][4];
  *Ee = *Ee+nsveb[4];//---total EE->*Ee 


  //------BB---------- 
  
  for(i=0;i<myid;i++)
    *Bb = *Bb + rcv_buf[i][5];//----caculate the beginning position for this proc BB--- 
  BB_old = BBHead->next;
  i = 1;
  while(BB_old)
    {
      p = ht_node_ptr->lookup(BB_old->key);
      NdTemp = (Node*)p;
      order = NdTemp->get_order();
      NdTemp->putglnum(i+*Bb+*Ee+*Vv+*Ss+*Nn);
      BB_old = BB_old->next;
      //i++; 
      i += pow(order-1, 2)*EQUATIONS;
    }
  *Bb = 0;
  for(i=0;i<myid;i++)
    *Bb = *Bb+rcv_buf[i][5];
  for(i=myid+1;i<numprocs;i++)
    *Bb = *Bb+rcv_buf[i][5];
  *Bb = *Bb+nsveb[5];//---total BB->*Bb 

}
