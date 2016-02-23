#include"../header/hpfem.h"
#include"../header/exvar.h"
#include"../header/post.h"


extern void MPI_New_Datatype();

void Pack_Sol(Neigh_Sol* tempsol,Neigh_Sol_Pack* el_Pack);

extern double* elem_sol(double* u, int* ix, HashTable* HT_Node_Ptr, HashTable* HT_Elem_Ptr,int myid, int numprocs,unsigned* key_el, double* Utemp);

extern void get_el_cord(HashTable* ht_node_ptr,Element* EmTemp,double* Xnod);


Neigh_Sol* prep_neigh_list(HashTable *ht_elem_ptr,HashTable *ht_node_ptr,
			   Neigh_Sol* neigh_list,int * elem_num,double* u,int* ix)
{ 
  Neigh_Sol* tempsol;

  int            i,j,k,mi,ijl,l;
  int            flag;
  int myid ,numprocs;
  Element*       EmTemp;
  Node*          NdTemp;
  unsigned       KeyTemp[KEYLENGTH];
  HashEntryPtr   entryp;
  HashEntry**    temp_entryp;
  unsigned*      keyP;
  unsigned*      neighbor;
  void*          p;
  int*           assocP;

  int elements = ht_elem_ptr->get_no_of_buckets(); 
  int nodes    = ht_node_ptr->get_no_of_buckets();
  
  
  FILE *fpt; 
  if((fpt = fopen("debug.txt","w"))==NULL)
  printf("error opening file\n");

  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  double Utemp[242]; //used to store temporary element solution

  //scan the hashtable, 


  for(i=0;i<numprocs;i++)
    elem_num[i] = 0;

  temp_entryp = ht_elem_ptr->getbucketptr();
  for(i=0;i<elements;i++)
    { 
      entryp = *(temp_entryp + i);
      while(entryp)
	{
	  
          EmTemp = (Element*)(entryp->value);
	  
	  if(!EmTemp->get_refined_flag())
	    {
	      assocP    = EmTemp->getassoc();
	      
	      for(mi=0;mi<4;mi++)/*--scan the sides--*/  
		{
		  if((*(assocP+mi) != -2) && (*(assocP+mi) != -1)&&(*(assocP+mi)!=myid))
		    { 
		      
	              tempsol= &neigh_list[assocP[mi]];
                      neighbor =  EmTemp->get_neighbors();
                      flag=0;
		      while(tempsol->next)
			{ 
			  int kk=0;
			  tempsol=tempsol->next;
			  for(l=0;l<KEYLENGTH;l++)
			    {                                     
			      if(tempsol->key[l] ==*(EmTemp->pass_key()+l))kk++;
			    }
			  if(kk==KEYLENGTH){flag=1;break;} //possibl bug
			  
			}
		      if(flag!=1)
			{
			  tempsol->next = new Neigh_Sol;
			  tempsol=tempsol->next;
			  elem_num[assocP[mi]]++; 
			  tempsol->key=EmTemp->pass_key();
			  
			  elem_sol(u,ix,ht_node_ptr,ht_elem_ptr,myid,numprocs,EmTemp->pass_key(),Utemp);  
			  for(int ii=0;ii<2;ii++)
			    for(int jj=0;jj<121;jj++)
			      {
				tempsol->solu[ii][jj] = Utemp[ii*121 +jj];
			      }
			  tempsol->nside =  mi>3 ? mi = mi-4:mi=mi; 
			  tempsol->norder = EmTemp->get_order();
			  tempsol->next =NULL;
			  get_el_cord(ht_node_ptr,EmTemp,tempsol->Xnod);
			}                 
		      
                    }
		  
		}
	      
	      
	    }// if(!EmTemp->get_refined_flag())
	  entryp =entryp->next;
	  
	}//end of while
      
    }//end of for
  
  
  fclose(fpt);  
  return neigh_list;
  
}//end of prep elem list





void send_neigh_sol(HashTable* ht_elem_ptr,HashTable* ht_node_ptr,
		    Neigh_Sol* neigh_list,Neigh_Sol_Pack** rec_el_pack,int *elem_num,int* rec_elem_num)
{
  int i,j,k;
  int flag;
  int no_migrate_elem;
  
  int myid,numprocs;
  Neigh_Sol* tempsol;
  
  Neigh_Sol_Pack* el_pack; //for packing element solution to be sent to diff proc
  
  MPI_Status status;
  MPI_Request   request;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  
  MPI_New_Datatype();
  
  int kk;
  for(i=0;i<numprocs;i++)       // loop to send interface elements
    {     
     if(i!=myid)
       {
	 if(*(elem_num +i) <=0)
	   { 
	     kk=-1;
	     MPI_Isend(&kk,1,MPI_INT,i,7, MPI_COMM_WORLD,&request);
	     MPI_Wait(&request, &status);
	   }
	 else
	   {
	     el_pack = new Neigh_Sol_Pack[*(elem_num + i)];
	     tempsol = neigh_list[i].next;
	     
	     j=0;
	     while(tempsol)
	       {
		 Pack_Sol(tempsol,&el_pack[j]);
		 tempsol =tempsol->next;
		 j++;
	       }
	     kk =*(elem_num+i);
	     MPI_Isend(&kk, 1,MPI_INT,i,7, MPI_COMM_WORLD,&request);
	     MPI_Wait(&request, &status);
	     
	     if(*(elem_num + i)>0)       
	       {
		 MPI_Isend(el_pack,*(elem_num + i),NSOLTYPE,i,8, MPI_COMM_WORLD,&request);
		 MPI_Wait(&request, &status);
	       }
	     
	     delete [] el_pack;
	     
	   }// end of else
	 
       }//end of outer else
    }//end of for
  
  
  //code to receive elements
  
  
  
  for(i=0;i<numprocs;i++)
    *(rec_elem_num +i) = -1;
    
  for(i=0;i<numprocs;i++) // loop to receive interface elements 
    {   
      
      if(i!=myid)
	{
	  
	  MPI_Irecv(&rec_elem_num[i],1,MPI_INT,i,7, MPI_COMM_WORLD, &request);
          MPI_Wait(&request, &status);
	  
          if(rec_elem_num[i]>0)
	    {
	      rec_el_pack[i] = new Neigh_Sol_Pack[rec_elem_num[i]];
	      no_migrate_elem += *(rec_elem_num+i);
	      MPI_Irecv(rec_el_pack[i],rec_elem_num[i],NSOLTYPE,i,8, MPI_COMM_WORLD,&request);
	      MPI_Wait(&request, &status);
	      
	    }
	  
	  
        }//end of  if(i!=myid)
    } // end of  for
  

  
  
  MPI_Barrier(MPI_COMM_WORLD); 
}// end of func








void Pack_Sol(Neigh_Sol* tempsol,Neigh_Sol_Pack* el_Pack)
{
  int i,j,k,ii,jj;
  
  for(i=0;i<KEYLENGTH;i++)
    {
      el_Pack->key[i] = tempsol->key[i];
    }
  
  for(i=0;i<MAX_ORDER;i++)
    {  
      el_Pack->norder[i] =  tempsol->norder[i];
      
    }
  
  el_Pack->nside = tempsol->nside;
  
  for( ii=0;ii<2;ii++)
    for(jj=0;jj<121;jj++)
      {
	el_Pack->solu[ii][jj] = tempsol->solu[ii][jj];
      } 
  
  for(ii=0;ii<18;ii++)
    {
      el_Pack->Xnod[ii] = tempsol->Xnod[ii];
    }    
  
}//end of func




