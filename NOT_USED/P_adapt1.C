#include  "../header/hpfem.h"
#include  "../header/enriched_neighbor_info.h"


extern void update_order_interp(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, 
				enriched_neighbor* refined_start, int myid, int numprocs);

extern void  data_com(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, 
		      int myid, int numprocs, int h_count);

extern void  htflush(HashTable*, HashTable*, int);


void P_adapt(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
	     int h_count, double TARGET)
{
  Element*  enriched[297200];/*maybe with linked list or creating new arrays while running*/
  int counter=0;
  int i, j, k;
  Element* EmTemp;
  Element* Neighbor;
  Element* EmTemp2;

  Node*  NdTemp;
  void* p;
  unsigned* MyKey;
  unsigned* NeiKey;
  unsigned* neighbors;
  int* order;
  int* neigh_proc;
  int myid, numprocs;
  HashEntryPtr entryp;
  double error;

  enriched_neighbor* enriched_start=new enriched_neighbor();
  enriched_neighbor* enriched_current=enriched_start;
  enriched_neighbor* enriched_new;

  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  htflush(HT_Elem_Ptr, HT_Node_Ptr, 1);/*-- new -> old*/

  int e_buckets=HT_Elem_Ptr->get_no_of_buckets();

  for(i=0;i<e_buckets;i++)/*-- increase the order of element seperately*/
    {
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
	{ 
	  EmTemp = (Element*)(entryp->value);	  
	  
	  if(!EmTemp->get_refined_flag())
	    {
	      error=sqrt(*(EmTemp->get_el_error()));
	      if(error > TARGET)
		{
		  enriched[counter]=EmTemp;
		  counter++;
		  order = EmTemp->get_order();
		  for(j=0;j<5;j++)              //-- increase order of four sides and bubble
		    if(*(order+j)<MAX_ORDER)		    
		      EmTemp->put_order(j, *(order+j)+1);		  
		}

	    }
	  entryp = entryp->next;
	}

    }

  /*cout<<myid<<" No. of Elements P-adapted: "<<counter<<endl<<flush;*/
  printf("%d No. of Elements P-adapted: %d\n",myid,counter); 

  //-- equivalize the order of the common edge of two elems
  for(i=0; i < counter; i++)
    {   

      EmTemp=enriched[i];

      /*################################################
	added by jp on oct12 
	purpose: update the error of enriched element
	###############################################*/
      int ord[5];
      int ord_max1;
      int ii, m;
      for(ii=0;ii<5;ii++)
	ord[ii] = *(EmTemp->get_order()+ii);
      max_order(ord, &ord_max1);
      
      
      order = EmTemp->get_order();
      neigh_proc = EmTemp->get_neigh_proc();
      neighbors  = EmTemp->get_neighbors();
      MyKey = EmTemp->pass_key();
      
      for(j=0;j<4;j++)
	{
	  int my_order = *(order+j);
	  if(*(neigh_proc+j) == -1)
	    {
	      p = HT_Node_Ptr->lookup(EmTemp->getNode()+(j+4)*KEYLENGTH);
	      NdTemp = (Node*)p;
	      NdTemp->put_order(*(order+j));		      
	    }
	  
	  if(*(neigh_proc+j) == myid)
	    {
	      
	      p = HT_Elem_Ptr->lookup(neighbors+j*KEYLENGTH);
	      Neighbor = (Element*)p;
	      int which_side = Neighbor->which_neighbor(MyKey);
	      if(which_side > 3) which_side +=-4;
	      int neigh_order = *(Neighbor->get_order()+which_side);	      
	      
	      int final_order = (my_order>=neigh_order) ? my_order:neigh_order;
	      
	      if(EmTemp->get_gen() <= Neighbor->get_gen())
		{		  	       
		  
		  if(EmTemp->get_gen() < Neighbor->get_gen())//smaller neighbor
		    {
		      assert(*(neigh_proc+j+4)==myid);
		      p = HT_Elem_Ptr->lookup(neighbors+(j+4)*KEYLENGTH);
		      Neighbor = (Element*)p;	
		      //which side should remain the same
		      //which_side = Neighbor->which_neighbor(MyKey);
		      assert(Neighbor->which_neighbor(MyKey)==which_side);

		      neigh_order = *(Neighbor->get_order()+which_side);
		      
		      if(neigh_order >= final_order)
			final_order=neigh_order;				     		
		      
		      else
			{
			  Neighbor->put_order(which_side, final_order);
			  NdTemp=(Node*)HT_Node_Ptr->lookup(Neighbor->getNode()+(which_side+4)*KEYLENGTH);
			  NdTemp->put_order(final_order);
			}
		      
		    }

		  Neighbor=(Element*)HT_Elem_Ptr->lookup(neighbors+j*KEYLENGTH);
		  Neighbor->put_order(which_side, final_order);
		  NdTemp=(Node*)HT_Node_Ptr->lookup(Neighbor->getNode()+(which_side+4)*KEYLENGTH);
		  NdTemp->put_order(final_order);

		  EmTemp->put_order(j, final_order);
		  NdTemp=(Node*)HT_Node_Ptr->lookup(EmTemp->getNode()+(j+4)*KEYLENGTH);
		  NdTemp->put_order(final_order);
		  
		}
	      else if(EmTemp->get_gen() > Neighbor->get_gen())
		//lookup the other small elements which has the same neighbor
		//get the max of the three
		{
		  int third_order=0;		  
		  int start; int end; int twin;

		  if(j<3) { start = j; end = j+1;}
		  else { start = 3; end = 0;}
		  
		  p = HT_Node_Ptr->lookup(EmTemp->getNode()+start*KEYLENGTH);
		  NdTemp = (Node*)p;
		  if(NdTemp->getinfo() == S_C_CON) twin = start - 1;
		  else twin = end;
		  if(twin<0) twin = 3;
		  
		  p = HT_Elem_Ptr->lookup(neighbors+twin*KEYLENGTH);
		  EmTemp2 = (Element*)p;
		  assert(EmTemp2);
		  assert(*(EmTemp2->get_neighbors()+j*KEYLENGTH)==*(Neighbor->pass_key()) && *(EmTemp2->get_neighbors()+j*KEYLENGTH+1)==*(Neighbor->pass_key()+1));
		  third_order=*(EmTemp2->get_order()+j);
		      
		
		  if(my_order>=neigh_order)final_order=my_order;
		  else final_order=neigh_order;
		  if(third_order>final_order)
		    final_order=third_order;
		  
		  Neighbor->put_order(which_side%4, final_order);//which_side can be > 3

		  EmTemp->put_order(j, final_order);
		  EmTemp2->put_order(j, final_order);
		  
		  NdTemp=(Node*)HT_Node_Ptr->lookup(Neighbor->getNode()+(which_side%4+4)*KEYLENGTH);
		  NdTemp->put_order(final_order);
		  
		  NdTemp=(Node*)HT_Node_Ptr->lookup(EmTemp->getNode()+(j+4)*KEYLENGTH);
		  NdTemp->put_order(final_order);

		  NdTemp=(Node*)HT_Node_Ptr->lookup(EmTemp2->getNode()+(j+4)*KEYLENGTH);
		  NdTemp->put_order(final_order);	  
		}
	    }
	  
	  else if((*(neigh_proc+j)>=0)&&(*(neigh_proc+j)!=myid))
	    {
	      if(*(EmTemp->get_neigh_gen()+j) < EmTemp->get_gen())
		{
		  
		  int start; int end; int twin;
		  if(j<3) { start = j; end = j+1;}
		  else { start = 3; end = 0;}
		  
		  p = HT_Node_Ptr->lookup(EmTemp->getNode()+start*KEYLENGTH);
		  NdTemp = (Node*)p;
		  if(NdTemp->getinfo() == S_C_CON) twin = start - 1;
		  else twin = end;
		  if(twin<0) twin = 3;
		  
		  p = HT_Elem_Ptr->lookup(neighbors+twin*KEYLENGTH);
		  Neighbor = (Element*)p;
		  //assert(Neighbor);
		  //if(!Neighbor) cout<<myid<<"  "<<*(neighbors+twin*KEYLENGTH)<<" Neighbor"<<"\n"<<flush;
		  
		  NeiKey = EmTemp->get_neighbors()+j*KEYLENGTH;
		  
		  int flag;
		  for(m=0;m<4;m++)
		    {
		      flag = 1;
		      for(k=0;k<KEYLENGTH;k++)
			{
			  if( *(Neighbor->get_neighbors()+m*KEYLENGTH+k)
			      != *(NeiKey+k))
			    {
			      flag = 0;
			      break;
			    }
			}
		      if(flag) break;
		    }

		  assert(m==j);
		  //assert(flag);
		  //if(!flag) cout<<myid<<"  "<<*(NeiKey+k)<<" flag"<<"\n"<<flush;
		  
		  int neigh_order = *(Neighbor->get_order()+m);
		  int my_order = *(order+j);
		  
		  int final_order = (my_order>=neigh_order) ? my_order:neigh_order;
		  
		  EmTemp->put_order(j, final_order);
		  Neighbor->put_order(m, final_order);
		  
		  NdTemp->put_order(final_order);//still the S_C_CON

		  NdTemp=(Node*)HT_Node_Ptr->lookup(EmTemp->getNode()+(j+4)*KEYLENGTH);
		  NdTemp->put_order(final_order);

		  NdTemp=(Node*)HT_Node_Ptr->lookup(Neighbor->getNode()+(m+4)*KEYLENGTH);
		  NdTemp->put_order(final_order);	  
		  
		  /*inserting the element in the linked list finally*/
		  /*there might be elements which are not necessary, if two small
		   elements at the interface of a big 1 were enriched but don't care*/
		  enriched_new=new enriched_neighbor();
		  enriched_current->next=enriched_new;

		  enriched_new->set_parameters(*(neigh_proc+j), NeiKey, MyKey, final_order);

		  enriched_current=enriched_new;

		}//if diff.proc and EmTemp is younger
	      else//if EmTemp and its neighbor a.) has same gen. b.)EmTemp is older
		{
		  
		  NeiKey = EmTemp->get_neighbors()+j*KEYLENGTH;
		  NdTemp=(Node*)HT_Node_Ptr->lookup(EmTemp->getNode()+(j+4)*KEYLENGTH);
		  NdTemp->put_order(*(order+j));
		  enriched_new=new enriched_neighbor();
		  enriched_current->next=enriched_new;
		  enriched_new->set_parameters(*(neigh_proc+j), NeiKey, MyKey, *(order+j));
		  enriched_current=enriched_new;
		  
		  
		}

	    }//if different proc

	}//for j=0 j<4
      
      /*################################################
	added by jp on oct12 
	purpose: update the error of enriched element
	###############################################*/
      for(ii=0;ii<5;ii++)
	ord[ii] = *(EmTemp->get_order()+ii);
      int ord_max2;
      max_order(ord, &ord_max2);
      
      int gn;
      double sl, er;
      gn = EmTemp->get_gen();
      sl = *EmTemp->get_el_solution();
      er = *EmTemp->get_el_error();      
      double hh = 1.0/pow(2,(gn+1));
      
      EmTemp->putel_sq(hh*pow(((double)ord_max1/(double)ord_max2)*sl,2),
		       hh*pow(((double)ord_max1/(double)ord_max2)*er,2));
      /*----------------------------------------------------------*/

    }//going through all the enriched elements
  

  MPI_Barrier(MPI_COMM_WORLD);

  update_order_interp(HT_Elem_Ptr, HT_Node_Ptr, enriched_start, myid, numprocs);
  htflush(HT_Elem_Ptr, HT_Node_Ptr, 2);
  //data_com(HT_Elem_Ptr, HT_Node_Ptr, myid, numprocs, h_count);
  

  /*this recalculation should be included in put_order of Element class*/
  int hash_size=HT_Elem_Ptr->get_no_of_buckets();
  for(i=0;i<hash_size;i++)
    {
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
	{ 
	  EmTemp = (Element*)(entryp->value);	
	  
	  if(!EmTemp->get_refined_flag())
	    EmTemp->update_ndof();

	  entryp = entryp->next;
	}
    }
}
