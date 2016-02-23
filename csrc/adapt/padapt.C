#include "../header/hpfem.h"

void  P_adapt(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, int h_count, 
	      double target, MatProps* matprops_ptr)
  /*-------------
    scanning element hashtable, if the error of an element is bigger than the target error
    1). check if it has been marked for refinement caused by other element refinement
        if not, store it in the 'refined' array. at the same time, counting the related
	refinement. if it is bigger than a limittation. refuse to do the refinement and
	remove the refinement mark of every related element.
    2). at the same time, if one refinement is required in an element along the interfaces, 
        check its neighbor, if this neighbor also need to be refinement, stop this round of
	checking and refuse to do the refinement. how to make the refinement infomation be
	passed through all the processors is very, very, very difficult.
    3). continue untill every element has been checked

    *------------------------------------------------------*/

{

  int k, j, i;
  HashEntryPtr   entryp;
  Element*       EmTemp;

  Element*       refined[297200];//maybe with linked list or creating new arrays while running
  int            count=0;
  int            myid;
  int            numprocs;
  int            ifg;//--
  int            refine_flag;
  unsigned       sent_buf[4*KEYLENGTH];//-- 1 source; 2 target; 3 refined flag; 4 generation
  unsigned       recv_buf[4*KEYLENGTH];//-- same
  int            htype = 101;//-- 101 is arbitary
  
  int hash_size=HT_Elem_Ptr->get_no_of_buckets();
  for(i=0;i<hash_size;i++)
    {
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
	{ 
	  EmTemp = (Element*)(entryp->value);
	  EmTemp->put_new_old(OLD);
	  entryp = entryp->next;
	}
    }
  

  for(k=0;k<297200;k++) refined[k] = 0;

  int h_begin = 1;
  int h_begin_type = 102;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  //clear out ghost elements
  delete_ghost_elms(HT_Elem_Ptr, myid);
  // determine which elements to refine and flag them for refinement
  double geo_target = element_weight(HT_Elem_Ptr, HT_Node_Ptr, myid, numprocs);
  //test_h_refine(HT_Elem_Ptr, myid, h_count);//--debugging

  //MPI_Barrier(MPI_COMM_WORLD); //-- every process arrive here

  int debug_ref_flag = 0;
  for(i=0; i<hash_size; i++)//-- every process begin to scan their own hashtable
    {
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
	{	  
	  EmTemp = (Element*)(entryp->value);
	  //-- this requirement is used to exclude the new elements
	  if(EmTemp->get_new_old()==OLD && !(EmTemp->get_refined_flag()))
	    {
	      if(*(EmTemp->get_el_error()) > geo_target+WEIGHT_ADJUSTER &&
		 EmTemp->get_gen() < REFINE_LEVEL)
		refine_flag = 1;
	      else  
		refine_flag = 0; 

	      // refinement indicator hack!!!
	      double* d_temp = EmTemp->get_el_solution();
	      double pile = 0;
	      for(j=0;j<4;j++)
		pile += d_temp[j]*d_temp[j];
	      if(pile > GEOFLOW_TINY && EmTemp->get_gen() < REFINE_LEVEL)
		refine_flag = 1;
	      else  
		refine_flag = 0; 

	      if(refine_flag != 0)
		{
		  //		  printf("refining right now!!!\n");
		  int sur = 0; int mii = 0;
		  while(refined[mii])
		    {  //-- KEYLENGTH should be considered
		      if(*(refined[mii]->pass_key()) == *(EmTemp->pass_key()) &&
			 *(refined[mii]->pass_key()+1) == *(EmTemp->pass_key()+1))
			{
			  sur = 1;
			  break;
			}
		      mii++;
		    }
		  if(!sur)
		    {
		      ifg = 1; int j = 0;
		      //-- check out the triggered refinement
		      depchk(EmTemp, HT_Elem_Ptr, &ifg, refined, &count);
		      if(ifg)
			{ 
			  //printf("%u %u is getting refined on proc %d\n",*(EmTemp->pass_key()), *(EmTemp->pass_key()+1), myid);
			  while(refined[j])
			    { 
			      if(!refined[j]->get_refined_flag())
				refine(refined[j], HT_Elem_Ptr, HT_Node_Ptr, matprops_ptr); 
			      j++; 
			    }
			}
		    }
		}
	      debug_ref_flag++;
	    }
	  entryp=entryp->next;
	}
    }

  /*-h_count for debugging*/
  update_neighbor_info(HT_Elem_Ptr, refined, count, myid, numprocs, HT_Node_Ptr, h_count);  

  //all_check(HT_Elem_Ptr, HT_Node_Ptr, myid, 3);

  //data_com(HT_Elem_Ptr, HT_Node_Ptr, myid, numprocs, h_count); //replaced by update_interproc.C
  /* transfer ghost elements to proper processors */
  
  move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr);

  return;
}



