#include "../header/hpfem.h"
#include "../header/geoflow.h"

double get_coef_and_eigen(HashTable* El_Table, HashTable* NodeTable,
			  MatProps* matprops_ptr, int ghost_flag)
{
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  double min_distance = 1000000, max_evalue = GEOFLOW_TINY, junk;

  int i,j,k, counter;
  double tiny = GEOFLOW_TINY, min_dx_dy_evalue = 10000000, hmax = 0;
  int el_counter = 0;
  double evalue = 1.0;  //might need to change this
  //-------------------go through all the elements of the subdomain and get 
  //-------------------the coefficients and eigenvalues and calculate the time step

  HashEntryPtr* buck = El_Table->getbucketptr();
  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i))
      {
	HashEntryPtr currentPtr = *(buck+i);
	while(currentPtr)
	  {
	    Element* Curr_El=(Element*)(currentPtr->value);
	    int refined_flag = Curr_El->get_refined_flag();
 	    if(refined_flag == 0 || (refined_flag == GHOST && ghost_flag == 1))//if this is a refined element don't involve!!!
	      {
		/* calculate hmax */
		if(hmax < *(Curr_El->get_state_vars()))
		  hmax = *(Curr_El->get_state_vars());
		double* d_uvec = Curr_El->get_d_state_vars();
		double * dx_ptr = Curr_El->get_dx();
#ifdef SUNOS
		gmfggetcoef_(Curr_El->get_state_vars(), d_uvec, (d_uvec+NUM_STATE_VARS), dx_ptr, 
				 &(matprops_ptr->bedfrict), &(matprops_ptr->intfrict), Curr_El->get_kactxy(),
				 (Curr_El->get_kactxy()+1), &tiny, &(matprops_ptr->epsilon));
		eigen_(Curr_El->get_state_vars(), (Curr_El->get_eigenvxymax()),
		       (Curr_El->get_eigenvxymax()+1), &evalue, &tiny, 
		       Curr_El->get_kactxy(), Curr_El->get_gravity());
#endif
#ifdef IBMSP
#endif
		//************************************************************
		//!!!!!!!!!!!!!!!!!!!!!!check dx & dy!!!!!!!!!!!!!!!!!!!!!!!!
		//************************************************************
		junk = c_dmin1(dx_ptr[0],dx_ptr[1]);
		if(junk/evalue < min_dx_dy_evalue) {
		  min_distance = junk;
		  max_evalue = evalue;
		}
		if(evalue > 1000000000.) {
		  printf(" eigenvalue is %e for procd %d momentums are %e %e for pile height %e\n",
			 evalue, myid, *(Curr_El->get_state_vars()+1), *(Curr_El->get_state_vars()+2),
			 *(Curr_El->get_state_vars()));		  
		}

		  min_dx_dy_evalue = c_dmin1(c_dmin1(dx_ptr[0],dx_ptr[1])/evalue, min_dx_dy_evalue);
	      }
	    currentPtr=currentPtr->next;      	    
	  }
      }
  double global_dt[2], dt[2];
  dt[0] = 0.5*min_dx_dy_evalue;
  dt[1] = -0.9*hmax;
  
  i = MPI_Allreduce(dt, global_dt, 2, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  dt[0] = 0.5*c_dmin1( global_dt[0], -global_dt[1]);
  if(myid == 0) {
    if(global_dt[0] < -global_dt[1])
      printf(" the time step is being set by the eigenvalues\n");
    else
      printf(" the time step is being set by the pile height\n");
  }
  printf("====== for proc %d time step should be %e from distance %e and evalue %e  =========\n",
	 myid, min_distance/max_evalue, min_distance, max_evalue);
    
  
  return dt[0];
}

