#include "../header/hpfem.h"

void calc_volume(HashTable* El_Table, HashTable* NodeTable, 
		 int myid, int nump, MatProps* matprops_ptr) 
{
  int i;
  double volume = 0, max_x=0, max_y=0, min_x=1.E30, min_y=1.E30;
  double ndcoord[2];
  Node* NodePtr;

  //-------------------go through all the elements of the subdomain and  
  //-------------------calculate the state variables at time .5*delta_t

  HashEntryPtr* buck = El_Table->getbucketptr();
  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i))
      {
	HashEntryPtr currentPtr = *(buck+i);
	while(currentPtr)
	  {
	    Element* Curr_El=(Element*)(currentPtr->value);
 	    if(Curr_El->get_refined_flag() == 0){
	      volume += Curr_El->calc_volume(NodeTable);
		// Abani introduces to get "flow extent"
	      //	if (*(Curr_El->get_el_solution()) > .005) {

	      //	  NodePtr=(Node*) (NodeTable->lookup(Curr_El->pass_key()));
	      //	  ndcoord[0]=*(NodePtr->get_coord());
	      //  ndcoord[1]=*(NodePtr->get_coord()+1);

	      //  double xx = 
	      //	    (ndcoord[0])*(matprops_ptr->LENGTH_SCALE);
	      ///	  double yy = 
		//    (ndcoord[1])*(matprops_ptr->LENGTH_SCALE);
		//  if (max_x < xx) max_x=xx;
		///	  if (max_y < yy) max_y=yy;
		//  if (min_y > yy) min_y=yy;
		//  if (min_x > xx) min_x=xx;
		//	} // end of extent

	    }
	    currentPtr=currentPtr->next;      	    
	  }
      }
  double gl_volume,gl_x_max=0,gl_y_max=0,gl_x_min=0,gl_y_min=0;
  MPI_Barrier(MPI_COMM_WORLD);
  i = MPI_Reduce(&volume, &gl_volume, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  //Abani introuduce
  // i = MPI_Reduce(&max_x, &gl_x_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  // i = MPI_Reduce(&max_y, &gl_y_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  // i = MPI_Reduce(&min_x, &gl_x_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  // i = MPI_Reduce(&min_y, &gl_y_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

  if(myid == 0) {
    gl_volume = gl_volume * (matprops_ptr->LENGTH_SCALE) * (matprops_ptr->LENGTH_SCALE) * (matprops_ptr->HEIGHT_SCALE);
    printf("volume is %e\n", gl_volume);
  }  

  return;
}

