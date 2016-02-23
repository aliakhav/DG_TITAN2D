#include "../header/hpfem.h"

void Assemble_mass_matrix(double el_stiff[], double el_load[], Element* Curr_El, 
			  HashTable* El_Table, HashTable* NodeTable, int myid)
{
  int i,j,k, counter;
  double end, start, schur_time = 0;

  assert(!Curr_El->get_refined_flag());
  //print out coordinate of element centroid
  /*	Node* nd = (Node*) NodeTable->lookup(Curr_El->pass_key());
	printf(" i is : %d x is %e y is %e\n",i, *(nd->get_coord()), *(nd->get_coord()+1));*/
                
  int ndof = Curr_El->get_no_of_dof();
  int* ix_el = new int[ndof];  //elm to subdof mapping
  int ix_counter = 0;
  
  for(j=0; j<8; j++)
    {

    }
  /*	Node* NodePtr=(Node*) (NodeTable->lookup(Curr_El->pass_key()));
	for(int k=*(NodePtr->getdof()); k<*((NodePtr->getdof())+1)+1; k++)
	{
	ix_el[ix_counter] = k;
	ix_counter++;
	}*/  //bubble functions get condensed out
  //element to subdomain stiffness matrix mapping done
  
  //----assemble the sub_domain stiffness-------
  
  int count1=0;
  int count2=0;
  Curr_El->get_stiffness(NodeTable,El_Table, el_stiff, 
			 el_load, Curr_El);//--get the element stiffness and rhs	

  return;
}
