#include "../header/hpfem.h"
#include <math.h>
#include "../header/blas.h"

//#define PRINT_GIS_ERRORS

#ifdef SUNOS
extern "C" void elemcom_(int*, int*, int*, int*, int*, int*, int*, 
			 double[][2], double[], double[], int*, double*);
extern "C" void masselem_(int*, int*, int*, int*, int*, double*,
			  double*, double*, double*, double*, double*, double*, int*);
#endif

#ifdef IBMSP
extern "C" void elemcom(int*, int*, int*, int*, int*, int*, int*, 
			double[][2], double[], double[], int*, double*);
extern "C" void masselem(int*, int*, int*, int*, int*, double*, double*,
			 double*, double*, double*, double*, double*, double*, int*);
#endif


/*  original element   */

Element::Element(unsigned nodekeys[][2], unsigned neigh[][2], int n_pro[], 
		 int mat, double* pile_height, int myid, int elm_loc_in[],
		 unsigned opposite_brother[])
{ 
  generation = 0;//--first generation 
  fluxIntegral= 0.;
  material=mat;
  int i, j;
  for(i=0; i<KEYLENGTH; i++)
    key[i]=nodekeys[8][i];//--using bubble key to represent the element
  
  for(i=0; i<8; i++)
    for(int j=0; j<KEYLENGTH; j++)
      node_key[i][j]=nodekeys[i][j];
  
  for(i=0; i<4; i++)
    {
      neigh_proc[i]=n_pro[i];
      neigh_proc[i+4] = -2;//-- -2 means regular element
      if(neigh_proc[i]!=-1)
	for(int j=0; j<KEYLENGTH; j++)
	  neighbor[i][j]=neighbor[i+4][j]=neigh[i][j];
    }
    
  order=ORDER; 

  for(i=0;i<EQUATIONS;i++)
    el_error[i]=0;
  
  for(i=0;i<8;i++)
    neigh_gen[i] = 0;
  
  ndof=calc_num_shape_functions(order)*EQUATIONS;
  
  refined=0;
  
  for(i=0;i<4;i++)
    for(j=0;j<KEYLENGTH;j++)
      brothers[i][j] = 0;

  myprocess = myid;
  
  new_old = OLD;
  elm_loc[0] = elm_loc_in[0];
  elm_loc[1] = elm_loc_in[1];
  calc_which_son();
  for(i=0;i<KEYLENGTH;i++) {
    brothers[which_son][i] = key[i];
    brothers[(which_son+2)%4][i] = opposite_brother[i];
  }
  switch(which_son) {
  case 0:
    for(i=0;i<KEYLENGTH;i++) {
      brothers[1][i] = neighbor[1][i];
      brothers[3][i] = neighbor[2][i];
    }
    break;
  case 1:
    for(i=0;i<KEYLENGTH;i++) {
      brothers[0][i] = neighbor[3][i];
      brothers[2][i] = neighbor[2][i];
    }
    break;
  case 2:
    for(i=0;i<KEYLENGTH;i++) {
      brothers[1][i] = neighbor[0][i];
      brothers[3][i] = neighbor[3][i];
    }
    break;
  case 3:
    for(i=0;i<KEYLENGTH;i++) {
      brothers[0][i] = neighbor[0][i];
      brothers[2][i] = neighbor[1][i];
    }
    break;
  }  
  opposite_brother_flag = 1;

  // initialize the solution values
  for(i=0;i<ELM_DOF;i++) {
    el_solution[i] = 0;
    prev_el_solution[i] = 0;
    temp_el_solution[i] = 0;

  }

  for(i=0;i<3;i++) {
    el_solution[i*EQUATIONS] = pile_height[i];
    prev_el_solution[i*EQUATIONS] = pile_height[i];
  }
  
  // if (el_solution[0] != 0) {
  //   el_solution[1] = 0.0001;
  //   prev_el_solution[1]=0.0001;
  //}   

  //printf("creating an original element\n");

}

//used for refinement
Element::Element(unsigned nodekeys[][2], unsigned neigh[][2], 
		 int n_pro[], int gen,  
		 int ord, int gen_neigh[], int mat, 
		 Element* fthTemp, HashTable* El_Table, 
		 HashTable* NodeTable, int myid, MatProps* matprops_ptr,
		 int elm_loc_in[])
{ 
  generation = gen;//--first generation
  fluxIntegral=0.;
  material = mat;
  int i;
  for(i=0; i<KEYLENGTH; i++)
    key[i]=nodekeys[8][i];//--using buble key to represent the element

  elm_loc[0] = elm_loc_in[0];
  elm_loc[1] = elm_loc_in[1];
  opposite_brother_flag = 1;
  for(i=0; i<8; i++)
    for(int j=0; j<KEYLENGTH; j++)
      node_key[i][j]=nodekeys[i][j];
  
  for(i=0; i<4; i++)
    {
      neigh_proc[i]=n_pro[i];
      neigh_proc[i+4] = -2;//-- -2 means regular element
      if(neigh_proc[i]!=-1)
	for(int j=0; j < KEYLENGTH; j++)
	  neighbor[i][j]=neighbor[i+4][j]=neigh[i][j];
      
      neigh_gen[i]=neigh_gen[i+4]=gen_neigh[i];
    }
  
  order=ord;
 
  for(i=0;i<EQUATIONS;i++)
    el_error[i]=0;

  ndof=calc_num_shape_functions(order)*EQUATIONS;

  lb_weight=0;
  
  refined=0;
  
  new_old = NEW;
  for(i=0;i<ELM_DOF;i++) {
    el_solution[i] = 0;
    prev_el_solution[i] = 0;
    temp_el_solution[i] = 0;
  }
  
  return;
}
/*********************************
 making a father element from its sons
*****************************************/
Element::Element(Element* sons[], HashTable* NodeTable,
		 HashTable* El_Table, MatProps* matprops_ptr) {
  int i, j;
  unsigned* son_nodes[4];
  opposite_brother_flag = 0;
  for(i=0;i<4;i++)
    son_nodes[i] = sons[i]->getNode();

  for(i=0;i<KEYLENGTH;i++) {
    node_key[0][i] = son_nodes[0][i];
    node_key[1][i] = son_nodes[1][KEYLENGTH+i];
    node_key[2][i] = son_nodes[2][2*KEYLENGTH+i];
    node_key[3][i] = son_nodes[3][3*KEYLENGTH+i];
    node_key[4][i] = son_nodes[0][KEYLENGTH+i];
    node_key[5][i] = son_nodes[1][2*KEYLENGTH+i];
    node_key[6][i] = son_nodes[2][3*KEYLENGTH+i];
    node_key[7][i] = son_nodes[3][i];
    key[i] = son_nodes[0][2*KEYLENGTH+i];
    elm_loc[i] = (*(sons[0]->get_elm_loc()+i))/2;
  }
  myprocess = sons[0]->get_myprocess();
  generation = sons[0]->get_gen() - 1;
  material = sons[0]->get_material();
  fluxIntegral=0.;
  for(i=0;i<EQUATIONS;i++)
    el_error[i]=0;
  refined=1; // not an active element yet!!!
  
  calc_which_son();
  //order information -- keep the highest order
  order = 0;
  for(i=0;i<4;i++)
    if(order < sons[i]->get_order())
      order = sons[i]->get_order();

  ndof = calc_num_shape_functions(order)*EQUATIONS;
  lb_weight=ndof;

  // neighbor information
  for(i=0;i<4;i++) {
    neigh_gen[i+4] = neigh_gen[i] = *(sons[i]->get_neigh_gen()+i);
    neigh_proc[i] = *(sons[i]->get_neigh_proc()+i);
    if(neigh_gen[i] > generation && neigh_proc[i] != -1) {
      neigh_proc[i+4] = neigh_proc[i];
      if(neigh_gen[i] <= generation-1)
	assert(neigh_gen[i] <= generation-1);
    }
    else
      neigh_proc[i+4] = -2;
  }
  for(i=0;i<4;i++) 
    for(j=0;j<KEYLENGTH;j++) 
      neighbor[i][j] = *(sons[i]->get_neighbors()+i*KEYLENGTH+j);

  for(j=0;j<KEYLENGTH;j++)
    neighbor[7][j] = *(sons[0]->get_neighbors()+7*KEYLENGTH+j);

  for(i=0;i<3;i++) 
    for(j=0;j<KEYLENGTH;j++) 
      neighbor[i+4][j] = *(sons[i+1]->get_neighbors()+i*KEYLENGTH+j);
   
  /* brother information -- requires that atleast one of this
     element's neighboring brothers is on this process in 
     order to get information on the brother that is not a neighbor */
  Element* EmTemp;
  switch(which_son) {
  case 0:
    for(i=0;i<KEYLENGTH;i++)
      brothers[0][i] = key[i];
    if(neigh_proc[1] == -1) {
      for(i=0;i<KEYLENGTH;i++)
	brothers[1][i] = 0;
    }
    else if(neigh_gen[1] == generation) {
      for(i=0;i<KEYLENGTH;i++)
	brothers[1][i] = neighbor[1][i];
    }
    else if(neigh_gen[1] == generation+1) {
      EmTemp = (Element*) El_Table->lookup(neighbor[1]);
      assert(EmTemp);
      unsigned* bro_key = EmTemp->getfather();
      for(i=0;i<KEYLENGTH;i++)
	brothers[1][i] = bro_key[i];
    }
    else
      assert(0);
    if(neigh_proc[2] == -1) {
      for(i=0;i<KEYLENGTH;i++)
	brothers[3][i] = 0;
    }
    else if(neigh_gen[2] == generation) {
      for(i=0;i<KEYLENGTH;i++)
	brothers[3][i] = neighbor[2][i];
    }
    else if(neigh_gen[2] == generation+1) {
      EmTemp = (Element*) El_Table->lookup(neighbor[2]);
      assert(EmTemp);
      unsigned* bro_key = EmTemp->getfather();
      for(i=0;i<KEYLENGTH;i++)
	brothers[3][i] = bro_key[i];
    }
    else
      assert(0);
    break;
  case 1:
    for(i=0;i<KEYLENGTH;i++)
      brothers[1][i] = key[i];
    if(neigh_proc[3] == -1) {
      for(i=0;i<KEYLENGTH;i++)
	brothers[0][i] = 0;
    }
    else if(neigh_gen[3] == generation) {
      for(i=0;i<KEYLENGTH;i++)
	brothers[0][i] = neighbor[3][i];
    }
    else if(neigh_gen[3] == generation+1) {
      EmTemp = (Element*) El_Table->lookup(neighbor[3]);
      assert(EmTemp);
      unsigned* bro_key = EmTemp->getfather();
      for(i=0;i<KEYLENGTH;i++)
	brothers[0][i] = bro_key[i];
    }
    else
      assert(0);
    if(neigh_proc[2] == -1) {
      for(i=0;i<KEYLENGTH;i++)
	brothers[2][i] = 0;
    }
    else if(neigh_gen[2] == generation) {
      for(i=0;i<KEYLENGTH;i++)
	brothers[2][i] = neighbor[2][i];
    }
    else if(neigh_gen[2] == generation+1) {
      EmTemp = (Element*) El_Table->lookup(neighbor[2]);
      assert(EmTemp);
      unsigned* bro_key = EmTemp->getfather();
      for(i=0;i<KEYLENGTH;i++)
	brothers[2][i] = bro_key[i];
    }
    else
      assert(0);
    break;
  case 2:
    for(i=0;i<KEYLENGTH;i++)
      brothers[2][i] = key[i];
    if(neigh_proc[0] == -1) {
      for(i=0;i<KEYLENGTH;i++)
	brothers[1][i] = 0;
    }
    else if(neigh_gen[0] == generation) {
      for(i=0;i<KEYLENGTH;i++)
	brothers[1][i] = neighbor[0][i];
    }
    else if(neigh_gen[0] == generation+1) {
      EmTemp = (Element*) El_Table->lookup(neighbor[0]);
      assert(EmTemp);
      unsigned* bro_key = EmTemp->getfather();
      for(i=0;i<KEYLENGTH;i++)
	brothers[1][i] = bro_key[i];
    }
    else
      assert(0);
    if(neigh_proc[3] == -1) {
      for(i=0;i<KEYLENGTH;i++) 
	brothers[3][i] = 0;
    }
    else if(neigh_gen[3] == generation) {
      for(i=0;i<KEYLENGTH;i++)
	brothers[3][i] = neighbor[3][i];
    }
    else if(neigh_gen[3] == generation+1) {
      EmTemp = (Element*) El_Table->lookup(neighbor[3]);
      assert(EmTemp);
      unsigned* bro_key = EmTemp->getfather();
      for(i=0;i<KEYLENGTH;i++)
	brothers[3][i] = bro_key[i];
    }
    else
      assert(0);
    break;
  case 3:
    for(i=0;i<KEYLENGTH;i++)
      brothers[3][i] = key[i];
    if(neigh_proc[0] == -1) {
      for(i=0;i<KEYLENGTH;i++)
	brothers[0][i] = 0;
    }
    else if(neigh_gen[0] == generation) {
      for(i=0;i<KEYLENGTH;i++)
	brothers[0][i] = neighbor[0][i];
    }
    else if(neigh_gen[0] == generation+1) {
      EmTemp = (Element*) El_Table->lookup(neighbor[0]);
      assert(EmTemp);
      unsigned* bro_key = EmTemp->getfather();
      for(i=0;i<KEYLENGTH;i++)
	brothers[0][i] = bro_key[i];
    }
    else
      assert(0);
    if(neigh_proc[1] == -1) {
      for(i=0;i<KEYLENGTH;i++)
	brothers[2][i] = 0;
    }
    else if(neigh_gen[1] == generation) {
      for(i=0;i<KEYLENGTH;i++)
	brothers[2][i] = neighbor[1][i];
    }
    else if(neigh_gen[1] == generation+1) {
      EmTemp = (Element*) El_Table->lookup(neighbor[1]);
      assert(EmTemp);
      unsigned* bro_key = EmTemp->getfather();
      for(i=0;i<KEYLENGTH;i++)
	brothers[2][i] = bro_key[i];
    }
    else 
      assert(0);
    break;
  }
  find_opposite_brother(El_Table);

  return;
}


void Element::get_stiffness(HashTable* NodeTable, HashTable* HT_Elem_Ptr,
			    double* el_stiffness,
			    double* el_rhs, Element* elmPtr) 

  //for ONE step H-refinement (icon)

{
  int i;
  double ndcoord[9][2];
  int ifg=2;
  int Nc=ndof;
  int Nelb[4];
  int icon[4];
  double bc_value[8];//--for elasticity
  Node* NodePtr;
  Element* ElemPtr;

  for(i=0; i<4; i++) 
    {
      Nelb[i] = 0; 
      icon[i] = 0;
    }
  for(i=0;i<8;i++)
    bc_value[i] = 0;

  if(generation)//filling up icon array
    {
      ElemPtr=(Element*) HT_Elem_Ptr->lookup(getfather());//if the father is not there it is not a constrained element
      
      if(ElemPtr)
	{
	  int j=0; //<---indicates which son it is
	  while(ElemPtr->son[j][0]!=key[0] || ElemPtr->son[j][1]!=key[1])//-- should use KEYLENGTH
	    {
	      j++;
	      if(j==4) 
		{
		  printf("error in get_el_stiffness\n\n");
		  //cerr<<"error in get_el_stiffness\n\n"<<flush;
		  exit(0);
		}
	    }
	  
	  int a=j-1;
	  if(a==-1)a=3;
	  NodePtr=(Node*)(NodeTable->lookup(node_key[a]));
	  assert(NodePtr);
	  
	  if(NodePtr->getinfo()==S_C_CON)
	    icon[a]=-1;
	  
	  a=j+1;
	  if(a==4) a=0;
	  NodePtr=(Node*)(NodeTable->lookup(node_key[a]));
	  assert(NodePtr);
	  
	  if(NodePtr->getinfo()==S_C_CON)
	    icon[j]=1;//-- j was changed to a
	  
	}
    }

  for(i=0; i<8; i++)
    {
      NodePtr=(Node*)(NodeTable->lookup(node_key[i]));
      assert(NodePtr);
      ndcoord[i][0]=NodePtr->coord[0];      
      ndcoord[i][1]=NodePtr->coord[1];
    }
  
  
  
  NodePtr=(Node*) (NodeTable->lookup(key));
  ndcoord[8][0]=NodePtr->coord[0];
  ndcoord[8][1]=NodePtr->coord[1];

/*  
#ifdef SUNOS
  elemcom_(&ifg, &no_of_eqns, &ndof, &Nc, order, Nelb, icon, ndcoord,
	   el_stiffness, el_rhs, &material, bc_value);
#endif

#ifdef IBMSP
  elemcom(&ifg, &no_of_eqns, &ndof, &Nc, order, Nelb, icon, ndcoord,
	   el_stiffness, el_rhs, &material, bc_value);
#endif

*/

}

int Element::which_neighbor(unsigned* FindNeigh)
{
  int i=0, j;
  int found = 0;
  while(found == 0)
    {
      found = compare_key(neighbor[i], FindNeigh);
      if(found != 0 && neigh_proc[i] <0)
        found = 0;
      if(found == 0)
        i++;
      assert(i<8);
    }
    
  return i;
}

void Element::change_neighbor(unsigned* newneighbs, int which_side, int proc, int reg)
{
  int j;
  switch(reg)
    {
    case 1: 
      j = 0;
    case 3: 
      assert(which_side<4);
      for(j=0; j<KEYLENGTH; j++)
	{
	  neighbor[which_side][j]=*(newneighbs+j);
	  neighbor[which_side+4][j]=*(newneighbs+KEYLENGTH+j);           
	}
      neigh_proc[which_side+4]=proc; //assuming no element movement
      neigh_gen[which_side]=neigh_gen[which_side+4]=neigh_gen[which_side]+1;
      break;
    case 4:
      j = 0;
    case 2:
      j = 0;
    case 5:
      for(j=0; j<KEYLENGTH; j++)	
	neighbor[which_side][j]=*(newneighbs+j);   
      neigh_gen[which_side]=neigh_gen[which_side]+1;
      break;
      
    case 6:
      for(j=0; j<KEYLENGTH; j++)	
	neighbor[which_side][j]=neighbor[which_side+4][j]=*(newneighbs+j);     
	neigh_gen[which_side]=neigh_gen[which_side+4]=neigh_gen[which_side]+1;
      break;

      /*Andrew's section called from update_interproc*/
    case 10://the refined element and old neighbor have the same gen.     
      assert(which_side<4);
      for(j=0; j<KEYLENGTH; j++)
	{
	  neighbor[which_side][j]=*(newneighbs+j);
	  neighbor[which_side+4][j]=*(newneighbs+KEYLENGTH+j);           
	}
      neigh_proc[which_side+4]=proc;
      
      neigh_gen[which_side]=neigh_gen[which_side+4]=neigh_gen[which_side]+1;
      break;

    case 11:
      for(j=0; j<KEYLENGTH; j++)
	neighbor[which_side][j]=neighbor[which_side+4][j]=*(newneighbs+j);
      neigh_gen[which_side]=neigh_gen[which_side+4]=neigh_gen[which_side]+1;
      break;


    }
}

void Element::update_ndof()
{
  if(order == 0)
    ndof = EQUATIONS;
  else if(order == 1)
    ndof = 4* EQUATIONS;
  else if(order == 2)
    ndof = 9* EQUATIONS;
  else
    assert(0);
}

void Element::get_nelb_icon(HashTable* NodeTable, HashTable* HT_Elem_Ptr,int* Nelb,int* icon) 

  //for ONE step H-refinement (icon)

{
  int i;
  int ifg=2;
  int Nc=ndof;
  double bc_value[4];//--for poisson equ
  Node* NodePtr;
  Element* ElemPtr;


  for(i=0; i<4; i++) 
    {
      Nelb[i] = 0; 
      icon[i] = 0;
      //bc_value[i] = .0;
    }//-- -1 may be better

  if(generation)//filling up icon array
    {
      ElemPtr=(Element*) HT_Elem_Ptr->lookup(getfather());//if the father is not there it is not a constrained element
      
      if(ElemPtr)
	{
	  int j=0; //<---indicates which son it is
	  while(ElemPtr->son[j][0]!=key[0] || ElemPtr->son[j][1]!=key[1])//-- should use KEYLENGTH
	    {
	      j++;
	      if(j==4) {cerr<<"error in get_el_stiffness\n\n"<<flush; exit(0);}
	    }
	  
	  int a=j-1;
	  if(a==-1)a=3;
	  NodePtr=(Node*)(NodeTable->lookup(node_key[a]));
	  assert(NodePtr);
	  
	  if(NodePtr->getinfo()==S_C_CON)
	    icon[a]=-1;
	  
	  a=j+1;
	  if(a==4) a=0;
	  NodePtr=(Node*)(NodeTable->lookup(node_key[a]));
	  assert(NodePtr);
	  
	  if(NodePtr->getinfo()==S_C_CON)
	    icon[j]=1;//-- j was changed to a
	  
	}
    }

}

Element:: ~Element()
{
/*  if(key[0] == (unsigned) 2501998986) {
    int mmmyid;
    MPI_Comm_rank(MPI_COMM_WORLD, &mmmyid);
    printf("?????????????????????????????????????????????????????? \n");
    printf("?????????????????????????????????????????????????????? \n");
    printf("deleting element %u %u on %d $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n",key[0], key[1], mmmyid);
    printf("?????????????????????????????????????????????????????? \n");
    printf("?????????????????????????????????????????????????????? \n");
    }*/
}


inline double max(double x, double y)
 { 
   if (x>=y) {return (x);}
   else {return(y);}
 }

inline double min(double x, double y)
 {
   if (x>=y) {return (y);}
   else {return(x);}
 }


/*
void Element::calc_gravity_vector(MatProps* matprops_ptr) 
{
  double max_slop = sqrt(zeta[0]*zeta[0]+zeta[1]*zeta[1]);
  double max_angle = atan(max_slope);
    
  double down_slope_gravity = 9.8*sin(max_angle);
  if(dabs(down_slope_gravity) > GEOFLOW_TINY) {
    gravity[0] = -down_slope_gravity*zeta[0]/max_slope;
    gravity[1] = -down_slope_gravity*zeta[1]/max_slope;
    gravity[2] = 9.8*cos(max_angle);
  }
  else {
    gravity[0] = 0;
    gravity[1] = 0;
    gravity[2] = 9.8;
  }

  for(int i=0;i<3;i++)
    gravity[i] = gravity[i]/matprops_ptr->GRAVITY_SCALE;
  
  return;
}
*/
int Element::determine_refinement(double target)
{
  int flag = 0, i;

/*  for(i=0;i<NUM_STATE_VARS*DIMENSION;i++)
    if(dabs(d_state_vars[i]) > target)
      flag = 1;
*/

  return flag;
}

// load-balancing stuff ...
void Element::put_lb_key(unsigned* in_key)
{
  int i;
  for(i=0;i<KEYLENGTH;i++)
    lb_key[i] = in_key[i];
  return;
}

void Element::copy_key_to_lb_key() {
  int i;
  for(i=0;i<KEYLENGTH;i++)
    lb_key[i] = key[i];
  return;
}
// only sends the elements and the necessary parents
int Element::BSFC_check_send_elements(HashTable* HT_Elem_Ptr) 
{
  int i = 0;
  while(i<4) {
    if(neigh_gen[i] < generation && neigh_proc[i] >= 0) {
      Element* EmTemp = (Element*) HT_Elem_Ptr->lookup(getfather());
      if(EmTemp->get_new_old() < 0) {  // first time we looked at the father
	EmTemp->put_new_old(myprocess);
	return 2;
      }
      else if(EmTemp->get_new_old() == myprocess ||
	      EmTemp->get_myprocess() == myprocess) // already getting sent to the correct proc
	return 1;
      else if(EmTemp->get_new_old() >= 0 &&
	      EmTemp->get_myprocess() == -1) { // father sent somewhere else but not to myprocess
	EmTemp->put_myprocess(myprocess);
	return 2;
      }
      else 
	assert(0);
    }
    i++;
  }

  return 1;
}
//send all of the ancestors
int Element::BSFC_check_send_elements() 
{
  return (generation+1);
}

/*
  project_sol takes projects the solution values from a father/son elements onto the son/father
*/
void Element::project_sol(HashTable* NodeTable, HashTable* HT_Elem_Ptr, int parent_son_flag) 
  // parent_son_flag is to determine whether we are coarsening or refining an element (0 is refining, 1 is unrefining)
{
  int i, j, ifg, l;
  double ndcoord[9][2];
  int Nc = ndof/EQUATIONS;
  double* el_stiffness = new double[Nc*Nc];
  double* el_rhs = new double[Nc];
  Node* NodePtr;
  int which_son2 = which_son;
  
  for(i=0; i<8; i++)
    {
      NodePtr=(Node*)(NodeTable->lookup(node_key[i]));
      assert(NodePtr);
      ndcoord[i][0]=NodePtr->coord[0];      
      ndcoord[i][1]=NodePtr->coord[1];
    }
  
  NodePtr=(Node*) (NodeTable->lookup(key));
  ndcoord[8][0]=NodePtr->coord[0];
  ndcoord[8][1]=NodePtr->coord[1];
  
  Element* EmTemp;
  if(parent_son_flag == 0) 
    EmTemp = (Element*) HT_Elem_Ptr->lookup(getfather());
  
  int Norder[5] = {order, order, order, order, order};
  double* sol_in = new double[Nc*4];
  double* sol_out = new double[Nc];
  double* father_sol0, *father_sol1, *father_sol2, *father_sol3;
  for(l=0;l<2;l++) {
    if(parent_son_flag == 0) {  // big to small element
      if(l==0)
	father_sol0 = EmTemp->get_el_solution();
      else
	father_sol0 = EmTemp->get_prev_el_solution();
    }
    else { //small to big element
      which_son2 = -1;
      if(l==0) {
	EmTemp = (Element*) HT_Elem_Ptr->lookup(son[0]);
	father_sol0 = EmTemp->get_el_solution();
	EmTemp = (Element*) HT_Elem_Ptr->lookup(son[1]);
	father_sol1 = EmTemp->get_el_solution();
	EmTemp = (Element*) HT_Elem_Ptr->lookup(son[2]);
	father_sol2 = EmTemp->get_el_solution();
	EmTemp = (Element*) HT_Elem_Ptr->lookup(son[3]);
	father_sol3 = EmTemp->get_el_solution();
      }
      else { 
	EmTemp = (Element*) HT_Elem_Ptr->lookup(son[0]);
	father_sol0 = EmTemp->get_prev_el_solution();
	EmTemp = (Element*) HT_Elem_Ptr->lookup(son[1]);
	father_sol1 = EmTemp->get_prev_el_solution();
	EmTemp = (Element*) HT_Elem_Ptr->lookup(son[2]);
	father_sol2 = EmTemp->get_prev_el_solution();
	EmTemp = (Element*) HT_Elem_Ptr->lookup(son[3]);
	father_sol3 = EmTemp->get_prev_el_solution();	
      }
    }
    for(i=0;i<EQUATIONS;i++) {
      int no_of_eqns = 1;
      double c = 1; 
      if(i==0 && l == 0)
	ifg = 0; //calculate el_mass
      else 
	ifg = 1;
      if(parent_son_flag == 0) {  // big to small element
	for(j=0;j<Nc;j++)
	  sol_in[j] = father_sol0[j*EQUATIONS+i];
      }
      else { //small to big element
	for(j=0;j<Nc;j++)
	  sol_in[j] = father_sol0[j*EQUATIONS+i];
	for(j=0;j<Nc;j++)
	  sol_in[j+Nc] = father_sol1[j*EQUATIONS+i];
	for(j=0;j<Nc;j++)
	  sol_in[j+2*Nc] = father_sol2[j*EQUATIONS+i];
	for(j=0;j<Nc;j++)
	  sol_in[j+3*Nc] = father_sol3[j*EQUATIONS+i];	
      }
      
#ifdef SUNOS
      masselem_(&ifg, Norder, &Nc, &no_of_eqns, &Nc, 
		&(ndcoord[0][0]), el_stiffness, el_rhs, sol_in, 
		(sol_in+Nc), (sol_in+2*Nc), (sol_in+3*Nc), &which_son2);
      if(i==0)
	tri_(el_stiffness, &Nc, &Nc);
      rhsub_(el_stiffness, el_rhs, &Nc, &Nc, sol_out);
#endif
#ifdef IBMSP
      masselem(&ifg, Norder, &Nc, &no_of_eqns, &Nc, 
	       &(ndcoord[0][0]), el_stiffness, el_rhs, sol_in, 
	       (sol_in+Nc), (sol_in+2*Nc), (sol_in+3*Nc), &which_son2);
      if(i==0)
	tri(el_stiffness, &Nc, &Nc);
      rhsub(el_stiffness, el_rhs, &Nc, &Nc, sol_out);
#endif
      if(l==0) {
	for(j=0;j<Nc;j++)
	  el_solution[j*EQUATIONS+i] = sol_out[j];
      }
      else {
	for(j=0;j<Nc;j++)
	  prev_el_solution[j*EQUATIONS+i] = sol_out[j];
      }
      
    }
  }
  delete []el_stiffness;
  delete []el_rhs;
  delete []sol_in;
  delete []sol_out;
  
  return;
}
/*
  calculates the volume of material in the element
*/
double Element::calc_volume(HashTable* NodeTable) 
{
  int i;
  double ndcoord[4][2], volume;
  Node* NodePtr;
  for(i=4; i<8; i++)
    {
      NodePtr=(Node*)(NodeTable->lookup(node_key[i]));
      assert(NodePtr);
      ndcoord[i-4][0]=NodePtr->coord[0];      
      ndcoord[i-4][1]=NodePtr->coord[1];
    }
  
  
  double d1 = sqrt((ndcoord[2][0]-ndcoord[0][0])*(ndcoord[2][0]-ndcoord[0][0])+
		    (ndcoord[2][1]-ndcoord[0][1])*(ndcoord[2][1]-ndcoord[0][1]));

  double d2 = sqrt((ndcoord[3][0]-ndcoord[1][0])*(ndcoord[3][0]-ndcoord[1][0])+
		    (ndcoord[3][1]-ndcoord[1][1])*(ndcoord[3][1]-ndcoord[1][1]));

  volume = d1*d2*el_solution[0];
  
  return(volume);
}


int Element::find_brothers(HashTable* El_Table, HashTable* NodeTable, 
			   double target, int myid, MatProps* matprops_ptr) 
{
  int i = 0, j;
  int unrefine_flag = 1;
  Element* bros[5];
  if(opposite_brother_flag == 0) {
    find_opposite_brother(El_Table);
    if(opposite_brother_flag == 0) 
      return 0;
  }
  while(i<4 && unrefine_flag == 1) {
    Element* EmTemp = (Element*) El_Table->lookup(&brothers[i][0]);
    if(EmTemp == NULL || EmTemp->refined != 0)
      return 0;
    bros[i+1] = EmTemp;
    unrefine_flag = EmTemp->check_unrefinement(target);
      i++;
  }
  
  if(unrefine_flag == 1) { // we want to unrefine this element...
    //first we create the father element
    //printf("==============================\n unrefining an element \n===============================\n");

    bros[0] = new Element((bros+1), NodeTable, El_Table, matprops_ptr);
    El_Table->add(bros[0]->pass_key(), bros[0]);
    assert(bros[0]);  // a copy of the parent should always be on the same process as the sons
    unrefine_elements(bros, El_Table, NodeTable);
  }
  
  return unrefine_flag;
}
int Element::check_unrefinement(double target) {
  int unrefine_flag = 0, i = 0;
  
  //  put in good element check here!!!
  //for(i=0;i<3;i++) { 
  //if( dabs(el_solution[0]-(*get_prev_el_solution())) > 0.0001*target)
  if(dabs(el_solution[0]) < 2*target  
     //||dabs(el_solution[3]) > 0.001*target || 
     //dabs(el_solution[6]) > 0.001*target*/) {
    && (dabs(el_solution[1])/c_dmax1(dabs(el_solution[0]),5.*target) < 0.1*target ||
        dabs(el_solution[2])/c_dmax1(dabs(el_solution[0]),5.*target) < 0.1*target )) {
    unrefine_flag = 1;
  }
  
  //myprocess=0;//CN oct 27
  //assert(myprocess==0);//CN oct 27

  i = 0;
  while(i<4 && unrefine_flag == 1) {  
    if((neigh_proc[i] != myprocess && neigh_proc[i] >= 0) || neigh_gen[i] > generation)  //maybe check on heigh_gen[i+4]...
      unrefine_flag = 0; 
    //   printf("myprocess%d,myid%d",myprocess,elm->myid);

    i++;
  }

  return unrefine_flag;
}

void Element::find_opposite_brother(HashTable* El_Table)
{
  /* brother information -- requires that atleast one of this
     element's neighboring brothers is on this process in 
     order to get information onthe brother that is not a neighbor */
  Element* EmTemp;
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);  
  int i;
  if(opposite_brother_flag == 1)
    return;
  switch(which_son) {
  case 0:
    if(neigh_proc[2] != -1 && neigh_proc[1] != -1) {
      if(neigh_proc[2] == myid) {
	EmTemp = (Element*) El_Table->lookup(neighbor[2]);
	if(*(EmTemp->get_neigh_gen()+1) == generation) {
	  for(i=0;i<KEYLENGTH;i++)
	    brothers[2][i] = *(EmTemp->get_neighbors()+KEYLENGTH+i);
	  opposite_brother_flag = 1;
      }
	else {
	  EmTemp = (Element*) El_Table->lookup(EmTemp->get_neighbors()+KEYLENGTH);
	  if(EmTemp != NULL && EmTemp->get_refined_flag() != GHOST && EmTemp->get_gen() == (generation+1)) {
	    unsigned* elm_father = EmTemp->getfather();
	    for(i=0;i<KEYLENGTH;i++)
	      brothers[2][i] = elm_father[i];
	    opposite_brother_flag = 1;
	  }
	}
      }
      else if(neigh_proc[1] == myid) {
	EmTemp = (Element*) El_Table->lookup(neighbor[5]);
	if(*(EmTemp->get_neigh_gen()+2) == generation) {
	  for(i=0;i<KEYLENGTH;i++)
	    brothers[2][i] = *(EmTemp->get_neighbors()+2*KEYLENGTH+i);
	  opposite_brother_flag = 1;
      }
	else {
	  EmTemp = (Element*) El_Table->lookup(EmTemp->get_neighbors()+6*KEYLENGTH);
	  if(EmTemp != NULL && EmTemp->get_refined_flag() != GHOST && EmTemp->get_gen() == (generation+1)) {
	    unsigned* elm_father = EmTemp->getfather();
	    for(i=0;i<KEYLENGTH;i++)
	      brothers[2][i] = elm_father[i];
	    opposite_brother_flag = 1;	    
	  }
	}
	
      }
    }
    break;
  case 1:
    if(neigh_proc[2] != -1 && neigh_proc[3] != -1) {
      if(neigh_proc[2] == myid) {
	EmTemp = (Element*) El_Table->lookup(neighbor[6]);
	if(*(EmTemp->get_neigh_gen()+3) == generation) {
	  for(i=0;i<KEYLENGTH;i++)
	    brothers[3][i] = *(EmTemp->get_neighbors()+3*KEYLENGTH+i);
	  opposite_brother_flag = 1;
	}
	else {
	  EmTemp = (Element*) El_Table->lookup(EmTemp->get_neighbors()+7*KEYLENGTH);
	  if(EmTemp != NULL && EmTemp->get_refined_flag() != GHOST && EmTemp->get_gen() == (generation+1)) {
	    unsigned* elm_father = EmTemp->getfather();
	    for(i=0;i<KEYLENGTH;i++)
	      brothers[3][i] = elm_father[i];
	    opposite_brother_flag = 1;
	  }
	}
      }
      else if(neigh_proc[3] == myid) {
	EmTemp = (Element*) El_Table->lookup(neighbor[3]);
	if(*(EmTemp->get_neigh_gen()+2) == generation) {
	  for(i=0;i<KEYLENGTH;i++)
	    brothers[3][i] = *(EmTemp->get_neighbors()+2*KEYLENGTH+i);
	  opposite_brother_flag = 1;
      }
	else {
	  EmTemp = (Element*) El_Table->lookup(EmTemp->get_neighbors()+2*KEYLENGTH);
	  if(EmTemp != NULL && EmTemp->get_refined_flag() != GHOST && EmTemp->get_gen() == (generation+1)) {
	    unsigned* elm_father = EmTemp->getfather();
	    for(i=0;i<KEYLENGTH;i++)
	      brothers[3][i] = elm_father[i];
	    opposite_brother_flag = 1;	    
	  }
	}
	
      }
    }
    break;
  case 2:
    if(neigh_proc[0] != -1 && neigh_proc[3] != -1) {
      if(neigh_proc[0] == myid) {
	EmTemp = (Element*) El_Table->lookup(neighbor[0]);
	if(*(EmTemp->get_neigh_gen()+3) == generation) {
	  for(i=0;i<KEYLENGTH;i++)
	    brothers[0][i] = *(EmTemp->get_neighbors()+3*KEYLENGTH+i);
	  opposite_brother_flag = 1;
	}
	else {
	  EmTemp = (Element*) El_Table->lookup(EmTemp->get_neighbors()+3*KEYLENGTH);
	  if(EmTemp != NULL && EmTemp->get_refined_flag() != GHOST && EmTemp->get_gen() == (generation+1)) {
	    unsigned* elm_father = EmTemp->getfather();
	    for(i=0;i<KEYLENGTH;i++)
	      brothers[0][i] = elm_father[i];
	    opposite_brother_flag = 1;
	  }
	}
      }
      else if(neigh_proc[3] == myid) {
	EmTemp = (Element*) El_Table->lookup(neighbor[3]);
	if(*(EmTemp->get_neigh_gen()) == generation) {
	  for(i=0;i<KEYLENGTH;i++)
	    brothers[0][i] = *(EmTemp->get_neighbors()+i);
	  opposite_brother_flag = 1;
	}
	else {
	  EmTemp = (Element*) El_Table->lookup(EmTemp->get_neighbors());
	  if(EmTemp != NULL && EmTemp->get_refined_flag() != GHOST && EmTemp->get_gen() == (generation+1)) {
	    unsigned* elm_father = EmTemp->getfather();
	    for(i=0;i<KEYLENGTH;i++)
	      brothers[0][i] = elm_father[i];
	    opposite_brother_flag = 1;	    
	  }
	}
	
      }
    }
    break;
  case 3:
    if(neigh_proc[4] != -1 && neigh_proc[1] != -1) {
      if(neigh_proc[4] == myid) {
	EmTemp = (Element*) El_Table->lookup(neighbor[4]);
	if(*(EmTemp->get_neigh_gen()+1) == generation) {
	  for(i=0;i<KEYLENGTH;i++)
	    brothers[1][i] = *(EmTemp->get_neighbors()+KEYLENGTH+i);
	  opposite_brother_flag = 1;
	}
	else {
	  EmTemp = (Element*) El_Table->lookup(EmTemp->get_neighbors()+KEYLENGTH);
	  if(EmTemp != NULL && EmTemp->get_refined_flag() != GHOST && EmTemp->get_gen() == (generation+1)) {
	    unsigned* elm_father = EmTemp->getfather();
	    for(i=0;i<KEYLENGTH;i++)
	      brothers[1][i] = elm_father[i];
	    opposite_brother_flag = 1;
	  }
	}
      }
      else if(neigh_proc[1] == myid) {
	EmTemp = (Element*) El_Table->lookup(neighbor[1]);
	if(*(EmTemp->get_neigh_gen()) == generation) {
	  for(i=0;i<KEYLENGTH;i++)
	    brothers[1][i] = *(EmTemp->get_neighbors()+i);
	  opposite_brother_flag = 1;
	}
	else {
	  EmTemp = (Element*) El_Table->lookup(EmTemp->get_neighbors());
	  if(EmTemp != NULL && EmTemp->get_refined_flag() != GHOST && EmTemp->get_gen() == (generation+1)) {
	    unsigned* elm_father = EmTemp->getfather();
	    for(i=0;i<KEYLENGTH;i++)
	      brothers[1][i] = elm_father[i];
	    opposite_brother_flag = 1;	    
	  }
	}	
      }
    }
    break;
  }

  return;
}

unsigned* Element::getfather()
{
  switch(which_son) {
  case 0:
    return node_key[2];
    break;
  case 1:
    return node_key[3];
    break;
  case 2:
    return node_key[0];
    break;
  case 3:
    return node_key[1];
    break;
  }
  printf("my key is %u %u in getfather on proc %d\n", key[0], key[1], myprocess);
  assert(0); // 0 <= which_son <= 3 !!!
}

/*
using elm_loc, which_son is calculated
*/
void Element::calc_which_son() {
  if(elm_loc[0] %2 == 0) {
    if(elm_loc[1] %2 == 0)
      which_son = 0;
    else
      which_son = 3;
  }
  else  {
    if(elm_loc[1] %2 == 0)
      which_son = 1;
    else
      which_son = 2;
  }

}

/* first approximation is that we calculate the time step according to
   average element values (i.e. the coefficients for the constant shape function) */
double Element::calc_time_step(HashTable* NodeTable, MatProps* matprops_ptr) 
{
  double dt = 0;
  int i, j;

  double xi[2] = {0.,0.};
  double kactxy[2];
  double ndcoords[18];
  double epsilon=(matprops_ptr->HEIGHT_SCALE)/(matprops_ptr->LENGTH_SCALE);
 
  Node* nd;
  for(i=0;i<8;i++) {
    nd = (Node*) NodeTable->lookup(node_key[i]);
    for(j=0;j<2;j++)
      ndcoords[i*2+j] = *(nd->get_coord()+j);
  }
  nd = (Node*) NodeTable->lookup(key);
  for(j=0;j<2;j++)
    ndcoords[16+j] = *(nd->get_coord()+j);

  int ndofe = ndof/EQUATIONS;
  double tiny = GEOFLOW_TINY;
  getkactxy_(&ndofe, &order, xi, (xi+1), ndcoords, &(matprops_ptr->intfrict), 
	     &(matprops_ptr->bedfrict), el_solution, &tiny, kactxy);
  double resolution;
  double xslope,yslope;
  j=Get_max_resolution(&resolution);
  i=Get_slope(resolution,ndcoords[16]*matprops_ptr->LENGTH_SCALE,ndcoords[17]*matprops_ptr->LENGTH_SCALE,&xslope,&yslope);
  double max_slope = sqrt((xslope)*(xslope)+(yslope)*(yslope));
  double max_angle = atan(max_slope); 
  double normal_gravity = 9.8*cos(max_angle)/matprops_ptr->GRAVITY_SCALE;

  if(el_solution[0] > GEOFLOW_TINY) {
    dt = c_dmax1(dabs(el_solution[1]/el_solution[0]), dabs(el_solution[2]/el_solution[0]));
    dt = dt + sqrt(epsilon* normal_gravity* kactxy[0] *el_solution[0]);		
  }
  else {
    dt = GEOFLOW_TINY;
  }
  double dx = c_dmax1(fabs(ndcoords[8] - ndcoords[12]), fabs(ndcoords[9] - ndcoords[13]));
  double dy = c_dmax1(fabs(ndcoords[10] - ndcoords[14]), fabs(ndcoords[11] - ndcoords[15]));
  dx = c_dmin1(dx, dy);
  dt = dx/dt;

  return(dt);
}
