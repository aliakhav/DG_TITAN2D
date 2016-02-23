#include "../header/hpfem.h"
#define APPLY_BC
#ifdef SUNOS
extern "C" void masselem_(int*, int*, int*, int*, int*, double*,
                          double*, double*, double*, double*, 
			  double*, double*, int*);
extern "C" void elbshal_(int*,double*,int*,int*,double*,double*,
			 double*,double*, double*,
			 double*,double*,double*,int*,double *,
			 double *, double*, double*,
			 double*, double*, double*,int*,int*);
extern "C" void getquad_(int*,double*,double*);
extern "C" void getkactxy_(int *,int *,double *,double *,
			   double *,double*,double*,double*,
			   double*,double*);
#endif	

double calc_time_step(HashTable*, HashTable*, MatProps*);

void slope_limit(HashTable* El_Table, HashTable* NodeTable, 
		 int rkstep, MatProps* matprops_ptr);

void update_el_solution(HashTable* El_Table);

void store_temp_el_solution(HashTable* El_Table);

void store_el_solution(HashTable* El_Table);

void flux_compute(double*, Element*, HashTable*, HashTable*,
		  int ,int ,double*, double*, double*,
		  double*, double*, double );

void flux_calcx(double* h, double* hp,double* grav,double* gravn, 
		double* kactxy, double* kactxyn, double* flux);
void flux_calcy(double* h, double* hp,double* grav,double* gravn, 
		double* kactxy, double* kactxyn, double* flux);

void step(HashTable* El_Table, HashTable* NodeTable, int myid, int nump,
	  double end_time, double* time_ptr, MatProps* matprops_ptr, 
	  int timestep, double* maxFluxIntegral) 
{
  //printf("========================================================\n");
    
  /* get coefficients, eigenvalues, hmax and calculate the time step */
  double dt = calc_time_step(El_Table, NodeTable, matprops_ptr);

  if (timestep<100)
    dt=dt*timestep/100;


  if(end_time < dt + *time_ptr)
    dt = end_time - *time_ptr;



  printf("for iteration %d time is %e and time step is %e\n",
	 timestep, *time_ptr, dt);

  *time_ptr += dt;
  double dt2 = .5*dt; // dt2 is set as dt/2 !

  int ii;
  if(timestep == 90)
    ii = 1;
  int i,j,k, counter;
  double tiny = GEOFLOW_TINY;
  double tmp=0.;
  double ndcoord[9][2];
  int nequ=EQUATIONS;
  //-------------------go through all the elements of the subdomain and 
  //  compute mass matrix, "element matrix(uprev*uprev)", fluxes, sources 
  //-------------------

  HashEntryPtr* buck = El_Table->getbucketptr();
  int num_buckets = El_Table->get_no_of_buckets();
  int num_rkstep, rkstep;  //runge-kutta step
  update_el_solution(El_Table); //make sure el_solution and prev_el_solution
                                //dont have negative pile height
  move_data(nump, myid, El_Table, NodeTable);

  if(ORDER == 0)
    num_rkstep = 1;
  else if(ORDER >= 1)
        num_rkstep = 2;
  // num_rkstep = 1;
  for(rkstep=0;rkstep<num_rkstep;rkstep++) {
    move_data(nump, myid, El_Table, NodeTable);
    // if (rkstep==0)
    for(ii=0; ii<num_buckets; ii++)
      if(*(buck+ii))
	{
	  HashEntryPtr currentPtr = *(buck+ii);
	  while(currentPtr)
	    {
	      Element* Curr_El=(Element*)(currentPtr->value);
	      if(Curr_El->get_refined_flag() == 0)
		{
		  //Node* nd = (Node*) NodeTable->lookup(Curr_El->pass_key());

		  if(*(Curr_El->get_elm_loc()) == 31 && 
		     *(Curr_El->get_elm_loc()+1) == 7)
		    i = 0;
		  
		  //get node coordinate
		  Node* NodePtr;
		  for(j=0; j<8; j++)
		    {
		      //NodePtr=(Node*)(NodeTable->lookup(node_key[j]));
		      NodePtr=(Node*)(NodeTable->lookup(Curr_El->getNode()+j*KEYLENGTH));
		      assert(NodePtr);
		      //ndcoord[j][0]=NodePtr->coord[0];      
		      ndcoord[j][0]=*(NodePtr->get_coord());      
		      ndcoord[j][1]=*(NodePtr->get_coord()+1);
		    }
		  NodePtr=(Node*) (NodeTable->lookup(Curr_El->pass_key()));
		  ndcoord[8][0]=*(NodePtr->get_coord());
		  ndcoord[8][1]=*(NodePtr->get_coord()+1);
		  int order=Curr_El->get_order(); 
		  int Nc = calc_num_shape_functions(order)*EQUATIONS;
		  
		  double * sol_in = new double[Nc*4];
		  double* el_mass = new double[Nc*Nc];
		  double* el_stiffness = new double[Nc*Nc];
		  double* el_rhs = new double[Nc];		
		  int which_son2 =-3;
		  int no_of_eqns=EQUATIONS;

		  for (i=0;i<4*Nc;i++)
		    sol_in[i]=0.;
		  
		  for (i=0;i<Nc;i++)
		    el_rhs[i]=0.;
		  for (j=0;j<Nc*Nc;j++)
		    {
		      el_stiffness[j]=0.;
		      el_mass[j]=0;
		    }
		  
		  int ifg=0;


#ifdef SUNOS
		  masselem_(&ifg, &order, &Nc, &no_of_eqns, &Nc, 
			    &(ndcoord[0][0]), el_mass, el_rhs, sol_in, 
			    (sol_in+Nc), (sol_in+2*Nc), (sol_in+3*Nc), 
			    &which_son2);
#endif

		  double* solptr = Curr_El->get_el_solution();
		  
		  for(j=0;j<Nc;j++)
		    sol_in[j] = *(solptr+j);
	
		  int ndofe=Nc/EQUATIONS;
		  //
		  // set gravity/curvature using linear interpolation of gravities corresponding ot nodal interpolation
		  //
		  double resolution, x, y,gravnd[12];
		  double  xslope,yslope,curvature[8];
		  j=Get_max_resolution(&resolution);
		  
		  for (j=0;j<4;j++)
		    {
		      x = ndcoord[j][0];
		      y = ndcoord[j][1];
		      i=Get_slope(resolution,
				  x*matprops_ptr->LENGTH_SCALE,
				  y*matprops_ptr->LENGTH_SCALE,
				  &xslope,&yslope);
		      i=Get_curvature(resolution, 
				      x*matprops_ptr->LENGTH_SCALE, 
				      y*matprops_ptr->LENGTH_SCALE, 
				      (curvature+j*2), 
				      (curvature+j*2+1));
		      
		      double max_slope = 
			sqrt((xslope)*(xslope)+(yslope)*(yslope));
		      //if (sol_in[0]>0)
		      //	printf("maxslope %e\n",max_slope);
		      double max_angle = atan(max_slope); 
		      double down_slope_gravity = 
			9.8*sin(max_angle)/matprops_ptr->GRAVITY_SCALE;
		      if(dabs(down_slope_gravity) > GEOFLOW_TINY) {
			gravnd[0+j*3] = -down_slope_gravity*(xslope)/max_slope;
			gravnd[1+j*3] = -down_slope_gravity*(yslope)/max_slope;
			gravnd[2+j*3] = 9.8*cos(max_angle)/
			  matprops_ptr->GRAVITY_SCALE;
		      }
		      else {
			gravnd[0+j*3] = 0;
			gravnd[1+j*3] = 0;
			gravnd[2+j*3] = 9.8/matprops_ptr->GRAVITY_SCALE;
		      }
		  }
		  for (i=0;i<8;i++)
		    curvature[i] = curvature[i] * (matprops_ptr->LENGTH_SCALE);
		  double* tempstorage = new double[3*ndofe+Nc*(Nc+2)]; 
		  //only used for temp storage in elbshal


#ifdef SUNOS
		  elbshal_(&Nc,&(ndcoord[0][0]),&order,&ndofe,sol_in,
			   el_stiffness,el_rhs,&tiny,gravnd,curvature, 
			   &(matprops_ptr->bedfrict), 
			   &(matprops_ptr->intfrict),&nequ,
			   &(matprops_ptr->LENGTH_SCALE),
			   &(matprops_ptr->HEIGHT_SCALE),tempstorage, 
			   (tempstorage+ndofe), 
			   (tempstorage+3*ndofe), 
			   (tempstorage+3*ndofe+Nc),
			   (tempstorage+3*ndofe+2*Nc),&rkstep,&timestep); 
		  
#endif
		  delete []tempstorage;
		  // compute fluxes and add
		  double * nodecord=&(ndcoord[0][0]);
		  double epsilon= (matprops_ptr->HEIGHT_SCALE)/
		    (matprops_ptr->LENGTH_SCALE);
		  double el_rhs2[9];
		  flux_compute(el_rhs, Curr_El,El_Table,NodeTable,
			       order,Nc,sol_in, nodecord,
			       gravnd, &(matprops_ptr->bedfrict), 
			       &(matprops_ptr->intfrict),epsilon);

		  //update maximum of the flux integral over the entire domain

		  tmp=c_dmax1(tmp, Curr_El->get_fluxIntegral());

		  //update maximum of the flux

		  // update
		  double* predictor=new double[Nc];
		  for(j=0;j<Nc;j++)
		    predictor[j]=0.;
		  double* solout= new double[Nc];
		  for (j=0;j<Nc;j++)
		    {
		      solout[j]=el_rhs[j];
		      el_rhs[j]=dt*el_rhs[j];
		    }
		  /* apply bc's */
#ifdef APPLY_BC
		  for(j=0;j<4;j++)
		    if(*(Curr_El->get_neigh_proc()+j) == INIT) 
		      {
			if (sol_in[0]>0){
			  printf("applying BC %d,%d\n",
				 *(Curr_El->get_elm_loc()),
				 *(Curr_El->get_elm_loc()+1));
			  for (k=0;k<Nc;k++)
			    {
			      el_rhs[k]=0.;
			    }
			}
		      }// this is a boundary!
#endif


#ifdef SUNOS
		  tri_(el_mass, &Nc, &Nc);
		  rhsub_(el_mass,el_rhs,&Nc,&Nc,predictor);
#endif


		  //intermediate Solution (predictor)--stored to temp_el_sol

		    for(j=0;j<Nc;j++)
		      predictor[j]= sol_in[j] + predictor[j];
		    Curr_El->store_temp_el_solution(predictor,Nc);
		  
		  delete []el_mass;
		  delete []el_rhs;
		  delete []el_stiffness;
		  delete []solout;
		  delete []predictor;
		  delete []sol_in;
		  
		} // end of if->refined
	      currentPtr=currentPtr->next;    
	    } 
	}// end of element loop

    move_data(nump, myid, El_Table, NodeTable);
    store_temp_el_solution(El_Table); //copy temp_el_solution to el_solution  
    slope_limit(El_Table, NodeTable, rkstep, matprops_ptr);
    double gl_maxFluxIntegral;
    i = MPI_Allreduce(&tmp, &gl_maxFluxIntegral, 1, 
		   MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    *maxFluxIntegral=gl_maxFluxIntegral;
    // printf(" maxFluxIntegral %e in step.C \n ", *maxFluxIntegral);
    //  exit(0);
 
   } //ends rk step loop

  /* finished step */
 

  //  HashEntryPtr* buck = El_Table->getbucketptr();
  // int num_buckets = El_Table->get_no_of_buckets();
  // int ii;

  for(ii=0; ii<num_buckets; ii++)
    if(*(buck+ii))
      {
	HashEntryPtr currentPtr = *(buck+ii);
	while(currentPtr)
	  {
	    Element* Curr_El=(Element*)(currentPtr->value);
	    int order=Curr_El->get_order();
	    int Nc = calc_num_shape_functions(order)*EQUATIONS; 
	    double* prevsol = Curr_El->get_prev_el_solution();
	    double* solptr = Curr_El->get_el_solution();
	    int* elm_loc=Curr_El->get_elm_loc();
	    double * sol_in = new double[Nc*4];

	    for(j=0;j<Nc;j++)
	      sol_in[j] = 0.5*(*(solptr+j)+*(prevsol+j));
	      //sol_in[j] = *(solptr+j);//for euler's

	    if( elm_loc[0]==11 && (elm_loc[1]== 17 || elm_loc[1]== 12)) {
	      printf("(in step.C) elmloc %d %d\n",elm_loc[0],elm_loc[1]);
	      for(j=0;j<Nc;j++)
		printf("solution[%d]= %e  ",j,sol_in[j]);
	      printf("(in step.C)\n");
	    }

	    if(sol_in[0] < 0.) {
	      printf("(in step.C) element %d %d has a negative pile height\n",
		     elm_loc[0], elm_loc[1]);
	      printf("myid sol_in, prev_sol %d %e %e",myid, *(solptr),*(prevsol));
	      assert(sol_in[0] >= 0.);
	    }
	    
	    if(Curr_El->get_refined_flag() == 0)
	      Curr_El->store_el_solution(sol_in,Nc);
	    delete []sol_in;
	    currentPtr = currentPtr->next;
	  }
      }
  
  calc_volume(El_Table, NodeTable,nump,myid, matprops_ptr);

  return;
}

void flux_compute(double* el_rhs, Element* Curr_El,HashTable* El_Table,
		  HashTable* NodeTable,int order,int Nc,double* sol_in, 
		  double *  ndcoord, double* gravnd, double* bedfriction, 
		  double* intfriction,double epsilon)
{
  // loop over edges
  // get my solutions and  neighbor solutions
  // go to FORTRAN to get quadrature point and shape function values 
  // run a loop over quadrature points -- evaluate fluxes and update appropriate cpt of el_rhs
  
  int i,j,nint,k,l;
  double xi[10],wa[10],phi[10],nphi[10];
  double xil,eta,xin,etan,rtx,rty,rnx,rny,ds,Xnod[3][2];
  double h[EQUATIONS],hp[EQUATIONS];
  Element* neigh_elm;
  Element* neigh_elm2;
  double flux_int[3][4]; //three components for h,uh,vh four sides 0,1,2,3
  double tmp=0;
  unsigned* neighbors,neighbor_on_this_side ;
  int*           neigh_gen;
  int            mygen; 
  int*           neighproc;
  int*           elm_loc;
  double *       sol_neig1=new double[Nc];
  double        dphi[20];
  double grav[3],gravn[3],fluxx[3]={0,0,0},fluxy[3]={0,0,0},flux[3]={0,0,0};
  double kactxy[2],kactxyn[2],ndcoordn[9][2];
  double tiny=GEOFLOW_TINY;
  order= Curr_El->get_order();

  neighbors= Curr_El->get_neighbors();
  neigh_gen = Curr_El->get_neigh_gen();
  mygen     = Curr_El->get_gen();
  neighproc    = Curr_El->get_neigh_proc();
  elm_loc=Curr_El->get_elm_loc();


  for (k=0;k<3;k++)
    for (l=0;l<4;l++)
      flux_int[k][l]=0.;

  tmp=0.;
  
  for (i=0;i<4;i++)
    if(neighproc[i] >= 0) {//only do flux calc if there is a neighbor on this side 
      
      if (i<=2)
	{
	  Xnod[0][0] = ndcoord [i*2];
	  Xnod[0][1] = ndcoord [i*2+1];
	  Xnod[2][0] = ndcoord[(i+1)*2];
	  Xnod[2][1] = ndcoord[(i+1)*2+1];
	}
        if (i==2)
	{
	   Xnod[0][0] = ndcoord[3*2];
	   Xnod[0][1] = ndcoord[3*2+1];
	   Xnod[2][0] = ndcoord[2*2];
	   Xnod[2][1] = ndcoord[2*2+1];
	   }
      if (i==3)
	{
	  Xnod[2][0] = ndcoord[0];
	  Xnod[2][1] = ndcoord[1];
	  Xnod[0][0] = ndcoord[3*2];
	  Xnod[0][1] = ndcoord[3*2+1]; 
	}
        
      Xnod[1][0] = ndcoord[(i+4)*2];
      Xnod[1][1] = ndcoord[(i+4)*2+1];

      neigh_elm = (Element*) El_Table->lookup(Curr_El->get_neighbors()+i*KEYLENGTH);
      neigh_elm2 = (Element*) El_Table->lookup(Curr_El->get_neighbors()+(i+4)*KEYLENGTH);

      //number of integration points is controlled by the polynomial order of 
      // the 2(or 3) elements sharing the edge
      if(order >= neigh_elm->get_order())
	nint = order+2;
      else
	nint = neigh_elm->get_order()+2;
      // if this is a big element next to 2 small elements, increase the number of integration
      // points and check the other neighbor polynomial order
      if(mygen < neigh_gen[i]) {
	Element* elm = (Element*) El_Table->lookup(Curr_El->get_neighbors()+(4+i)*KEYLENGTH);
	if(order < elm->get_order() && neigh_elm->get_order() < elm->get_order())
	  nint = elm->get_order()+2;
	nint++;
      }
      
      /* use only even amounts of quadrature points so that 
	 numerical integration is never done at x=0 or y=0 */
      if(nint%2 == 1)
	nint++;

      	nint=4;

#ifdef SUNOS
      getquad_(&nint,xi,wa);
#endif
      //      int jj;
      //     for (jj=0;jj<4;jj++)
      //	  if (*(Curr_El->get_elm_loc())==53 && *(Curr_El->get_elm_loc()+1)==22 && i==2)
      //	      printf(" xij %e\n",xi[jj]); 


        if (Curr_El->get_gen()>neigh_elm->get_gen())
      	nint=2;//NCC 10 Nov 2003



      //  printf("order %d,nint %d,xi %e\n",order,nint,xi[1]); 
      //  printf("elmloc %d %d,nint %d,xi %e\n",elm_loc[1],elm_loc[2],nint,xi[1]); 
      for (j=0;j<nint;j++)
	{
	  switch(i) {
	  case 0:
	    eta=-1.;
	    xil=xi[j];
	    break;
	  case 1:
	    xil=1.;
	    eta=xi[j];
	    break;
	  case 2:
	    eta=1.;
	    xil=xi[j];
	    break;
	  case 3:
	    xil=-1.;
	    eta=xi[j];
	    break;
	  }

	  //   int jj;
	  //   for (jj=0;jj<4;jj++)
	  //	  if (*(Curr_El->get_elm_loc())==53 && *(Curr_El->get_elm_loc()+1)==22 && i==2)
	  //	      printf(" xijj %e\n",xi[jj]); 


	  
	  /*
	   * WARNING:  the code below REQUIRES that the mesh is a mapped mesh with 
	   *           each element having the same orientation (e.g. each element's
	   *           'zero' node is the bottom left node)
	   */
	  int neigh_side=neigh_elm->which_neighbor(Curr_El->pass_key());
	  if(Curr_El->get_gen() == neigh_elm->get_gen()) {  
	    // same generation
	    switch(neigh_side) {
	    case 0:
	      etan=-1.;
	      xin=xi[j];
	      break;
	    case 1:
	      xin=1.;
	      etan=xi[j];
	      break;
	    case 2:
	      etan=1.;
	      xin=xi[j];
	      break;
	    case 3:
	      xin=-1.;
	      etan=xi[j];
	      break;
	    }
	  }
	  
	  if(Curr_El->get_gen() < neigh_elm->get_gen()) {
	    // neighbor is smaller than I am ("big to small")
 //first find out which neighbor to use
	    //side 0
	    if( i==0 && xil<0 ) {
	      neigh_elm = 
		(Element*) El_Table->lookup(Curr_El->get_neighbors()+i*KEYLENGTH);
	      neigh_side=Curr_El->which_neighbor(neigh_elm->pass_key());
	    }
	    else if( i==0 && xil>0 ) {
	      neigh_elm = 
		(Element*) El_Table->lookup(Curr_El->get_neighbors()+(i+4)*KEYLENGTH);
	      neigh_side=Curr_El->which_neighbor(neigh_elm->pass_key());
	    }
	    //side 1
	    if( i==1 && eta < 0.) {
	      neigh_elm = 
		(Element*) El_Table->lookup(Curr_El->get_neighbors()+i*KEYLENGTH);
	      neigh_side=Curr_El->which_neighbor(neigh_elm->pass_key());      
	    }
	    else if( i==1 && eta > 0.) {
	      neigh_elm = 
		(Element*) El_Table->lookup(Curr_El->get_neighbors()+(i+4)*KEYLENGTH);
	      neigh_side=Curr_El->which_neighbor(neigh_elm->pass_key());      
	    }
	    //side 2	
	    if( i==2 && xil > 0) {
	      neigh_elm = 
		(Element*) El_Table->lookup(Curr_El->get_neighbors()+i*KEYLENGTH);
	      neigh_side=Curr_El->which_neighbor(neigh_elm->pass_key());	
	    }
	    else if( i==2 && xil < 0) {
	      neigh_elm = 
		(Element*) El_Table->lookup(Curr_El->get_neighbors()+(i+4)*KEYLENGTH);
	      neigh_side=Curr_El->which_neighbor(neigh_elm->pass_key());	
	    }	
	    //side 3
	    if( i==3 && eta > 0.) {
	      neigh_elm = 
		(Element*) El_Table->lookup(Curr_El->get_neighbors()+i*KEYLENGTH);
	      neigh_side=Curr_El->which_neighbor(neigh_elm->pass_key());
	    }
	    else if( i==3 && eta < 0.) {
	      neigh_elm = 
		(Element*) El_Table->lookup(Curr_El->get_neighbors()+(i+4)*KEYLENGTH);
	      neigh_side=Curr_El->which_neighbor(neigh_elm->pass_key());
	    }
	    //end "find out which neighbor to use"	   

	    //	    if (*(Curr_El->get_elm_loc())==52 && *(Curr_El->get_elm_loc()+1)==23 && i==1)
	    //	      printf("neigh_side %d, neigh_elm_loc %d %d integration point %d xil %e\n",neigh_side,*(neigh_elm->get_elm_loc()),*(neigh_elm->get_elm_loc()+1),j,xil); 

	    int nint1=2;double xi1[10];
	    getquad_(&nint1,xi1,wa);
	    double xij;
	    switch(j){
	    case 0:
	      xi1[0]=xi1[0];
	      break;
	    case 1:
	      xi1[1]=xi1[1];
	      break;
	    case 2:
	      xi1[2]=xi1[0];
	      break;
	    case 3:
	      xi1[3]=xi1[1];
	      break;
	    }

	    switch(j){
	    case 0:
	      wa[0]=0.5*wa[0];
	      break;
	    case 1:
	      wa[1]=0.5*wa[1];
	      break;
	    case 2:
	      wa[2]=0.5*wa[0];
	      break;
	    case 3:
	      wa[3]=0.5*wa[1];
	      break;
	    }
	      
	    switch(neigh_side) {
	    case 0:
	      //  neigh_elm = (Element*) El_Table->lookup(Curr_El->get_neighbors()+2*KEYLENGTH);
	      eta=-1.;
	      xil=0.5*xi1[j]-0.5;
	      etan=1.;
	      xin=xi1[j];
	      break;
	    case 1:
	      //  neigh_elm = (Element*) El_Table->lookup(Curr_El->get_neighbors()+3*KEYLENGTH);
	      xil=1.;
	      eta=0.5*xi1[j]-0.5;
	      xin=-1.;
	      etan=xi1[j];
	      break;
	    case 2:
	      //  neigh_elm = (Element*) El_Table->lookup(Curr_El->get_neighbors()+0*KEYLENGTH);
	      eta=1.;
	      xil=0.5*xi1[j]+0.5;
	      etan=-1.;
	      xin=xi1[j];
	      break;
	    case 3:
	      //  neigh_elm = (Element*) El_Table->lookup(Curr_El->get_neighbors()+1*KEYLENGTH);
	      xil=-1.;
	      eta=0.5*xi1[j]+0.5;
	      xin=1.;
	      etan=xi1[j];
	      break;
	    case 4:
	      //  neigh_elm = (Element*) El_Table->lookup(Curr_El->get_neighbors()+6*KEYLENGTH);
	      eta=-1.;
	      xil=0.5*xi1[j]+0.5;
	      etan=1.;
	      xin=xi1[j];
	      break;
	    case 5:
	      //  neigh_elm = (Element*) El_Table->lookup(Curr_El->get_neighbors()+7*KEYLENGTH);
	      xil=1.;
	      eta=0.5*xi1[j]+0.5;
	      xin=-1.;
	      etan=xi1[j];
	      break;
	    case 6:
	      //  neigh_elm = (Element*) El_Table->lookup(Curr_El->get_neighbors()+4*KEYLENGTH);
	      eta=+1.;
	      xil=0.5*xi1[j]-0.5;
	      etan=-1.;
	      xin=xi1[j];
	      break;
	    case 7:
	      // neigh_elm = (Element*) El_Table->lookup(Curr_El->get_neighbors()+5*KEYLENGTH);
	      xil=-1.;
	      eta=0.5*xi1[j]-0.5;
	      xin=1.;
	      etan=xi1[j];
	      break;
	    }
	  }
	  // end ("big to small")

	  if (Curr_El->get_gen() > neigh_elm->get_gen()) {
	    // neighbor is bigger than I am (" small to big")
	    // first figure out which map to use...


	    int nint1=2;double xi1[2];double wa1[2];
	    getquad_(&nint1,xi1,wa);
	    double xij;
	    switch(j){
	    case 0:
	      xi1[0]=xi1[0];
	      break;
	    case 1:
	      xi1[1]=xi1[1];
	      break;
	    case 2:
	      xi1[2]=xi1[0];
	      break;
	    case 3:
	      xi1[3]=xi1[1];
	      break;
	    }
	    /*	    switch(j){
	    case 0:
	      wa[0]=wa[0];
	      break;
	    case 1:
	      wa[1]=wa[1];
	      break;
	    case 2:
	      wa[2]=wa[0];
	      break;
	    case 3:
	      wa[3]=wa[1];
	      break;
	      }*/


	    //	    int jj;
	    //	    for (jj=0;jj<4;jj++)
	    //	      printf("xi= %e \n",xi1[j]);

	    

	    //side 0
	    if( i==0 && neigh_elm->which_neighbor(Curr_El->pass_key())<4 ) {
	      xil=xi1[j];
	      eta=-1.;
	      xin=0.5*xi1[j]+0.5;
	      etan=1.;
	      //	      printf("neighbor % d \n",neigh_elm->which_neighbor(Curr_El->pass_key()));
	    }
	    else if( i==0 && neigh_elm->which_neighbor(Curr_El->pass_key())>=4)  {
	      xil=xi1[j];
	      eta=-1.;
	      xin=0.5*xi1[j]-0.5;
	      etan=1.;
	      //	      printf("neighbor % d \n",neigh_elm->which_neighbor(Curr_El->pass_key()) );
	    }

	    //side 1
	    if( i==1 && neigh_elm->which_neighbor(Curr_El->pass_key())<4) {
	      xil=1.;
	      eta=xi1[j];
	      xin=-1.;
	      etan=0.5*xi1[j]+0.5;
	      //	      printf("neighbor % d \n",neigh_elm->which_neighbor(Curr_El->pass_key()));
	    }
	    else if( i==1 && neigh_elm->which_neighbor(Curr_El->pass_key())>=4 )  {
	      xil=1.;
	      eta=xi1[j];
	      xin=-1.;
	      etan=0.5*xi1[j]-0.5;
	      //	      printf("neighbor % d \n",neigh_elm->which_neighbor(Curr_El->pass_key()));
	    }

	    //side 2	
	    if( i==2 && neigh_elm->which_neighbor(Curr_El->pass_key())<4 ) {
	      xil=xi1[j];
	      eta=1.;
	      xin=0.5*xi1[j]-0.5;
	      etan=-1.;
	      //	      printf("neighbor % d \n",neigh_elm->which_neighbor(Curr_El->pass_key()));
	    }
	    else if( i==2 && neigh_elm->which_neighbor(Curr_El->pass_key())>=4)  {
	      xil=xi1[j];
	      eta=1.;
	      xin=0.5*xi1[j]+0.5;
	      etan=-1.;
	      //	      printf("neighbor % d \n",neigh_elm->which_neighbor(Curr_El->pass_key()));
	    }	

	    //side 3
	    if( i==3 && neigh_elm->which_neighbor(Curr_El->pass_key()) <4) {
	      xil=-1.;
	      eta=xi1[j];
	      xin=1.;
	      etan= 0.5*xi1[j]-0.5;
	      //	      printf("neighbor % d \n",neigh_elm->which_neighbor(Curr_El->pass_key()));
	    }
	    else if ( i==3 && neigh_elm->which_neighbor(Curr_El->pass_key())>=4 )  {
	      xil=-1.0;
	      eta=xi1[j];
	      xin=1.;
	      etan=0.5*xi1[j]+0.5;
	      //	      printf("neighbor % d \n",neigh_elm->which_neighbor(Curr_El->pass_key()));
	    }

	    //end "find out which neighbor to use"
	  }
	  
	  // get neighbor solution, gravity at quad point
	  double* neigh_solptr;
	  
	  // got to find which side of neighbor to get correct quadrature point in neighbor element
	  neigh_solptr= neigh_elm->get_el_solution();
	  for(k=0;k<Nc;k++)
	    sol_neig1[k] = neigh_solptr[k];
	  
	  int neigh_order=neigh_elm->get_order();
	  
	  //	  if (fabs(xin)>1 || fabs(etan)>1) 
	  //printf("bad qaudrature %u %d %d %d %e, %e %e %e\n",Curr_El->pass_key(),i,neigh_side,j,xil,eta,xin,etan);

	  int nord = 1+Curr_El->get_order(); 
#ifdef SUNOS
	  shape2dg_(&order,&xil,&eta,phi);
	  // dshap2dg_(&order,&xil,&eta,dphi);
#endif
	  
	  // get my solution at quadrature point
	  h[0]=0.;h[1]=0.;h[2]=0.;
	  hp[0]=0.;hp[1]=0.;hp[2]=0.;
	  for (k=0;k<Nc/EQUATIONS;k++)
	    {
	      h[0]= h[0]+sol_in[k*EQUATIONS]*phi[k];
	      h[1]= h[1]+sol_in[k*EQUATIONS+1]*phi[k];
	      h[2]= h[2]+sol_in[k*EQUATIONS+2]*phi[k];
	    }
	  
#ifdef SUNOS	  
	  shape2dg_(&neigh_order,&xin,&etan,nphi);
#endif

	  for (k=0;k<Nc/EQUATIONS;k++)
	    {
	      hp[0]= hp[0]+sol_neig1[k*EQUATIONS]*nphi[k];
	      hp[1]= hp[1]+sol_neig1[k*EQUATIONS+1]*nphi[k];
	      hp[2]= hp[2]+sol_neig1[k*EQUATIONS+2]*nphi[k];
	    }
	  
	  if (h[0]<-tiny ||hp[0]<-tiny) {
	    //printf("trouble (h<0) h hp %e %e side %d neigh side %d, quadrature %e %e qudrature neighbor %e %e, %e %e %e\n",h[0],hp[0],i,neigh_side,xin,etan,Xnod[0][0],Xnod[0][1],nphi[0],nphi[1],nphi[2]);
	  }
	  
	  for(k=0;k<3;k++) 
	    grav[k]=0;
	  
	  double phi_grav[4];  // using integrated hierarchical Legendre shape functions for gravity
	  int order_grav[5] = {1,1,1,1,1};
	  shape2_(order_grav, &xil, &eta, phi_grav);
	  for(k=0;k<4;k++)
	    {
	      grav[0]=grav[0]+gravnd[k*3]*phi_grav[k];
	      grav[1]=grav[1]+gravnd[k*3+1]*phi_grav[k];
	      grav[2]=grav[2]+gravnd[k*3+2]*phi_grav[k];
	    }
	  for(k=0;k<3;k++){gravn[k]=grav[k];}
	  
	  // get kactxy and neighbor kactxy
	  int ndofe=Nc/EQUATIONS;
	  getkactxy_(&ndofe,&order,&xil,&eta,ndcoord,intfriction,bedfriction,&(sol_in[0]),&tiny,kactxy);
	  kactxy[0]=kactxy[0]*epsilon;
	  kactxy[1]=kactxy[1]*epsilon;
	  
	  Node* NodePtr;
	  for(k=0; k<8; k++)
	    {
	      NodePtr=(Node*)(NodeTable->lookup(neigh_elm->getNode()+k*KEYLENGTH));
	      assert(NodePtr);
	      ndcoordn[k][0]=*(NodePtr->get_coord());      
	      ndcoordn[k][1]=*(NodePtr->get_coord()+1);
	    }
	  
	  NodePtr=(Node*) (NodeTable->lookup(neigh_elm->pass_key()));
	  ndcoordn[8][0]=*(NodePtr->get_coord());
	  ndcoordn[8][1]=*(NodePtr->get_coord()+1);
	  
	  double* ndcoordp = &(ndcoordn[0][0]);
	  getkactxy_(&ndofe,&neigh_order,&xin,&etan,ndcoordp,intfriction,bedfriction,sol_neig1,&tiny,kactxyn);
	  kactxyn[0]=kactxyn[0]*epsilon;
	  kactxyn[1]=kactxyn[1]*epsilon;
	  

	  if (hp[0]>tiny || h[0]>tiny) {
	    if(i==0 || i == 3) {
	      flux_calcx(hp,h,gravn,grav,kactxyn,kactxy,fluxx);
	      flux_calcy(hp,h,gravn,grav,kactxyn,kactxy,fluxy);
	    }
	    else {
	      flux_calcx(h,hp,grav,gravn,kactxy,kactxyn,fluxx);
	      flux_calcy(h,hp,grav,gravn,kactxy,kactxyn,fluxy);
	    }
	  }
	  else 
	    for(k=0;k<3;k++) 
	      {
		fluxx[k] = 0;
		fluxy[k] = 0;
	      }
	  //normals
	  rtx = Xnod[0][0]*(xi[j]-0.5)-Xnod[1][0]*2.*xi[j]+
	    Xnod[2][0]*(xi[j]+0.5);
	  rty = Xnod[0][1]*(xi[j]-0.5)-Xnod[1][1]*2.*xi[j]+
	    Xnod[2][1]*(xi[j]+0.5);
	  ds = sqrt(rtx*rtx+rty*rty);

	  rtx = rtx/ds;
	  rty = rty/ds;
	  rnx = rty;
	  rny = -rtx;

	  double signx,signy;
	  switch(i) 
	    {
	    case 0:
	      signx=0;
	      signy=1;
	      break;
	    case 1:
	      signx=-1;
	      signy=0;
	      break;
	    case 2:
	      signx=0;
	      signy=-1;
	      break;
	    case 3:
	      signx=1;
	      signy=0;
	    }

	  if (h[0]>tiny || hp[0]>tiny) 
	    for (k=0;k<3;k++) 
	      {
		flux[k]=ds*wa[j]*(fluxx[k]*signx+fluxy[k]*signy);
	      }
	  else { 
	    flux[0]=0.;
	    flux[1]=0;
	    flux[2]=0;
	  }
	  
	  // accumulate
	  int in;

	  if (h[0]>tiny || hp[0]>tiny)
	    {
	      for (k=0;k<Nc/EQUATIONS;k++) 
	       	{ 
		  for(in=0;in<3;in++)
		    { 
		      el_rhs[k*EQUATIONS+in]=
			el_rhs[k*EQUATIONS+in]+phi[k]*flux[in];
		      flux_int[in][i]= flux_int[in][i]+ dabs(phi[k]*flux[in]);
	  //flux_int is the integral of fluxes -- for debugging purpose 
		    }
	      	}
	    }
	}// end of integration point loop
       //calculate maximum of the integrated flux 
      tmp=tmp+dabs(flux_int[0][i])+dabs(flux_int[1][i])+dabs(flux_int[2][i]);
    }

  Curr_El->put_fluxIntegral(tmp);
  
  delete []sol_neig1;
  return;
}


void flux_calcx(double* h, double* hp,double* grav,double* gravn, 
		double* kactxy, double* kactxyn, double* flux)
  /* Function to calculate Davis fluxes give edge states of element and neighbor
     h -> my edge state, hp -> my neighbor state, kactxy -> my k active/passive, kactxyn -> neighbor k active/passive
     grav -> gravity 
  */
{
  double a,ap,s,sp;
  if((hp[0] <= GEOFLOW_TINY)  && (h[0] <= GEOFLOW_TINY )) {
    flux[0] = 0;
    flux[1] = 0;
    flux[2] = 0;
  }
  
  if ((hp[0] > GEOFLOW_TINY)  && (h[0] > GEOFLOW_TINY) ) {
    a=sqrt(kactxy[1]*grav[2]*h[0]);
    ap=sqrt(kactxyn[1]*gravn[2]*hp[0]);
    s= c_dmin1(0.,c_dmin1(h[1]/h[0]-a,hp[1]/hp[0]-ap));
    sp= c_dmax1(0.,c_dmax1(h[1]/h[0]+a,hp[1]/hp[0]+ap));
    
    flux[0]=((sp*h[1]-s*hp[1]+s*sp*(hp[0]-h[0]))/(sp-s));
    flux[1]=(sp*(h[1]*h[1]/h[0]+0.5*h[0]*a*a)
	     -s*(hp[1]*hp[1]/hp[0]+0.5*hp[0]*ap*ap)+
	     s*sp*(hp[1]-h[1]))/(sp-s);
    flux[2]=(sp*h[1]*h[2]/h[0]-s*hp[1]*hp[2]/hp[0]+
	     s*sp*(hp[2]-h[2]))/(sp-s);
  }
  
  if (hp[0] > GEOFLOW_TINY  && h[0] < GEOFLOW_TINY ) {
    a=0;
    ap=sqrt(kactxyn[1]*gravn[2]*hp[0]);
    h[0]=0,h[1]=0,h[2]=0;
    s= c_dmin1(0,c_dmin1(hp[1]/hp[0]-2*ap,hp[1]/hp[0]-ap));
    sp= c_dmax1(0,c_dmax1(hp[1]/hp[0]-2*ap,hp[1]/hp[0]+ap));
    
    flux[0]=((sp*h[1]-s*hp[1]+s*sp*(hp[0]-h[0]))/
	     (sp-s));
    flux[1]=((sp*(0+0.5*h[0]*a*a)
	      -s*(hp[1]*hp[1]/hp[0]+0.5*hp[0]*ap*ap)+
	      s*sp*(hp[1]-h[1]))/(sp-s));
    flux[2]=((0-s*hp[1]*hp[2]/hp[0]+
	      s*sp*(hp[2]-h[2]))/(sp-s));
  }
  
  if (hp[0] < GEOFLOW_TINY  && h[0] > GEOFLOW_TINY ) {
    a=sqrt(kactxy[1]*grav[2]*h[0]);
    ap=0;
    hp[0]=0,hp[1]=0,hp[2]=0;
    s= c_dmin1(0,c_dmin1(h[1]/h[0]-a,h[1]/h[0]+2*a));
    sp= c_dmax1(0,c_dmax1(h[1]/h[0]+a,h[1]/h[0]+2*a));
    
    flux[0]=((sp*h[1]-s*hp[1]+s*sp*(0-h[0]))/
	     (sp-s));
    flux[1]=((sp*(h[1]*h[1]/h[0]+0.5*h[0]*a*a)
	      -s*(0+0.5*hp[0]*ap*ap)+
	      s*sp*(hp[1]-h[1]))/(sp-s));
    flux[2]=((sp*h[1]*h[2]/h[0]+
	      s*sp*(hp[2]-h[2]))/(sp-s));
  }
  return;
}

void flux_calcy(double* h, double* hp,double* grav,double* gravn, 
		double* kactxy, double* kactxyn, double* flux)
  /*
    Function to calculate Davis fluxes give edge states of element and neighbor
  h -> my edge state, hp -> my neighbor state, kactxy -> my k active/passive, kactxyn -> neighbor k active/passive
  grav -> gravity 
  */
{
  double a,ap,s,sp;
  
  if((hp[0] <= GEOFLOW_TINY)  && (h[0] <= GEOFLOW_TINY )) {
    flux[0] = 0;
    flux[1] = 0;
    flux[2] = 0;
  }
  
  if ((hp[0] > GEOFLOW_TINY)  && (h[0] > GEOFLOW_TINY) ) {
    a=sqrt(kactxy[1]*grav[2]*h[0]);
    ap=sqrt(kactxyn[1]*gravn[2]*hp[0]);
    s= c_dmin1(0,c_dmin1(h[2]/h[0]-a,hp[2]/hp[0]-ap));
    sp= c_dmax1(0,c_dmax1(h[2]/h[0]+a,hp[2]/hp[0]+ap));
    
    flux[0]=((sp*h[2]-s*hp[2]+s*sp*(hp[0]-h[0]))/
	     (sp-s));
    flux[1]=((sp*h[1]*h[2]/h[0]-s*hp[1]*hp[2]/hp[0]+
	      s*sp*(hp[1]-h[1]))/(sp-s));
    flux[2]=((sp*(h[2]*h[2]/h[0]+0.5*h[0]*a*a)
	      -s*(hp[2]*hp[2]/hp[0]+0.5*hp[0]*ap*ap)+
	      s*sp*(hp[2]-h[2]))/(sp-s));
  }
  
  if (hp[0] > GEOFLOW_TINY  && h[0] < GEOFLOW_TINY ) {
    a=0;
    ap=sqrt(kactxyn[1]*gravn[2]*hp[0]);
    h[0]=0,h[1]=0,h[2]=0;
    s= c_dmin1(0,c_dmin1(hp[2]/hp[0]-2*ap,hp[2]/hp[0]-ap));
    sp= c_dmax1(0,c_dmax1(hp[2]/hp[0]-2*ap,hp[2]/hp[0]+ap));
    
    flux[0]=((sp*h[2]-s*hp[2]+s*sp*(hp[0]-h[0]))/
	     (sp-s));
    flux[1]=((0-s*hp[1]*hp[2]/hp[0]+
	      s*sp*(hp[1]-h[1]))/(sp-s));
    flux[2]=((sp*(0+0.5*h[0]*a*a)
	      -s*(hp[2]*hp[2]/hp[0]+0.5*hp[0]*ap*ap)+
	      s*sp*(hp[2]-h[2]))/(sp-s));
  }
  
  if (hp[0] < GEOFLOW_TINY  && h[0] > GEOFLOW_TINY ) {
    a=sqrt(kactxy[1]*grav[2]*h[0]);
    ap=0;
    hp[0]=0,hp[1]=0,hp[2]=0;
    s= c_dmin1(0,c_dmin1(h[2]/h[0]-a,h[2]/h[0]+2*a));
    sp= c_dmax1(0,c_dmax1(h[2]/h[0]+a,h[2]/h[0]+2*a));
    
    flux[0]=((sp*h[2]-s*hp[2]+s*sp*(hp[0]-h[0]))/
	     (sp-s));
    flux[1]=((sp*h[1]*h[2]/h[0]-0+
	      s*sp*(hp[1]-h[1]))/(sp-s));
    flux[2]=((sp*(h[2]*h[2]/h[0]+0.5*h[0]*a*a)
	      -s*(0+0.5*hp[0]*ap*ap)+
	      s*sp*(hp[2]-h[2]))/(sp-s));
  }
  return;
}

void update_el_solution(HashTable* El_Table) {
  // first go through all of the elements and copy the current solution (el_solution)
  // into the previous solution (el_prev_solution)
  HashEntryPtr* buck = El_Table->getbucketptr();
  int num_buckets = El_Table->get_no_of_buckets();
  int ii;
  for(ii=0; ii<num_buckets; ii++)
    if(*(buck+ii))
      {
	HashEntryPtr currentPtr = *(buck+ii);
	while(currentPtr)
	  {
	    Element* Curr_El=(Element*)(currentPtr->value);
 	    if(Curr_El->get_refined_flag() == 0)
	      Curr_El->update_el_solution();
	    currentPtr = currentPtr->next;
	  }
      }

  return;
}

void store_temp_el_solution(HashTable* El_Table) {
  // first go through all of the elements and copy the current solution (el_solution)
  // into the previous solution (el_prev_solution)
  HashEntryPtr* buck = El_Table->getbucketptr();
  int num_buckets = El_Table->get_no_of_buckets();
  int ii;
  for(ii=0; ii<num_buckets; ii++)
    if(*(buck+ii))
      {
	HashEntryPtr currentPtr = *(buck+ii);
	while(currentPtr)
	  {
	    Element* Curr_El=(Element*)(currentPtr->value);
 	    if(Curr_El->get_refined_flag() == 0)
	      Curr_El->store_temp_el_solution();
	    currentPtr = currentPtr->next;
	  }
      }

  return;
}


//given an order of approximation, calculate the number of shape
//functions needed to achieve that order
int calc_num_shape_functions(int order) {
  if(order == 0)
    return(1);
  else if(order == 1)
    return(3);
  else if(order ==2)
    return(6);
  else {
    printf("can only do orders of 0, 1, or 2\n");
    assert(0);
  }
    
}


double calc_time_step(HashTable* El_Table, HashTable* NodeTable, MatProps* matprops_ptr) {
  double timestep = .0001;
  HashEntryPtr* buck = El_Table->getbucketptr();
  int num_buckets = El_Table->get_no_of_buckets();
  int ii, i;
  double proc_info[2];

  proc_info[0]= 999999999999.;
  proc_info[1]=0.;

  for(ii=0; ii<num_buckets; ii++)
    if(*(buck+ii))
      {
	HashEntryPtr currentPtr = *(buck+ii);
	while(currentPtr)
	  {
	    Element* Curr_El=(Element*)(currentPtr->value);
	    if(*(Curr_El->get_elm_loc()) == 16 && *(Curr_El->get_elm_loc()+1) == 7)
	      int rkstep = 0;
	    if(Curr_El->get_refined_flag() == 0)
	      { 
		double dt = Curr_El->calc_time_step(NodeTable, matprops_ptr);
		if(dt < proc_info[0])
		  proc_info[0] = dt;
		if(proc_info[1] < *(Curr_El->get_el_solution()))
		  proc_info[1] = *(Curr_El->get_el_solution());
	      }
	    currentPtr=currentPtr->next;
	  }
      }


 
  proc_info[0] = 0.5*proc_info[0];
  proc_info[1] = -0.9 * proc_info[1];

  double gl_proc_info[2];
  i = MPI_Allreduce(proc_info, gl_proc_info, 2, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

/*  if((-gl_proc_info[1]) < gl_proc_info[0])
    printf("the time step is being set by the pile height\n");
  else
    printf("the time step is being set by the eigenvalues\n");
    
   printf("maximum pile height is %e \n",(gl_proc_info[1]/(-0.9))*matprops_ptr->HEIGHT_SCALE);   
*/
    timestep = 0.005 * c_dmin1(gl_proc_info[0], -gl_proc_info[1]);
  
   //timestep =0.00005;

  return(timestep);
}

