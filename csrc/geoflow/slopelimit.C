#include "../header/hpfem.h"

double min_mod(double, double, double);

void slope_limit(HashTable* El_Table, HashTable* NodeTable, int rkstep, MatProps* matprops_ptr) {
  HashEntryPtr* buck = El_Table->getbucketptr();
  int num_buckets = El_Table->get_no_of_buckets();
  int ii;

  // do slope limiting
  for(ii=0; ii<num_buckets; ii++)
    if(*(buck+ii))
      {
	HashEntryPtr currentPtr = *(buck+ii);
	while(currentPtr)
	  {
	    Element* Curr_El=(Element*)(currentPtr->value);
	    if(*(Curr_El->get_elm_loc()) == 31 && *(Curr_El->get_elm_loc()+1) == 8)
	      rkstep = 0;
	    // only look at active elements with higher order shape functions
	    if(Curr_El->get_refined_flag() == 0 && Curr_El->get_order() > 0)
	      { 
		Curr_El->slopelimit(NodeTable, El_Table, matprops_ptr );
	      }
	    currentPtr=currentPtr->next;
	  }
      }

  return;
}

void Element::slopelimit(HashTable* NodeTable, HashTable* El_Table, MatProps* matprops_ptr) {
  int i, j, k, direction,in;
  if(elm_loc[0] == 13 && elm_loc[1] == 6)
    i = 0;
  assert(order);  //order needs to be greater than 0

   assert(prev_el_solution[0] >= 0.);

   if(prev_el_solution[0] < 0.) {
    printf("element %d %d on process %i has a negative pile height\n",elm_loc[0], elm_loc[1], myprocess);
      /*assert(prev_el_solution[0] >= 0.);*/
    }
  // if pile height is 0 (prev_el_solution[0] = 0) set the momentums to zero...
  if(prev_el_solution[0] < 0) {
    for(i=0;i<ELM_DOF;i++)
      prev_el_solution[i] = 0;
  }

  for(direction=0;direction<2;direction++) 
    if(neigh_proc[1+direction] >= 0 && neigh_proc[3-3*direction] >= 0) {
      Node* np = (Node*) NodeTable->lookup(node_key[5+direction]);
      Node* nn = (Node*) NodeTable->lookup(node_key[7-direction*3]);
      Node* n = (Node*) NodeTable->lookup(key);
      Element* ep = (Element*) El_Table->lookup(neighbor[1+direction]);
      Element* ep2 = NULL;
      if(neigh_gen[1+direction] > generation)
	ep2 = (Element*) El_Table->lookup(neighbor[5+direction]);
      Element* en = (Element*) El_Table->lookup(neighbor[3-direction*3]);
      Element* en2 = NULL;
      if(neigh_gen[3-direction*3] > generation)
	en2 = (Element*) El_Table->lookup(neighbor[7-direction*3]);
      double dx = *(np->get_coord()+direction) - *(nn->get_coord()+direction);
      Node* np_neigh = (Node*) NodeTable->lookup(ep->pass_key());
      double dxp = *(np_neigh->get_coord()+direction) - *(n->get_coord()+direction);
      Node* nn_neigh = (Node*) NodeTable->lookup(en->pass_key());
      double dxn = *(n->get_coord()+direction) - *(nn_neigh->get_coord()+direction);
  
      int order = get_order(); 
      int Nc = calc_num_shape_functions(order)*EQUATIONS;
      int order_n = en->get_order(); 
      int Ncn = calc_num_shape_functions(order_n)*EQUATIONS;
      int order_p = ep->get_order(); 
      int Ncp = calc_num_shape_functions(order_p)*EQUATIONS;

      double xi[2],detj,xi0,xi1;
      double phi[9],dphi[9][2];
      double phin[9];
      double phip[9];
      double phin2[9];
      double phip2[9];
      double gpsi[9],dgpsi[9][2];
      double rjac[2][2],rjacinv[2][2];
      double ndcoord[9][2];
      double U[3],Un[3],Up[3],Un2[3],Up2[3],dUdx[3],dUdy[3];
      double r[3][3],rinv[3][3];

      double myslope2[3],slope_p2[3],slope_n2[3];

      xi[0]=xi[1]=0.;
      xi0=0.;
      xi1=0.;

      for(i=0;i<3;i++)
	for (j=0;j<3;j++){
	  r[i][j]=0.;
	  rinv[i][j]=0.;
	}
      
      //  printf("check%e\n",matprops_ptr->GRAVITY_SCALE);

      gshape_(&xi[0], gpsi, &dgpsi[0][0]);
      // for (i=0;i<9;i++)
      // printf("gpsi[i] %e",gpsi[i]);

      Node* NodePtr;
      for(j=0; j<8; j++)
	{
	  NodePtr=(Node*)(NodeTable->lookup(node_key[j]));
	  //	  NodePtr=(Node*)(getNode()+j*KEYLENGTH);
	  assert(NodePtr);
	  //	  ndcoord[j][0]=NodePtr->coord[0];      
	  ndcoord[j][0]=*(NodePtr->get_coord());      
	  ndcoord[j][1]=*(NodePtr->get_coord()+1);
	}


      NodePtr=(Node*) (NodeTable->lookup(pass_key()));
      ndcoord[8][0]=*(NodePtr->get_coord());
      ndcoord[8][1]=*(NodePtr->get_coord()+1);

      //     dx=sqrt((ndcoord[0][0]-ndcoord[1][0])*(ndcoord[0][0]-ndcoord[1][0])+
      //	      (ndcoord[0][1]-ndcoord[1][1])*(ndcoord[0][1]-ndcoord[1][1]));

      for(i=0;i<2;i++)
	for (j=0;j<2;j++) {
	  rjac[i][j]=0.;
	  for (k=0;k<9;k++){
	    rjac[i][j]=rjac[i][j]+dgpsi[k][j]*ndcoord[k][i];
	    //    printf("j= %d i = %d rjac[j][i]= %e ",i,j,rjac[j][i]);
	  }
	}
      
      //     ...inverse jacobi matrix rjac(i,j) to get rjacinv(i,j)
      detj=rjac[0][0]*rjac[1][1]-rjac[1][0]*rjac[0][1];

      rjacinv[0][0]=rjac[1][1]/detj;
      rjacinv[1][1]=rjac[0][0]/detj;
      rjacinv[1][0]=-rjac[1][0]/detj;
      rjacinv[0][1]=-rjac[0][1]/detj; 

#ifdef SUNOS
      shape2dg_(&order,&xi0,&xi1,phi);
      dshap2dg_(&order,&xi0,&xi1,&dphi[0][0]);
      shape2dg_(&order_n,&xi0,&xi1,phin);
      shape2dg_(&order_p,&xi0,&xi1,phip);
#endif
      int in;
      double vel,kactx,kacty;
      double cosphi,tandel;

      cosphi= cos(matprops_ptr->intfrict); 
      tandel= tan(matprops_ptr->bedfrict);

      // get normal gravity component

      double resolution;
      double xslope,yslope;
      j=Get_max_resolution(&resolution);
      i=Get_slope(resolution,(ndcoord[8][0])*matprops_ptr->LENGTH_SCALE,
		  (ndcoord[8][1])*matprops_ptr->LENGTH_SCALE,&xslope,&yslope);

      double max_slope = sqrt((xslope)*(xslope)+(yslope)*(yslope));
      double max_angle = atan(max_slope); 
      double normal_gravity = 9.8*cos(max_angle)/matprops_ptr->GRAVITY_SCALE;


      for (i=0;i<3;i++){
	  dUdx[i]=0.;
	  dUdy[i]=0.;
      }      

      for (in=0;in<3;in++){
	Up[in]=0.;
	for (i=0;i<Ncp/EQUATIONS;i++)
	  Up[in]=Up[in]+(*(ep->get_prev_el_solution()+i*EQUATIONS+in))*phip[i];
	    }

      if(ep2 != NULL){
	int order_p2=ep2->get_order(); 
	shape2dg_(&order_p2,&xi0,&xi1,phip2);
	int Ncp2 = calc_num_shape_functions(order_p2)*EQUATIONS;
	for (in=0;in<3;in++){
	  Up2[in]=0.;
	  for (i=0;i<Ncp2/EQUATIONS;i++)
	    Up2[in]=Up2[in]+
	      (*(ep2->get_prev_el_solution()+i*EQUATIONS+in))*phip2[i];
	}
      }
      
      for (in=0;in<3;in++){
	U[in]=0.;	
	for (i=0;i<Nc/EQUATIONS;i++){
	  U[in] =U[in]+ (*(get_prev_el_solution()+i*EQUATIONS+in))*phi[i];
	}	
      }

      for (in=0;in<3;in++){
	dUdx[in]=0.;
	dUdy[in]=0.;
	for (i=0;i<Nc/EQUATIONS;i++){
	  dUdx[in]=dUdx[in]+(*(get_prev_el_solution()+i*EQUATIONS+in))*
	    (dphi[i][0]*rjacinv[0][0]+dphi[i][1]*rjacinv[0][1]);
	  dUdy[in]=dUdy[in]+(*(get_prev_el_solution()+i*EQUATIONS+in))*
	    (dphi[i][0]*rjacinv[1][0]+dphi[i][1]*rjacinv[1][1]);
	}	  
      }
      
      if (U[0] > GEOFLOW_TINY){
	vel=dUdx[1]/U[0] - U[1]*dUdx[0]/(U[0]*U[0])+
	  dUdy[2]/U[0] - U[2]*dUdy[0]/(U[0]*U[0]);
      }	  
      else{
	vel=0.;
      }
      
      if (U[0]>GEOFLOW_TINY){
	kactx=(2.0/(cosphi*cosphi))*(1.0-c_sgn(vel)*
	      sqrt(1.0-(1.0+tandel*tandel)*cosphi*cosphi))-1.0;
	kacty=(2.0/(cosphi*cosphi))*(1.0-c_sgn(vel)*
	      sqrt(1.0-(1.0+tandel*tandel)*cosphi*cosphi))-1.0;
	// if there is no yielding...
	if (dabs(U[1]/U[0])< GEOFLOW_TINY ||
	    dabs(U[2]/U[0])< GEOFLOW_TINY) {
	  kactx=1.0;
	  kacty=1.0;
	}
      }
      else {
	vel = 0.0;
	kactx = 1.0;
	kacty = 1.0;
      }
      
      kactx=kactx*matprops_ptr->HEIGHT_SCALE/matprops_ptr->LENGTH_SCALE;
      kacty=kacty*matprops_ptr->HEIGHT_SCALE/matprops_ptr->LENGTH_SCALE;
      
      for (in=0;in<3;in++){
	Un[in]=0.;
	for (i=0;i<Ncn/EQUATIONS;i++)
	  Un[in]=Un[in]+(*(en->get_prev_el_solution()+i*EQUATIONS+in))*phin[i];
      }
      
      if(en2 != NULL){
	int order_n2=en2->get_order();
	shape2dg_(&order_n2,&xi0,&xi1,phin2);
	int Ncn2 = calc_num_shape_functions(order_n2)*EQUATIONS;
	for (in=0;in<3;in++){
	  Un2[in]=0.;
	  for (i=0;i<Ncn2/EQUATIONS;i++)
	    Un2[in]=Un2[in]+
	      (*(en2->get_prev_el_solution()+i*EQUATIONS+in))*phin2[i];
	}
      }

      double slope_p[3],slope_n[3],myslope[3],temp;

	//	dxp=1.;
	//	dxn=1.;
	
      for(i=0;i<3;i++) {
	slope_p[i] = (Up[i] - U[i])/dxp;
	
	if(ep2 != NULL ) {
	  temp = (Up2[i] - U[i])/dxp;
	  /*if(temp*slope_p[i] >0){
	    if(dabs(temp)<dabs(slope_p[i]))
	    slope_p[i] = temp;
	    }  
	    else if (temp*slope_p[i] <0) 
	    slope_p[i]=0;*/
	  slope_p[i]=0.5*(temp+slope_p[i]);
	}
	
	slope_n[i] = (U[i] - Un[i])/dxn;
	
	if(en2 != NULL) {
	  temp = (U[i] - Un2[i])/dxn;
	  /*if(temp*slope_n[i] >0 ){
	    if(dabs(temp)<dabs(slope_n[i]))
	    slope_n[i] = temp;
	    }  
	    else if (temp*slope_n[i] <0) 
	    slope_n[i]=0;*/
	  slope_n[i]=0.5*(temp+slope_n[i]);
	}
	
	 if (direction==0)
	   myslope[i]=dUdx[i];
	 if (direction==1)
	   myslope[i]=dUdy[i];
	
	//	myslope[i]=prev_el_solution[EQUATIONS+direction*EQUATIONS+i];

      }
      
      //  if(elm_loc[0] == 52 && elm_loc[1] == 23)
      //  printf("slopes before limiting %d %e %e \n",i,dUdx[i],dUdy[i]);
      
      if (direction==0) {//begin if loop

	double c=sqrt(kactx*normal_gravity*prev_el_solution[0]);	
	
	if (prev_el_solution[0]> 0.) {
	  
	  rinv[0][0]=(prev_el_solution[1]/prev_el_solution[0]+c)/(2.*c);
	  rinv[0][1]=-1./(2.*c);
	  rinv[0][2]=0.;
	  rinv[1][0]=-prev_el_solution[2]/prev_el_solution[0];
	  rinv[1][1]=0.;
	  rinv[1][2]=1.;
	  rinv[2][0]=(-prev_el_solution[1]/prev_el_solution[0]+c)/(2.*c);
	  rinv[2][1]=1./(2.*c);
	  rinv[2][2]=0.;
	  
	  r[0][0]=1.;
	  r[0][1]=0.;
	  r[0][2]=1.;
	  r[1][0]=prev_el_solution[1]/prev_el_solution[0]-c;
	  r[1][1]=0.;
	  r[1][2]=prev_el_solution[1]/prev_el_solution[0]+c;
	  r[2][0]=prev_el_solution[2]/prev_el_solution[0];
	  r[2][1]=1.;
	  r[2][2]=prev_el_solution[2]/prev_el_solution[0];
	
	}
	
	for (i=0;i<3;i++){
	  myslope2[i]=0.;
	  slope_p2[i]=0.;
	  slope_n2[i]=0.;
	  for (j=0;j<3;j++){
	    myslope2[i]+=rinv[i][j]*myslope[j];
	    slope_p2[i]+=rinv[i][j]*slope_p[j];
	    slope_n2[i]+=rinv[i][j]*slope_n[j];
	  }
	}
	
	for (i=0;i<3;i++)
	  myslope2[i]=min_mod(myslope2[i], slope_p2[i], slope_n2[i]);

	for (i=0;i<3;i++){	    
	  myslope[i]=0.;
	  for (j=0;j<3;j++)
	    myslope[i]+=r[i][j]*myslope2[j];
	}

	//   for (i=0;i<3;i++)
	//   myslope[i]=min_mod(myslope[i], slope_p[i], slope_n[i]);

	double tmp=0;
	for (i=0;i<3;i++)	
	  tmp=tmp+dabs(prev_el_solution[EQUATIONS+direction*EQUATIONS+i]- 
		       myslope[i]/rjacinv[0][0]);

	if (tmp > 0.00000001)
	  for (i=9;i<ELM_DOF;i++)
	    prev_el_solution[i]=0;
	
	for (i=0;i<3;i++){
	  prev_el_solution[EQUATIONS+direction*EQUATIONS+i]=
	    myslope[i]/rjacinv[0][0];	      
	}
      }
      
      if (direction==1){//begin if loop
	double c=sqrt(kactx*normal_gravity*prev_el_solution[0]);
	if (prev_el_solution[0] > 0.){
	  rinv[0][0]=(prev_el_solution[2]/prev_el_solution[0]+c)/(2.*c);
	  rinv[0][1]=0.;
	  rinv[0][2]=-1./(2.*c);
	  rinv[1][0]=-prev_el_solution[1]/prev_el_solution[0];
	  rinv[1][1]=1.;
	  rinv[1][2]=0.;
	  rinv[2][0]=(-prev_el_solution[2]/prev_el_solution[0]+c)/(2.*c);
	  rinv[2][1]=0.;
	  rinv[2][2]=1./(2.*c);
	
	  r[0][0]=1.;
	  r[0][1]=0.;
	  r[0][2]=1.;
	  r[1][0]=prev_el_solution[1]/prev_el_solution[0];
	  r[1][1]=1.;
	  r[1][2]=prev_el_solution[1]/prev_el_solution[0];
	  r[2][0]=prev_el_solution[2]/prev_el_solution[0]-c;
	  r[2][1]=0.;
	  r[2][2]=prev_el_solution[2]/prev_el_solution[0]+c;
	}

	for (i=0;i<3;i++){
	  myslope2[i]=0.;
	  slope_p2[i]=0.;
	  slope_n2[i]=0.;
	  for (j=0;j<3;j++){
	    myslope2[i]+=rinv[i][j]*myslope[j];
	    slope_p2[i]+=rinv[i][j]*slope_p[j];
	    slope_n2[i]+=rinv[i][j]*slope_n[j];
	  }
	}
	
	for (i=0;i<3;i++)
	  myslope2[i]=min_mod(myslope2[i], slope_p2[i], slope_n2[i]);

	for (i=0;i<3;i++){	    
	  myslope[i]=0.;
	  for (j=0;j<3;j++)
	    myslope[i]+=r[i][j]*myslope2[j];
	}
	
	//for (i=0;i<3;i++)
	//	  myslope[i]=min_mod(myslope[i], slope_p[i], slope_n[i]);

	double tmp=0;
	for (i=0;i<3;i++)	
	  tmp=tmp+dabs(prev_el_solution[EQUATIONS+direction*EQUATIONS+i]- 
		       myslope[i]/rjacinv[1][1]);

	if (tmp > 0.00000001)
	  for (i=9;i<ELM_DOF;i++)
	    prev_el_solution[i]=0;	    
	
	for (i=0;i<3;i++){	  	  
	  prev_el_solution[EQUATIONS+direction*EQUATIONS+i]=
	    myslope[i]/rjacinv[1][1];	    
	}	
      }
    }

  // printf("equations %d",EQUATIONS);
  
  // check the four nodes to see if the pile height is less than 0 there...
  double x[2] = {-1.,-1.};
  int Nc = calc_num_shape_functions(order)*EQUATIONS;
  double* f = new double[Nc/EQUATIONS];

  assert(prev_el_solution[0] >= 0.);
  shape2dg_(&order, x, (x+1), f);
  double pile[4] = {0,0,0,0};
  int limitflag = 0;
  for(i=0;i<Nc/EQUATIONS;i++) 
    pile[0] += prev_el_solution[i*EQUATIONS]*f[i];
  
  if(pile[0] < 0.) {
    limitflag = 1;
    pile[0] = 0;
  }
  
  x[0] = 1.;
  shape2dg_(&order, x, (x+1), f);
  for(i=0;i<Nc/EQUATIONS;i++)
    pile[1] += prev_el_solution[i*EQUATIONS]*f[i];
  
  if(pile[1] < 0.) {
    limitflag = 1;
    pile[1] = 0;
  }
      
  x[1] = 1.;
  shape2dg_(&order, x, (x+1), f);
  for(i=0;i<Nc/EQUATIONS;i++)
    pile[2] += prev_el_solution[i*EQUATIONS]*f[i];
  
  if(pile[2] < 0.) {
    limitflag = 1;
    pile[2] = 0;
  }

  x[0] = -1.;
  shape2dg_(&order, x, (x+1), f);
  for(i=0;i<Nc/EQUATIONS;i++)
    pile[3] += prev_el_solution[i*EQUATIONS]*f[i];
  
  if(pile[3] < 0.) {
    limitflag = 1;
    pile[3] = 0;
  }
  
  //   get rid of negative pile height...
   if(limitflag != 0) 
    remove_neg_pile_height(pile, prev_el_solution[0], (prev_el_solution+EQUATIONS), (prev_el_solution+2*EQUATIONS));
  
  x[0] = x[1] = -1;
  shape2dg_(&order, x, (x+1), f);
  double pile2 = 0;
  limitflag = 0;
  for(i=0;i<Nc/EQUATIONS;i++) 
    pile2 += prev_el_solution[i*EQUATIONS]*f[i];
  
  if(pile2 < 0.) {
    limitflag = 1;
  }
  
  x[0] = 1.;
  pile2 = 0;
  shape2dg_(&order, x, (x+1), f); 
  for(i=0;i<Nc/EQUATIONS;i++)
    pile2 += prev_el_solution[i*EQUATIONS]*f[i];
  
  if(pile2 < 0.) {
    limitflag = 1;
  }
      
  pile2 = 0;
  x[1] = 1.;
  shape2dg_(&order, x, (x+1), f);
  for(i=0;i<Nc/EQUATIONS;i++)
    pile2 += prev_el_solution[i*EQUATIONS]*f[i];
  
  if(pile2 < 0.) {
    limitflag = 1;
  }

  pile2 = 0;
  x[0] = -1.;
  shape2dg_(&order, x, (x+1), f);
  for(i=0;i<Nc/EQUATIONS;i++)
    pile2 += prev_el_solution[i*EQUATIONS]*f[i];
  
  if(pile2 < 0.) {
    limitflag = 1;
  }

  delete []f;
  return;
}


double min_mod(double s1, double s2, double s3) {
  double slope = 0.;
  if(s1*s2 < 0 || s1*s3 < 0 || s2*s3 < 0) {
    //printf("set to zero a slope\n");
    return(slope);
  }

  slope = s1;
  if(dabs(s2) < dabs(slope))
    slope = s2;

  if(dabs(s3) < dabs(slope))
    slope = s3;

  /*if(slope != s1)
    printf("reduced the slope\n");*/

  return(slope);
}

void remove_neg_pile_height(double pile_in[], double pile_ave, double *pilex, double *piley) { 
  double pile_height[2];
  pile_height[0] = .25*(pile_in[1]+pile_in[2]-pile_in[0]-pile_in[3]);
  pile_height[1] = .25*(pile_in[2]+pile_in[3]-pile_in[1]-pile_in[0]);
  *pilex = pile_height[0];
  *piley = pile_height[1];
  if((dabs(pile_height[0])+dabs(pile_height[1])) > pile_ave) {
    double reduction = pile_ave/(dabs(pile_height[0])+dabs(pile_height[1])+0.000001*GEOFLOW_TINY);
    *pilex = reduction*pile_height[0];
    *piley = reduction*pile_height[1];
  }

  return;
}
