#include "../header/hpfem.h"
#include "../header/blas.h"

void mat_vec(double*, double*, double*, double*, int*, int*, 
	     int*, int*, int, int, int, int*);
void mat_vec(double*, double*, double*, int*, int*, int*, int*, int, int);

//diagonal preconditioner and storing all of search vectors on all processors
void gmres_sol(int myid, int subdof, int totdof, int bb, int Bb, 
	       HashTable* BT_Node_Ptr, HashTable* BT_Elem_Ptr, 
	       int num_blocks, int* rpntr, 
	       int* bpntr, int* bindx, int* indx, double* load, 
	       double* stiff, int* ix, double* u)

{
  if(myid ==0)
    printf("in gmres with global vector storage \n");
  int i=0, j=0, m1=0, m2=0;
  double rm_zm, p_A_p, rmnew_zmnew, alpha, beta, sum;
  double tol = .0000000001;
  double* y;
  
  //  convergence parameters
  double zm0 = 1000000;
  double reduction = .000001;  //reduce preconditioned residual zm0 by reduction
  double resid = zm0*2+1;
  int restart_steps = 50;  //amount of gmres iterations before restart
  if(restart_steps >= totdof)
    restart_steps = totdof - 1; 
  int max_iter = 100;
  double resid2;
  int one = 1;

  for(i=0;i<totdof;i++)  u[i] = 0;

  //givens rotation variables
  double nu = 0;
  double* rcos = new double[restart_steps];
  double* rsin = new double[restart_steps];
  //setup GMRES
  double start = MPI_Wtime();
  //create diagonal preconditioner
  //vel. prec.
  int numrows, k;
  double* store = new double[totdof];
  for(i=0;i<totdof;i++)  *(store+i) = 0;
  for(i=0;i<num_blocks;i++)  //goes through all block rows
    {
      numrows = *(rpntr+i+1) - *(rpntr+i);  //number of rows in this block
      for(j=*(bpntr+i);j<*(bpntr+i+1);j++) //find diag block in this row
	if(*(bindx+j) == i)	  
	  for(k=0;k<numrows;k++)
	    *(store+*(ix+k+*(rpntr+i))) = *(stiff+*(indx+j)+k+k*numrows);
    } 

 
  double* diag_glob = new double[totdof];
  i=MPI_Allreduce(store, diag_glob, totdof, MPI_DOUBLE,
                  MPI_SUM,MPI_COMM_WORLD);  //sum diagonals
  
  for(i=0;i<totdof;i++)  
    *(diag_glob+i) = 1/(*(diag_glob+i));
  //done with diag prec
  
  //get global rhs
  for(i=0;i<totdof;i++)  store[i] = 0;
  for(i=0;i<subdof;i++)
    store[ix[i]] = load[i];
  double* rhs = new double[totdof];
  i=MPI_Allreduce(store, rhs, totdof, MPI_DOUBLE,
                  MPI_SUM,MPI_COMM_WORLD);  
  delete []store;
  
  double* gg = new double[restart_steps+1];
  double* v = new double[(totdof)*(restart_steps+1)];  //orthogonal vectors are stored in rows
  double* h = new double[(1+restart_steps)*restart_steps];
  // start of gmres
  m1 = 0; 

  while(m1<max_iter && resid >zm0*reduction)
    {
      //calculate preconditioned residual
      store = new double[totdof]; 
      mat_vec(u, v, store, stiff, rpntr, bpntr, bindx,
	      indx, num_blocks, subdof, totdof, ix); 
      delete []store;

      for(i=0;i<totdof;i++)
	v[i] = rhs[i] - v[i];

      sum = 0;
      for(i=0;i<totdof;i++)
	{
	  v[i] = v[i]*diag_glob[i];
	  sum += v[i]*v[i];
	}
      gg[0] = sqrt(sum);
      sum = 1/gg[0];
      for(i=0;i<totdof;i++)
	v[i] = v[i]*sum;


      resid = gg[0];
      if(m1 == 0)
	zm0 = resid;  //initial residual
      m2 = 0;
/*      if(myid==0)
	printf("outer loop residual is %e zm0 is %e\n",resid, zm0); */
      while(m2< restart_steps && resid > zm0*reduction)
	{
	  store = new double[totdof];
	  mat_vec((v+m2*totdof), (v+(m2+1)*totdof), store, stiff, rpntr, bpntr, bindx,
		  indx, num_blocks, subdof, totdof, ix); 
	  delete []store;

	  for(i=0;i<totdof;i++)
	    v[totdof*(m2+1)+i] = v[totdof*(m2+1)+i]*diag_glob[i];

	  for(i=0;i<=m2;i++)
	    {
#ifdef SUNOS
	      h[i*restart_steps+m2] = ddot_(&totdof, (v+totdof*(m2+1)), &one, (v+i*totdof), &one);
#endif
#ifdef IBMSP
	      h[i*restart_steps+m2] = ddot(&totdof, (v+totdof*(m2+1)), &one, (v+i*totdof), &one);
#endif
#ifdef CRAY
	      h[i*restart_steps+m2] = DDOT(&totdof, (v+totdof*(m2+1)), &one, (v+i*totdof), &one);
#endif
	      for(j=0;j<totdof;j++)
		v[totdof*(m2+1)+j] = v[totdof*(m2+1)+j] - h[i*restart_steps+m2]*v[j+i*totdof];
	    }
#ifdef SUNOS
	  h[(m2+1)*restart_steps+m2] = ddot_(&totdof, (v+totdof*(m2+1)), &one,(v+totdof*(m2+1)), &one);
#endif
#ifdef IBMSP
	  h[(m2+1)*restart_steps+m2] = ddot(&totdof, (v+totdof*(m2+1)), &one,(v+totdof*(m2+1)), &one);
#endif
#ifdef CRAY
	  h[(m2+1)*restart_steps+m2] = DDOT(&totdof, (v+totdof*(m2+1)), &one,(v+totdof*(m2+1)), &one);
#endif
	  h[(m2+1)*restart_steps+m2] = sqrt(h[(m2+1)*restart_steps+m2]);
	  double denom = 1/h[(m2+1)*restart_steps+m2];
	  for(i=0;i<totdof;i++)
	    v[totdof*(m2+1)+i] = v[totdof*(m2+1)+i]*denom;
	  
	  //Givens rotation
	  for(i=0;i<m2;i++)  //finish Givens rotation from previous iterations
	    {
	      sum = h[i*restart_steps+m2];
	      h[i*restart_steps+m2] = rcos[i]*sum+rsin[i]*h[(i+1)*restart_steps+m2];
	      h[(i+1)*restart_steps+m2] = -rsin[i]*sum+rcos[i]*h[(i+1)*restart_steps+m2];
	    }
	  nu = sqrt(pow(h[m2*restart_steps+m2],2)+pow(h[(m2+1)*restart_steps+m2],2));
	  rcos[m2] = h[m2*restart_steps+m2]/nu; 
	  rsin[m2] = h[(m2+1)*restart_steps+m2]/nu;  
	  h[m2*restart_steps+m2] = rcos[m2]*h[m2*restart_steps+m2] + rsin[m2]*h[(m2+1)*restart_steps+m2]; 
	  h[(m2+1)*restart_steps+m2] = 0;
	  gg[m2+1] = -gg[m2]*rsin[m2];
	  gg[m2] = gg[m2]*rcos[m2];
	  	  
	  resid = fabs(gg[m2+1]);
	  m2++;
/*	  if(myid ==0)
	    printf("inner loop residual is %e zm0 is %e m2 is %d \n",resid, zm0, m2); */
	  

	  //checks for orthogonality of v vectors
/*	  for(i=0;i<m2;i++)
	    {
#ifdef SUNOS
	      sum = ddot_(&totdof, (v+totdof*m2), &one, (v+i*totdof), &one);
#endif
#ifdef IBMSP
	      sum = ddot(&totdof, (v+totdof*m2), &one, (v+i*totdof), &one);
#endif
#ifdef CRAY
	      sum = DDOT(&totdof, (v+totdof*m2), &one, (v+i*totdof), &one);
#endif
	      if(sum > .0001)
		printf("nonorthogonal!!!  sum = %e \n", sum);
	    }*/

	}
      y = new double[m2];
      
      for(i=m2-1;i>-1;i--)
	{
	  for(j=m2-1;j>i;j--)
	    gg[i] = gg[i] - h[i*restart_steps+j]*y[j];
	  y[i] = gg[i]/h[i*restart_steps+i];
	}
	    
	
      for(i=0;i<m2;i++)
	for(j=0;j<totdof;j++)
	  u[j] += v[j+i*totdof]*y[i];


      //checks for accuracy of the residual
/*      store = new double[totdof];
      mat_vec(u, v, store, stiff, rpntr, bpntr, bindx,
	      indx, num_blocks, subdof, totdof, ix); 
      delete []store;

      sum = 0;
      for(i=0;i<totdof;i++)
	{
	  v[i] = (rhs[i] - v[i])*diag_glob[i];
	  sum +=v[i]*v[i];
	}
      sum = sqrt(sum);
      printf("real preconditioned residual is %e \n",sum);*/
      //done checking residual accuracy

      delete []y;
      m1++;
    }
  if(myid ==0)
    if(m1 == max_iter)
      printf("gmres did not converge!!!!!!!!!! \n");
  delete []gg;
  delete []v;
  delete []h;
  double end = MPI_Wtime();  //timer
  end = end -start;

  
  if(myid == 0)
    {  
      printf("diagonal preconditioner\n");
      printf("iter = %d", m2+(m1-1)*restart_steps);
      printf("   resid = %e\n", resid);
      printf("corner and edge dof = %d\n", totdof);
      printf("   solution time = %e\n", end); 
    }
  MPI_Barrier(MPI_COMM_WORLD); 

  delete []diag_glob;
  delete []rhs;
  delete []rcos;
  delete []rsin;

  delete []stiff;
  delete []load;
  delete []rpntr;
  delete []bpntr;
  delete []bindx;
  delete []indx;  
  //reconstruct bubble functions solution
  start = MPI_Wtime();
  Assemble_bubble(totdof, subdof, Bb, bb, u, ix, BT_Elem_Ptr, BT_Node_Ptr, myid);
  end = MPI_Wtime();  //timer
  end = end -start;
  printf("bubble reconstruction time = %e for proc %d\n", end, myid);
  return;
}


//**********************************************************************************
//diagonal preconditioner and storing local portion of search vectors on each processor
void gmres_sol(int subdof, int totdof, int bb, int Bb,
	       int* ix, int nump, int myid, double* u,
	       double* load, double* stiff, int* rpntr, int* bpntr, 
	       int* bindx, int* indx, int num_blocks, HashTable* BT_Node_Ptr,
	       HashTable* BT_Elem_Ptr, int NNSS, int nnss, int vvee)
{
  if(myid ==0)
    printf("in gmres with local vector storage \n");
  int i=0, j=0, m1=0, m2=0;
  double rm_zm, p_A_p, rmnew_zmnew, alpha, beta, sum;
  double tol = .0000000001;
  
  //  convergence parameters
  double zm0 = 1000000;
  double reduction = .000001;  //reduce preconditioned residual zm0 by reduction
  double resid = 2*zm0+1;
  int restart_steps = 50;  //amount of gmres iterations before restart
  if(restart_steps >= totdof)
    restart_steps = totdof - 1;
  int max_iter = 100;
  double resid2;
  int one = 1;

  for(i=0;i<totdof;i++)  u[i] = 0;

  //givens rotation variables
  double nu = 0;
  double* rcos = new double[restart_steps];
  double* rsin = new double[restart_steps];

  //setup which dof belong to which proc(used for global inner products)
  int* loc2 = new int[NNSS];
  for(i=0;i<NNSS;i++)
    loc2[i]= nump;
  
  for(i=0;i<nnss;i++)
    loc2[ix[i]] = myid;

  int* loc = new int[NNSS];

  i=MPI_Allreduce(loc2, loc, NNSS, MPI_INT, MPI_MIN,MPI_COMM_WORLD);
  delete []loc2;
  char* c_loc = new char[nnss];  //char is used to minimize memory usage(only stores 0 & 1)
  for(i=0;i<nnss;i++)
    {
      c_loc[i] = 0;
      if(loc[ix[i]] == myid)
	c_loc[i] = 1;  //if c_loc[i] == 1, dof belongs to this proc
    }
  delete []loc;

  //setup GMRES
  double start = MPI_Wtime();

  //create diagonal preconditioner for this processors subdomain
  int numrows, k;
  double* diag_subdomain = new double[subdof];
  for(i=0;i<subdof;i++)  *(diag_subdomain+i) = 0;
  for(i=0;i<num_blocks;i++)  //goes through all block rows
    {
      numrows = *(rpntr+i+1) - *(rpntr+i);  //number of rows in this block
      for(j=*(bpntr+i);j<*(bpntr+i+1);j++) //find diag block in this row
	if(*(bindx+j) == i)	  
	  for(k=0;k<numrows;k++)
	    *(diag_subdomain+k+*(rpntr+i)) = *(stiff+*(indx+j)+k+k*numrows);
    } 
  int array_length = NNSS; 
  if(NNSS < (restart_steps+1))  //array_length should usually be NNSS but for 1 proc will be restart
    array_length = restart_steps+1;

  double* store2 = new double[array_length]; 
  double* store = new double[array_length];
  for(i=0;i<NNSS;i++)  store[i] = 0;
  for(i=0;i<nnss;i++)
    store[ix[i]] = diag_subdomain[i];

  i=MPI_Allreduce(store, store2,NNSS, MPI_DOUBLE,
                  MPI_SUM,MPI_COMM_WORLD);  //sum interface diagonals
  for(i=0;i<nnss;i++)
    diag_subdomain[i] = 1/store2[ix[i]];  //diag is in subdomain ordering
  for(i=nnss;i<subdof;i++)  
    *(diag_subdomain+i) = 1/(*(diag_subdomain+i));
  //done with diagonal prec

  //make sure that each processor has global value for the rhs(load)
  for(i=0;i<NNSS;i++)  store[i] = 0;
  for(i=0;i<nnss;i++)
    store[ix[i]] = load[i];

  i=MPI_Allreduce(store, store2,NNSS, MPI_DOUBLE,
                  MPI_SUM,MPI_COMM_WORLD);  //sum interface
  for(i=0;i<nnss;i++)
    load[i] = store2[ix[i]];  //load is in subdomain ordering
  //done with global sum for rhs

  double* gg = new double[restart_steps+1];
  double* v = new double[(subdof)*(restart_steps+1)];  //orthogonal vectors are stored in rows
  double* h = new double[(1+restart_steps)*restart_steps];
  // start of gmres
  m1 = 0; 

  while(m1<max_iter && resid >zm0*reduction)
    {
      //put global u into local store
      double* u_local = new double[subdof];
      for(i=0;i<subdof;i++)
	u_local[i] = u[ix[i]];
      
      //matrix vector multiply
      mat_vec(u_local, v, stiff, rpntr, bpntr, bindx,
	      indx, num_blocks, subdof); 
      delete []u_local;
      for(i=0;i<NNSS;i++)  store[i] = 0;
      for(i=0;i<nnss;i++)
	store[ix[i]] = v[i];
      i=MPI_Allreduce(store, store2,NNSS, MPI_DOUBLE,
		      MPI_SUM,MPI_COMM_WORLD);  //sum interface
      for(i=0;i<nnss;i++)
	v[i] = store2[ix[i]];  //v is in subdomain ordering

      for(i=0;i<subdof;i++)
	v[i] = load[i] - v[i];

      sum = 0;
      for(i=0;i<nnss;i++)
	{
	  v[i] = v[i]*diag_subdomain[i];
	  if(c_loc[i] ==1)  //dof belongs to this proc
	    sum += v[i]*v[i];
	}
      for(i=nnss;i<subdof;i++)
	{
	  v[i] = v[i]*diag_subdomain[i];
	  sum += v[i]*v[i];
	}
      double sum2;
      i=MPI_Allreduce(&sum, &sum2, 1,MPI_DOUBLE,
		      MPI_SUM,MPI_COMM_WORLD);

      gg[0] = sqrt(sum2);
      sum = 1/gg[0];
      for(i=0;i<subdof;i++)
	v[i] = v[i]*sum;

      resid = gg[0];
      if(m1 == 0)
	zm0 = resid;  //initial residual
      m2 = 0;
/*      if(myid==0)
	printf("outer loop residual is %e zm0 is %e m1 is %d proc %d\n",resid, zm0, m1, myid);*/
      while(m2< restart_steps && resid > zm0*reduction)
	{
	  mat_vec((v+m2*subdof), (v+(m2+1)*subdof), stiff, rpntr, bpntr, bindx,
		  indx, num_blocks, subdof); 

	  for(i=0;i<NNSS;i++)  store[i] = 0;
	  for(i=0;i<nnss;i++)
	    store[ix[i]] = v[i+(m2+1)*subdof];
	  i=MPI_Allreduce(store, store2,NNSS, MPI_DOUBLE,
			  MPI_SUM,MPI_COMM_WORLD);  //sum interface
	  for(i=0;i<nnss;i++)
	    v[i+(m2+1)*subdof] = store2[ix[i]];  //v is in subdomain ordering	  

	  for(i=0;i<subdof;i++)
	    v[subdof*(m2+1)+i] = v[subdof*(m2+1)+i]*diag_subdomain[i];

	  //right so far
	  for(i=0;i<=m2;i++)
	    {
#ifdef SUNOS
	      store[i] = ddot_(&vvee, (v+subdof*(m2+1)+nnss), &one, (v+i*subdof+nnss), &one);
#endif
#ifdef IBMSP
	      store[i] = ddot(&vvee, (v+subdof*(m2+1)+nnss), &one, (v+i*subdof+nnss), &one);
#endif
#ifdef CRAY
	      store[i] = DDOT(&vvee, (v+subdof*(m2+1)+nnss), &one, (v+i*subdof+nnss), &one);
#endif
	      for(j=0;j<nnss;j++)
		if(c_loc[j] == 1)
		  store[i] += v[subdof*(m2+1)+j]* v[i*subdof+j];
	    }
	  i = MPI_Allreduce(store, store2, m2+1,MPI_DOUBLE,
			    MPI_SUM,MPI_COMM_WORLD);
	  for(i=0;i<=m2;i++)
	    h[i*restart_steps+m2] = store2[i];

	  for(i=0;i<=m2;i++)
	    {
	      for(j=0;j<subdof;j++)
		v[subdof*(m2+1)+j] = v[subdof*(m2+1)+j] - h[i*restart_steps+m2]*v[j+i*subdof];
	    }
#ifdef SUNOS	  
	  sum = ddot_(&vvee,(v+subdof*(m2+1)+nnss),&one,(v+subdof*(m2+1)+nnss),&one);
#endif
#ifdef IBMSP	  
	  sum = ddot(&vvee,(v+subdof*(m2+1)+nnss),&one,(v+subdof*(m2+1)+nnss),&one);
#endif
#ifdef CRAY
	  sum = DDOT(&vvee,(v+subdof*(m2+1)+nnss),&one,(v+subdof*(m2+1)+nnss),&one);
#endif
	  for(i=0;i<nnss;i++)
	    if(c_loc[i] ==1)
	      sum += v[subdof*(m2+1)+i]*v[subdof*(m2+1)+i];
	  i = MPI_Allreduce(&sum, (h+(m2+1)*restart_steps+m2), 1,MPI_DOUBLE,
			    MPI_SUM,MPI_COMM_WORLD);	  
	  h[(m2+1)*restart_steps+m2] = sqrt(h[(m2+1)*restart_steps+m2]);
	  double denom = 1/h[(m2+1)*restart_steps+m2];
	  for(i=0;i<subdof;i++)
	    v[subdof*(m2+1)+i] = v[subdof*(m2+1)+i]*denom;
	  
	  //Givens rotation
	  for(i=0;i<m2;i++)  //finish Givens rotation from previous iterations
	    {
	      sum = h[i*restart_steps+m2];
	      h[i*restart_steps+m2] = rcos[i]*sum+rsin[i]*h[(i+1)*restart_steps+m2];
	      h[(i+1)*restart_steps+m2] = -rsin[i]*sum+rcos[i]*h[(i+1)*restart_steps+m2];
	    }
	  nu = sqrt(pow(h[m2*restart_steps+m2],2)+pow(h[(m2+1)*restart_steps+m2],2));
	  rcos[m2] = h[m2*restart_steps+m2]/nu; 
	  rsin[m2] = h[(m2+1)*restart_steps+m2]/nu;  
	  h[m2*restart_steps+m2] = rcos[m2]*h[m2*restart_steps+m2] + rsin[m2]*h[(m2+1)*restart_steps+m2]; 
	  h[(m2+1)*restart_steps+m2] = 0;
	  gg[m2+1] = -gg[m2]*rsin[m2];
	  gg[m2] = gg[m2]*rcos[m2];
	  	  
	  resid = fabs(gg[m2+1]);
	  m2++;
/*	  if(myid ==0)
	    printf("inner loop residual is %e zm0 is %e m2 is %d %e %e \n",
		   resid, zm0, m2,rcos[m2-1],rsin[m2-1]); */
	  
	}//done with inner loop of gmres
      double* y = new double[m2];
     
      for(i=m2-1;i>-1;i--)
	{
	  for(j=m2-1;j>i;j--)
	    gg[i] = gg[i] - h[i*restart_steps+j]*y[j];
	  y[i] = gg[i]/h[i*restart_steps+i];
	}
	    
      double* store_global = new double[totdof];

      if(myid != 0) // first proc stores solution from previous outer iteration
	for(i=0;i<totdof;i++) store_global[i] = 0;
      else
	for(i=0;i<totdof;i++) store_global[i] = u[i];

      for(i=0;i<m2;i++)
	{
	  for(j=0;j<nnss;j++)
	    if(c_loc[j] == 1)
	      store_global[ix[j]] += v[j+i*subdof]*y[i];
	  for(j=nnss;j<subdof;j++)
	    store_global[ix[j]] += v[j+i*subdof]*y[i];
	}
      i = MPI_Allreduce(store_global, u, totdof,MPI_DOUBLE,
			MPI_SUM,MPI_COMM_WORLD);
      delete []store_global;




      delete []y;
      m1++;
    }
  if(myid ==0)
    if(m1 == max_iter)
      printf("gmres did not converge!!!!!!!!!! \n");
  delete []store;
  delete []store2;
  delete []gg;
  delete []v;
  delete []h;
  double end = MPI_Wtime();  //timer
  end = end -start;

  
  if(myid == 0)
    {  
      printf("diagonal preconditioner\n");
      printf("iter = %d", m2+(m1-1)*restart_steps);
      printf("   resid = %e\n", resid);
      printf("corner and edge dof = %d\n", totdof);
      printf("   solution time = %e\n", end); 
    }
  MPI_Barrier(MPI_COMM_WORLD); 

  delete []diag_subdomain;
  delete []c_loc;
  delete []rcos;
  delete []rsin;

  delete []stiff;
  delete []load;
  delete []rpntr;
  delete []bpntr;
  delete []bindx;
  delete []indx;  
  //reconstruct bubble functions solution
  start = MPI_Wtime();
  Assemble_bubble(totdof, subdof, Bb, bb, u, ix, BT_Elem_Ptr, BT_Node_Ptr, myid);
  end = MPI_Wtime();  //timer
  end = end -start;
  printf("bubble reconstruction time = %e for proc %d\n", end, myid);
  return;
}
