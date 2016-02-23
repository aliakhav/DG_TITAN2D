#include "../header/hpfem.h"


#ifdef SUNOS
extern "C" void herror_(int* nequ, int* Norder, double* Xnod,
			double* Utemp, double* errorsq, double* solsq, 
			int* matid );
extern "C" void elemerrp_(int* nequ, int* Norder, double* Xnod,
			double* Utemp, double* errorsq, double* solsq);

extern void pre_flex(HashTable* ht_elem_ptr,HashTable* ht_node_ptr,double* u,int* ix);
#endif

#ifdef IBMSP
extern "C" void herror(int* nequ, int* Norder, double* Xnod,
		       double* Utemp, double* errorsq, double* solsq, 
		       int* matid );
#endif

#ifdef CRAY
extern "C" void HERROR(int* nequ, int* Norder, double* Xnod,
		       double* Utemp, double* errorsq, double* solsq, 
		       int* matid );
#endif


extern void gross_process(double* u, int* ix, 
			  HashTable* HT_Node_Ptr, HashTable* HT_Elem_Ptr);

extern void reconstruct (HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr);
void Post_process1(double* u, int* ix, 
		  HashTable* HT_Node_Ptr, HashTable* HT_Elem_Ptr,
		  NNLink* NNHead, VVLink* VVHead, SSLink* SSHead,
		  EELink* EEHead, BBLink* BBHead,
		  int myid, int numprocs, int h_count)
{
  int i, j, k,kk;
  HashEntryPtr   entryp;
  Element*       EmTemp;
  Node*          NdTemp;
  unsigned       KeyTemp[KEYLENGTH];
  unsigned*      keyP;
  void*          p;
  int*           assocP;
  int*           dofP;
  int            ndof;
  int            power;
  int            interf;
  int            offset;
  double*        sol;
  double*        coord;
  double*        Utemp;
  double         Xnod[9][2]; 

  double	 debug[242];

  char           filename[17] = "solutionxxxx.dat";
  FILE*          fp;
  int		 material;

  filename[8] = 48+POWER%10;
  filename[9] = 48+myid/10;
  filename[10] = 48+myid%10;
  filename[11] = 48+h_count;

  int e_buckets=HT_Elem_Ptr->get_no_of_buckets();

  ///////////////////////////////////////////////
  gross_process(u, ix, HT_Node_Ptr, HT_Elem_Ptr);
  ///////////////////////////////////////////////
  reconstruct(HT_Elem_Ptr, HT_Node_Ptr);
  ///////////////////////////////////////////////

  int count = 1;

  for(i=0;i<e_buckets;i++){ 
    entryp = *(HT_Elem_Ptr->getbucketptr() + i); 
    while(entryp)
      { 
	EmTemp = (Element*)(entryp->value);
	if(!EmTemp->get_refined_flag())
	  {
	    ndof = EmTemp->get_no_of_dof();
	    //Utemp = new double[2*ndof+64];	    
	    //-- +64 is used to avoid S_C_CON situation
	    Utemp = new double[242];
	    for (kk = 0; kk < 242; kk ++) // initialize
	      Utemp[kk] = 0.0;
	    offset = 0;

	    keyP = EmTemp->getNode();
	    //--corner and edge nodes
	    for(j=0;j<8;j++){  
	      p = HT_Node_Ptr->lookup(keyP+j*KEYLENGTH);
	      assert(p); //--watching window
	      NdTemp = (Node*)p;
	      dofP = NdTemp->getdof();
	      int nodeinfo = NdTemp->getinfo();
	      sol = NdTemp->getsol();

	      if ( (j < 4) && ( nodeinfo == S_C_CON) ) {
		for( k = *(dofP+1)-*dofP+1; k < *(dofP+1)-*dofP+1+EQUATIONS; k+=EQUATIONS) {
		  for(kk=0;kk<EQUATIONS;kk++)
		    *(Utemp+kk+offset) = *(sol+k+kk);
		  offset += EQUATIONS;
		}
	      }
	      else {
		for(k=0;k<(*(dofP+1)-*dofP+1);k+=EQUATIONS) {
		  for(kk=0;kk<EQUATIONS;kk++)
		    *(Utemp+kk+offset) = *(sol+k+kk);
		  offset += EQUATIONS;
		}  
	      }
	      coord = NdTemp->get_coord();
	      Xnod[j][0] = *coord;
	      Xnod[j][1] = *(coord+1);
	    }
	     
	    //--bubble node
	    keyP = EmTemp->pass_key();
	    p = HT_Node_Ptr->lookup(keyP);
	    NdTemp = (Node*)p;
	    dofP = NdTemp->getdof(); 
	    sol = NdTemp->getsol();
	    
	    for(k=0;k<(*(dofP+1)-*dofP+1);k+=EQUATIONS)
	      {	
		for(kk=0;kk<EQUATIONS;kk++)
		  *(Utemp+kk+offset) = *(sol+k+kk);
		offset += EQUATIONS;
	      }
	    
	    coord = NdTemp->get_coord();
	    Xnod[j][0] = *coord;
	    Xnod[j][1] = *(coord+1);
	    
	    int      nequ = EQUATIONS;
	    double   errorsq=0, solsq=1; // returned values from estimator
	    
	    int      Norder[5];
	    for(j=0;j<5;j++) Norder[j] = *(EmTemp->get_order() +j);

	    int Nc = EmTemp->get_no_of_dof();
	    material=EmTemp->get_material();


	    for(int o=0; o<242; o++)
		debug[o]=*(Utemp+o);

#ifdef SUNOS
	    herror_(&nequ, Norder, &Xnod[0][0], Utemp, &errorsq, &solsq, 
		    &material); 

#endif

#ifdef IBMSP
	    herror(&nequ, Norder, &Xnod[0][0], Utemp, &errorsq, &solsq, 
		   &material); 
#endif

#ifdef CRAY
	    HERROR(&nequ, Norder, &Xnod[0][0], Utemp, &errorsq, &solsq, 
		   &material); 
#endif

	    
	    EmTemp->putel_sq(solsq, errorsq); 
	    
	    delete []Utemp;	
	    count++;
	  }
	
	entryp = entryp->next;
	
      }
    
  }

//  pre_flex(HT_Elem_Ptr, HT_Node_Ptr, u, ix); 
}
  
