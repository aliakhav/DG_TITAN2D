////////////////////////////////////////
// Estimate the error of each element///
////////////////////////////////////////


#include "../header/hpfem.h"

extern "C" void elemerrp_(int* nequ, int* Norder, 
			  double* Xnod, double* Utemp, 
			  double* errorsq, double* solsq);


void err_estimate( HashTable* HT_Node_Ptr, HashTable* HT_Elem_Ptr )
{
  int i, j, k;
  HashEntryPtr   entryp;
  Element*       EmTemp;
  Node*          NdTemp;
  unsigned*      keyP;
  int*           dofP;
  int            ndof;
  int            offset;
  double*        sol;
  double*        Utemp;
  int e_buckets=HT_Elem_Ptr->get_no_of_buckets();


  for(i=0;i<e_buckets;i++){ 
    entryp = *(HT_Elem_Ptr->getbucketptr() + i); 
    while(entryp)
      { 
	EmTemp = (Element*)(entryp->value);
	if(!EmTemp->get_refined_flag())
	  {
	    ndof = EmTemp->get_no_of_dof();
	    Utemp = new double[2*ndof+64];
	    //-- +64 is used to avoid S_C_CON situation
	    offset = 0;
	    keyP = EmTemp->getNode();//--corner and edge nodes
	    for(j=0;j<8;j++)
	      {  
		p = HT_Node_Ptr->lookup(keyP+j*KEYLENGTH);
		NdTemp = (Node*)p;
		dofP = NdTemp->getdof();
		int nodeinfo = NdTemp->getinfo();	      

		if ( (nodeinfo == S_C_CON) && (j<4) ) 
		  {
		    for(k=*(dofP+1)-*dofP+1; k<*(dofP+1)-*dofP+1+1*EQUATIONS; k++) {
		      *(Utemp+2*offset) = *(sol+k);
		      offset++;
		    }
		  }

		else 
		  {
		    for(k=0;k<(*(dofP+1)-*dofP+1);k++){
		      *(Utemp+2*offset) = *(sol+k);
		      //--Utemp[2*(k+j)]=Utemp(1,k+j) in Fortran
		      //--k, kth equ of a node; j, jth node
		      //--U(2,..) is invalid for single equation probl
		      offset++;
		    }
		  }  
	      }
	    
	    keyP = EmTemp->pass_key();//--bubble node
	    p = HT_Node_Ptr->lookup(keyP);
	    NdTemp = (Node*)p;
	    dofP = NdTemp->getdof(); 

	    for(k=0;k<(*(dofP+1)-*dofP+1);k++)
	      {	
		Utemp[2*offset] = sol[k];
		offset++;
	      }

	    int      nequ = EQUATIONS;
	    double   errorsq[2] = {0.0, 0.0};
	    double   solsq[2] = {0.0, 0.0};
	    
	    int      Norder[5];
	    for(j=0;j<5;j++) Norder[j] = *(EmTemp->get_order() +j);

	    elemerrp_(&nequ, Norder, &Xnod[0][0], Utemp, errorsq, solsq);
	    
	    EmTemp->putel_sq(solsq, errorsq); 
	    
	    delete [] Utemp;	

	  }

	entryp = entryp->next;
	
      }
  }
  
  
}
