#include "../header/hpfem.h"

extern void gross_process(double* u, int* ix, 
			  HashTable* HT_Node_Ptr, HashTable* HT_Elem_Ptr);

extern void reconstruct (HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr);



double* elem_sol(double* u, int* ix, HashTable* HT_Node_Ptr, HashTable* HT_Elem_Ptr,int myid, int numprocs,unsigned* key_el, double* Utemp)
{
  int i, j, k;
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

 
        gross_process(u, ix, HT_Node_Ptr, HT_Elem_Ptr);
        reconstruct(HT_Elem_Ptr, HT_Node_Ptr);

        int count =1;
        EmTemp = (Element*)HT_Elem_Ptr->lookup(key_el);
	if(!EmTemp->get_refined_flag())
	  {
	    ndof = EmTemp->get_no_of_dof();
	   
	
	    for (int kk = 0; kk < 242; kk ++) // initialize
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
		for( k = *(dofP+1)-*dofP+1; k < *(dofP+1)-*dofP+1+1*EQUATIONS; k+=EQUATIONS) {
		  *(Utemp+offset) = *(sol+k);
		  *(Utemp+1+offset) = *(sol+k+1);
		  //*(Utemp+2*offset) = *(sol+k);
		  offset += EQUATIONS;
		}
	      }
	      else {
		for(k=0;k<(*(dofP+1)-*dofP+1);k+=EQUATIONS) {
		  *(Utemp+offset) = *(sol+k);
		  *(Utemp+1+offset) = *(sol+k+1);
		
		  offset += EQUATIONS;
		}  
	      }
	      
	     
	    }
	     
	    //--bubble node
	    keyP = EmTemp->pass_key();
	    p = HT_Node_Ptr->lookup(keyP);
	    NdTemp = (Node*)p;
	    dofP = NdTemp->getdof(); 
	    sol = NdTemp->getsol();
	    
	    for(k=0;k<(*(dofP+1)-*dofP+1);k+=EQUATIONS)
	      {	
		Utemp[offset] = sol[k];
		Utemp[1+offset] = sol[k+1];
		//Utemp[2*offset] = sol[k];
		offset += EQUATIONS;
	      }
	    
	    count++;
	  }
	


  return Utemp;


}
  
