//////////////////////////////////////////////////////////
/* associate the solution in solution vector u with node
   allocates solution memory in the node member for all 
   types of nodes except S_S_CON.  routine also inserts
   all non-constrained solution values (including
   side solutions of larger element at S_C_CON)         */
//////////////////////////////////////////////////////////
#include "../header/hpfem.h"


void gross_process(double* u, int* ix, 
		  HashTable* HT_Node_Ptr, HashTable* HT_Elem_Ptr)
{
  int i, j, k;
  HashEntryPtr   entryp;
  Node*          NdTemp;
  int*           dofP;
  double*        sol;
  int            nodeinfo;

  // scan node table
  
  int nodes = HT_Node_Ptr->get_no_of_buckets();

  for(i=0;i<nodes;i++){ 
    entryp = *(HT_Node_Ptr->getbucketptr() + i); 
    while(entryp)
      { 
	NdTemp = (Node*)(entryp->value);

	dofP = NdTemp->getdof();
	nodeinfo = NdTemp->getinfo();
	if ( (nodeinfo != S_S_CON) && (nodeinfo > 0) && (*dofP >= 0) ) 
	  {
	    if ( nodeinfo == S_C_CON ) 
	      {
		sol = new double[*(dofP+1)-*dofP+1 +EQUATIONS];
		for ( k = *(dofP+1)-*dofP+1; k < *(dofP+1)-*dofP+1 +EQUATIONS; k ++ )
		  sol[k] = 0.0; // initialization
		// EQUATIONS is used for reconstructed solution
	      }
	    else 
	      sol = new double[*(dofP+1)-*dofP+1];
	    
	    for(k=0;k<(*(dofP+1)-*dofP+1);k++)
	      *(sol+k) = *(u+*(ix+*dofP+k)); 
	    //--solution = *(u+*(ix+*dofP))

	    NdTemp->putsol(sol);
	    NdTemp->put_sol_deleted(0);
	  }
	else { }
	
	entryp = entryp->next;
	
      }
  }  
  
}
