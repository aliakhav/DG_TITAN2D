//---------------------------------------------------------------------
//   routine name       - solelm
//   latest revision    - September, 1991
//   purpose            - routine gives solution at all degree freedoms
//                        of an element
//   arguments :
//     in:    Nel,Kind  - element number and kind
//     out:   Sol       - local dof of the element: Nel, kind 
//---------------------------------------------------------------------

#include "../header/hpfem.h"

extern void eval1 ( int, double, double[2][11], double[2] );

#ifdef SUNOS
extern "C" void eval2_(int*, double*, int*, double*);
#endif

#ifdef IBMSP
extern "C" void eval2(int*, double*, int*, double*);
#endif

#ifdef CRAY
extern "C" void EVAL2(int*, double*, int*, double*);
#endif


void reconstruct (HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr) 
{

  // variables used in eval1  
  double     Sol[EQUATIONS][121];
  double     Nvalues[EQUATIONS][11], val1[EQUATIONS];

  // variables used in table scanning
  Element*   EmTemp;
  Node*      NdTemp;
  unsigned*  nodekey;
  HashEntryPtr  entryp;
  int        i, j, k, l;
  int        Norder;

  

  // initialize variables
  for ( i=0; i<EQUATIONS; i++)
    {
      val1[i] = 0.0;
      for ( j=0; j<11; j++ )
	Nvalues[i][j] = 0.0;
    }

  // find myid
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);


  // scan element table, do eval1, eval2
  int elements = HT_Elem_Ptr->get_no_of_buckets();
  for(int iii=0;iii<elements;iii++){ 
    entryp = *(HT_Elem_Ptr->getbucketptr() + iii); 
    while(entryp)
      { 
	EmTemp = (Element*)(entryp->value);
	if(!EmTemp->get_refined_flag())
	  {
	    nodekey = EmTemp->getNode();
	    unsigned* neigh_key = EmTemp->get_neighbors();
	    int  gen  =  EmTemp->get_gen();
	    int* neigh_proc = EmTemp->getassoc();

	    for ( i=0; i<8; i++ ) 
	      {		
		NdTemp = (Node*) ( HT_Node_Ptr->lookup( nodekey + i*KEYLENGTH ) );
		if ( i < 4 ) // vertices
		  {
		    if ( (NdTemp->getinfo() == S_C_CON) && (!NdTemp->get_reconstructed()) )
		      {
			// find out nod1 and nod2
			// first, get father
			Element* father = (Element*) ( HT_Elem_Ptr->lookup( EmTemp->getfather() ) );
			assert ( father );  // a watching window for repartition
			
			// secondly, find side
			int  neigh_gen;
			if ( *(neigh_proc + i) == myid )
			  {
			    Element* neighbor = (Element*) ( HT_Elem_Ptr->lookup( neigh_key + i*KEYLENGTH ) ); 
			    neigh_gen = neighbor->get_gen(); 
			  }			  
			else
			  neigh_gen = *( EmTemp->get_neigh_gen() + i );

			int  side;
			if ( neigh_gen < gen ) 
			  { side = i;}
			else 
			  {
			    if ( i-1<0 ) side = 3;
			    else side = i - 1;
			  }
			
			// nod1 and nod2
			unsigned* fathernod = father->getNode();
			// nod1
			Node* nod1 = (Node*)( HT_Node_Ptr->lookup( fathernod + side*KEYLENGTH ) );
			assert( nod1->getinfo() == CORNER );
			
			// nod2
			if ( side == 3 ) side = 0;
			else side = side + 1;
			Node* nod2 = (Node*)( HT_Node_Ptr->lookup( fathernod + side*KEYLENGTH ) );
			assert( nod2->getinfo() == CORNER );
			
			// nod3 nod3 = NdTemp
			
			// fill out Nvalues (Dof)
			for(k=0;k<EQUATIONS;k++)
			  {
			    Nvalues[k][0] = *( nod1->getsol()+k);
			    Nvalues[k][1] = *( nod2->getsol()+k);
			  }

			Norder = NdTemp->get_order();
			int* dof = NdTemp->getdof();
			for ( k=0; k<*(dof+1)-*dof+1; k += EQUATIONS ) 
			  {	
			    for(l=0;l<EQUATIONS;l++)
			      Nvalues[l][k/EQUATIONS+2] = *(NdTemp->getsol()+k+l);   
			  }
			
			// eval1, get the physical solution of constrained nodes
			double Eta = 0.0;
			eval1( Norder, Eta, Nvalues, val1);
			
			// val1[EQUATION]: physical solution
			// put physical solution into NdTemp
			// an S_C_CON has extra places in its sol for storing
			// the interplated value e.g.
			// *dof ~ *(dof+1)-*dof is used for higher order value
			// *(dof+1)-*dof+1 ~ *(dof+1)-*dof+1+1*EQUATIONS-1 is used for the linear (physical) value
			double* sol = NdTemp->getsol();
			for ( k=0; k<EQUATIONS; k++)
			  *(sol+k+ *(dof+1)- *dof +1) = val1[k];
			
			// mark nod3, avoid doing the same construction
			NdTemp->put_reconstructed ( 1 );
		      }		  
		    
		  }
		else  // mid-nodes
		  {
		    if (NdTemp->getinfo() == S_S_CON)
		      {
			// reconize iposi
			int jk = i - 3;
			if ( jk == 4 ) jk = 0;
			int  iposi  = 0;
			Node* nod3 = (Node*)(HT_Node_Ptr->lookup( nodekey + jk*KEYLENGTH ) );

			if ( nod3->getinfo() == S_C_CON)  
			  {
			    if(i < 6)
			      iposi = 1;
			    else
			      iposi = 2;
			  }
			else  
			  {
			    if(i < 6)
			      iposi = 2;
			    else
			      iposi = 1;
			    nod3 = (Node*)(HT_Node_Ptr->lookup( nodekey + (i-4)*KEYLENGTH ) );
			    assert ( nod3->getinfo() == S_C_CON );
			  }

			Norder = nod3->get_order();
			int*  dof = nod3->getdof();
			double*  sol = nod3->getsol();

			// fill in Nvalues2 (Dof)
			double* Nvalues2 = new double[11*EQUATIONS];
			for ( k=0; k<*(dof+1)-*dof+1; k++)
			  Nvalues2[k] = sol[k];	
			
			double* val_s_s_con = new double[11*EQUATIONS];
			// call eval2_
			// val_s_s_con brings back the values 
			// (polynomial coeff associated with S_S_CON
#ifdef SUNOS		       
			eval2_(&Norder, Nvalues2, &iposi, val_s_s_con);
#endif

#ifdef IBMSP
			eval2(&Norder, Nvalues2, &iposi, val_s_s_con);
#endif

#ifdef CRAY
			EVAL2(&Norder, Nvalues2, &iposi, val_s_s_con);
#endif
			
			delete []Nvalues2;
			// put val_s_s_con to S_S_CON sol
			int nodeorder = *( EmTemp->get_order() + i%4 );
			double* sol_SS = new double[(nodeorder-1)*EQUATIONS];
			for ( k=0; k<*(dof+1)-*dof+1; k+=EQUATIONS ) {
			  for(l=0;l<EQUATIONS;l++)
			    *(sol_SS+k+l) = val_s_s_con[k+l];
			}
			NdTemp->putsol( sol_SS );
			NdTemp->put_sol_deleted(0);

		      }
		  }
	      
	      }

	  }

	else {}

	entryp = entryp->next;
      }

  }

}
		    
