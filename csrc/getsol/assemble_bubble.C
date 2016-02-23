#include "../header/hpfem.h"
#include "../header/blas.h"

void Assemble_bubble(int tot_dof, int subdof, int Bb, int bb, double u[], 
		     int* ix, 
		     HashTable* El_Table, HashTable* NodeTable, int myid) 
/* first a dof chain is created for each element of the subdomain
the contribution of each element is added to the sub_stiffness matrix*/ 

{
  int i,j,k;

  /*char  filename[14] = "ChanOutxx.dat";
  filename[7] = 48 + myid/10;  
  filename[8] = 48 + myid%10;
  ofstream eout;
  eout.open(filename, ios::out);  
  */
//  printf("in assemble_bubble %d\n", myid);  MPI_Barrier(MPI_COMM_WORLD);
  //-------------------go through all the elements of the subdomain------------------------
  //-------------------and create an element to subdomain mapping--------------------------

  HashEntryPtr* buck = El_Table->getbucketptr();
  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i))
      {
	HashEntryPtr currentPtr = *(buck+i);
	while(currentPtr)
	  {
	    Element* Curr_El=(Element*)(currentPtr->value);

 	    if(!(Curr_El->get_refined_flag()))//if this is a refined element don't involve!!!
	      {
		int      ndof = Curr_El->get_no_of_dof();
		double* el_stiff = new double[ndof*ndof];
		double* el_load  = new double[ndof];
		int* ix_el = new int[ndof];  //elm to subdof mapping
		int ix_counter = 0;

		for(j=0; j<8; j++)
		  {
		    Node* NodePtr=(Node*) (NodeTable->lookup(Curr_El->getNode()+KEYLENGTH*j));
		    if((NodePtr->getinfo()==CORNER)||(NodePtr->getinfo()==SIDE)||(NodePtr->getinfo()==BUBBLE))
		      {

			for(k=*(NodePtr->getdof()); k<*((NodePtr->getdof())+1)+1; k++)
			  {
			    ix_el[ix_counter] = k;
			    ix_counter++;				
			  }
		      }
		    else if(NodePtr->getinfo()==S_C_CON && j<4)
		      {

			
			Element* fatter=(Element*)El_Table->lookup(Curr_El->getfather());
			Node* TargetNode=(Node*)(NodeTable->lookup(fatter->getNode()+j*KEYLENGTH));//see notebook
			//looking up the appropriate corner node of the father
			//both of them have to have 1*EQUATIONS associated
			for(int k=*(TargetNode->getdof()); k<*((TargetNode->getdof())+1)+1; k++)
			  {
			    ix_el[ix_counter] = k;
			    ix_counter++;
			  }
			
		      }
		    else if(NodePtr->getinfo()==S_C_CON && j>3)//if it is an S_C_CON but it is a big element
		      {
			for(int k=*(NodePtr->getdof()); k<*((NodePtr->getdof())+1)+1; k++)
			  {
			    ix_el[ix_counter] = k;
			    ix_counter++;
			  }
		      }


		    else if(NodePtr->getinfo()==S_S_CON)//the side nodes of the constrained element is added to the constarined node
		      {
			Node* TargetNode=(Node*)(NodeTable->lookup(Curr_El->getNode()+(j-4)*KEYLENGTH));//see note
			if(TargetNode->getinfo()!=S_C_CON)
			  
			  if(j!=7) TargetNode=(Node*)(NodeTable->lookup(Curr_El->getNode()+(j-3)*KEYLENGTH));
			  else TargetNode=(Node*)(NodeTable->lookup(Curr_El->getNode()));
			  			    
			for(int k=*(TargetNode->getdof()); k<*((TargetNode->getdof())+1)+1; k++)
			  {
			    ix_el[ix_counter] = k;
			    ix_counter++;
 			  }
		    		  		    
		      }
		  }
		Node* NodePtr=(Node*) (NodeTable->lookup(Curr_El->pass_key()));
		for(int k=*(NodePtr->getdof()); k<*((NodePtr->getdof())+1)+1; k++)
		  {
		    ix_el[ix_counter] = k;
		    ix_counter++;
		  }
 //element to subdomain stiffness matrix mapping done
 //----assemble the sub_domain stiffness-------

		int count1=0;
		int count2=0;
		//--get the element stiffness and rhs
		Curr_El->get_stiffness(NodeTable, El_Table, 
				       el_stiff, el_load, Curr_El);
		int int_dof = (Curr_El->get_order()) - 1; 
		int_dof = int_dof*int_dof*2; //# of bubble functions
		if(int_dof < 0)
		  int_dof = 0;
		for(k=0;k<ndof-int_dof;k++)
		  for(j=ndof-int_dof;j<ndof;j++)
		    *(el_load+j) -= *(el_stiff+k*ndof+j) * *(u+ix[ix_el[k]]);
		    
		double* int_el_stiff = new double[int_dof*int_dof];
		for(j=0;j<int_dof;j++)
		  for(k=0;k<int_dof;k++)
		    *(int_el_stiff+j*int_dof+k) = 
		      *(el_stiff+(j+ndof-int_dof)*ndof+k+ndof-int_dof);

		if(int_dof > 0)
		  {
#ifdef SUNOS		
		    tri_(int_el_stiff,&int_dof, &int_dof);
#endif
		    
#ifdef IBMSP
		    tri(int_el_stiff,&int_dof, &int_dof);
#endif
		    
#ifdef CRAY
		    TRI(int_el_stiff,&int_dof, &int_dof);
#endif
		  }    
		delete []el_stiff;
		double* u_bubble = new double[int_dof];
		if(int_dof > 0)
		  {
#ifdef SUNOS
		    //printf("rhsub ndof %d int_dof %d i %d proc %d \n",ndof, int_dof, i, myid);
		    //cout<<*(currentPtr->key)<<" "<<*(currentPtr->key+1)<<endl<<flush;
		    rhsub_(int_el_stiff, (el_load+ndof-int_dof),&int_dof, &int_dof, u_bubble); 
		    //		printf("done rhsub %d \n",myid);
#endif
		    
#ifdef IBMSP
		    rhsub(int_el_stiff, (el_load+ndof-int_dof),&int_dof, &int_dof, u_bubble);  
#endif
		    
#ifdef CRAY
		    RHSUB(int_el_stiff, (el_load+ndof-int_dof),&int_dof, &int_dof, u_bubble);  
#endif
		  }
		for(j=0;j<int_dof;j++)  //insert u_bubble into global solution
		  { 
		    u[ix[ix_el[ndof-int_dof+j]]] = *(u_bubble+j);
		  }
		  
		    
		delete []u_bubble;
		delete []el_load;
		delete []ix_el;
	      }
	    currentPtr=currentPtr->next;      
	  }
      }
//  printf("done with assemble_bubble\n");
}
