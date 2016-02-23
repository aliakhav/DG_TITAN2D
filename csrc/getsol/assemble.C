#include "../header/hpfem.h"

extern "C" void dgeco_(double*,int*, int*, int*, double*, double*);
void Assemble_stiffness(int tot_dof, double stiff[], double lod[], 
			HashTable* El_Table, HashTable* NodeTable, int myid, 
			int* rpntr, int* bpntr, int* bindx, int* indx, 
			int* num_blocks, int stiff_size, int NSVE, int* ix) 
/* first a mapping is created for each element to the subdomain
the contribution of each element is added to the sub_stiffness matrix*/ 

{
  int i,j,k, counter;
  double end, start, schur_time = 0;

  for(i=0; i<tot_dof; i++)
    lod[i]=0;
  
  for(i=0;i<stiff_size;i++)
    *(stiff+i) = 0.;

  int* block_loc = new int[tot_dof];
  int* sub_loc = new int[tot_dof];
  
  for(i=0;i<*num_blocks;i++)
    {
      k=0;
    for(j=*(rpntr+i);j<*(rpntr+1+i);j++)
      {
	*(block_loc+j) = i; //mapping to rpntr from row/column
	*(sub_loc+j) = k;  //mapping to sub matrix row/column
	k++;
      }
    }

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
               //print out coordinate of element centroid
	/*	Node* nd = (Node*) NodeTable->lookup(Curr_El->pass_key());
		printf(" i is : %d x is %e y is %e\n",i, *(nd->get_coord()), *(nd->get_coord()+1));*/
                
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
	/*	Node* NodePtr=(Node*) (NodeTable->lookup(Curr_El->pass_key()));
		for(int k=*(NodePtr->getdof()); k<*((NodePtr->getdof())+1)+1; k++)
		{
		ix_el[ix_counter] = k;
		ix_counter++;
		}*/  //bubble functions get condensed out
 //element to subdomain stiffness matrix mapping done

 //----assemble the sub_domain stiffness-------

		int count1=0;
		int count2=0;
		Curr_El->get_stiffness(NodeTable,El_Table, el_stiff, 
		       el_load, Curr_El);//--get the element stiffness and rhs	
		int int_dof=(Curr_El->get_order()) - 1; 
		int_dof = int_dof*int_dof*2; //# of bubble functions
		if(int_dof > 0)
		  {
		    start = MPI_Wtime();
		    //fact_matrix(el_stiff, el_load, ndof,int_dof);
		    end = MPI_Wtime();
		    schur_time += end - start;
		  }
		else
		  int_dof = 0;
		for(count1=0;count1<ndof-int_dof;count1++)
		  {
		    lod[ix_el[count1]] += el_load[count1];//ix_el[count1] is row index
		    
		    for(count2=count1;count2<ndof-int_dof;count2++)
		      {		    
			//vbr block storage
			int rblock = *(block_loc+ix_el[count1]);
			int r_remain = *(sub_loc+ix_el[count1]);
			int rblock2 = *(bpntr+rblock);
			int cblock = *(block_loc+ix_el[count2]);
			int c_remain = *(sub_loc+ix_el[count2]);
			counter = rblock2;
			while(cblock != *(bindx+counter))
			  counter++;  
			//counter is location in val of this block
			counter = *(indx+counter);
			counter += r_remain+c_remain*
			  (*(rpntr+rblock+1)-*(rpntr+rblock));

			//put in transpose of elm stiff to local stiff(fortran->C++)
			*(stiff+counter) +=*(el_stiff+ndof*count2+count1);
			
			if(count1 != count2)
			  {
			    rblock = *(block_loc+ix_el[count2]);
			    r_remain = *(sub_loc+ix_el[count2]);
			    rblock2 = *(bpntr+rblock);
			    cblock = *(block_loc+ix_el[count1]);
			    c_remain = *(sub_loc+ix_el[count1]);
			    counter = rblock2;
			    while(cblock != *(bindx+counter))
			      counter++;  
			    //counter is location in val of this block
			    counter = *(indx+counter);
			    counter += r_remain+c_remain*
			      (*(rpntr+rblock+1)-*(rpntr+rblock));

			    *(stiff+counter) +=*(el_stiff+ndof*count1+count2);
			  }
			
		      }
		    
		  }

		delete []el_stiff;
		delete []el_load;
		delete []ix_el;

	      }

	    currentPtr=currentPtr->next;      

	  }
	
      }

  printf("bubble elimination time = %e for proc %d\n", schur_time, myid);
  delete []block_loc;
  delete []sub_loc;
}
