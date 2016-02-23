#include "../header/hpfem.h"

void insert_blockdof(dofPtr*, int, int);
void deletingdof2(dofPtr);
int make_block(int sub_dof, int** rpntr, int** bpntr, int** bindx, int** indx, 
	       int* num_blocks, int* ix, HashTable* El_Table, HashTable* NodeTable,
	       int nn, int ss, int vv, int ee)
  /*routine to determine the amount of blocks and their size
  block structure will be symmetric
  make_block.gi.C requires less computation but for large problems make_block.gi.C 
  will require more memory than the storage needed for the subdomain stiffness matrix

  make_block.C determines the sparsity structure of the stiffness matrix and fills in the 
  integer arrays rpntr, bpntr, bindx, indx and size of stiffness matrix

  blocks are stored by columns

  variable block sizes are determined from the dof at each corner and edge node
  each block row contains a linked list holding the non-zero pattern for that row

  for a detailed description of the VBR sparse storage format, go to:
  http://wings.buffalo.edu/eng/mae/acm2e/afeapi_download.html
*/
{
		int myid;
		MPI_Comm_rank(MPI_COMM_WORLD, &myid); 
  int i,j,k;
  int* block_size = new int[sub_dof];
  int* block_loc = new int[sub_dof];
  for(i=0;i<sub_dof;i++)
    {
      *(block_size+i) = 0;
      *(block_loc+i) = 0;
    }
  HashEntryPtr* buck = El_Table->getbucketptr();
  int buckets = El_Table->get_no_of_buckets();
  for(i=0; i<buckets; i++)
    if(*(buck+i))
      {
	HashEntryPtr currentPtr = *(buck+i);
	while(currentPtr)
	  {  
	    Element* Curr_El=(Element*)(currentPtr->value);
	    if(!(Curr_El->get_refined_flag()))//if this is a refined element don't involve!!!
	      {
		for(j=0; j<8; j++)
		  {
		    Node* NodePtr=(Node*) (NodeTable->lookup(Curr_El->getNode()+KEYLENGTH*j));
		    if((NodePtr->getinfo()==CORNER)||(NodePtr->getinfo()==SIDE)
		       ||NodePtr->getinfo()==S_C_CON) 
		      {
			int start = *(NodePtr->getdof());
			int end = *((NodePtr->getdof())+1);
			int dof = end+1-start;  //amount of unknowns for this node(block size)
			if(dof > 0)
			  {
			    *(block_size+start) = dof;
			    for(k=start;k<end+1;k++)
			      *(block_loc+k) = start;  //points to beginning of block
			  }
		      }
		  }
		int* jjjj = Curr_El->get_order();
		assert(jjjj);
	      }

	    currentPtr=currentPtr->next;   
	  }
      }
  *num_blocks = 0;
  for(i=0;i<sub_dof;i++)
    if(*(block_size+i))       
      (*num_blocks)++;
  
  *rpntr= new int[*num_blocks+1];  //note that cpntr is same as rpntr
  *(*rpntr) = 0;
  int counter=0;
  for(i=0;i<sub_dof;i++)
    {
      if(*(block_size+i))
	{
	  counter++;
	  *(*rpntr+counter) = *(*rpntr+counter-1)+*(block_size+i);
	  for(j=i;j<(i+*(block_size+i));j++)
	    *(block_loc+j) = counter-1;  //mapping to rpntr from row/column
	}
    }

  dofPtr* row_block_dofPtr = new dofPtr[*num_blocks];
  for(i=0;i<(*num_blocks);i++)  row_block_dofPtr[i] = NULL;
  
  //-------------------go through all the elements of the subdomain------------------------
  //-------------------and create an element to subdomain mapping--------------------------
  buck = El_Table->getbucketptr();
  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i))
      {
	HashEntryPtr currentPtr = *(buck+i);
	int icounter = 0;
	while(currentPtr)
	  {
	    icounter++;
	    Element* Curr_El=(Element*)(currentPtr->value);

 	    if(!(Curr_El->get_refined_flag()))//if this is a refined element don't involve!!!
	      {
		int      ndof = Curr_El->get_no_of_dof();
		int* ix_el = new int[ndof];  //elm to subdof mapping
		int ix_counter = 0;

		for(j=0; j<8; j++)
		  {
		    Node* NodePtr=(Node*) (NodeTable->lookup(Curr_El->getNode()+KEYLENGTH*j));
		    int abcd = S_S_CON;
		    if(icounter ==15)
		      abcd = NodePtr->getinfo();
		    if((NodePtr->getinfo()==CORNER)||(NodePtr->getinfo()==SIDE)||
		       (NodePtr->getinfo()==BUBBLE))
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
			  		
			int kstart = *(TargetNode->getdof());
			int kend = *((TargetNode->getdof())+1)+1;
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
		int int_dof=*(Curr_El->get_order()+4) - 1; 
		int_dof = int_dof*int_dof*2; //# of bubble functions
		if(int_dof < 0) 
		  int_dof = 0;
		for(count1=0;count1<ndof-int_dof;count1++)		
		  for(count2=0;count2<ndof-int_dof;count2++)
		    insert_blockdof(row_block_dofPtr,*(block_loc+ix_el[count1]),
				    *(block_loc+ix_el[count2]));

		delete []ix_el;				
	      }
	    currentPtr=currentPtr->next;      

	  }
	
      }
  //-----linked chain finished------
/*  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  for(i=0;i<*num_blocks;i++)
    {
      printf("\n chain: %d\n",i);      //print out the chain to file.
      dofPtr tempdof=row_block_dofPtr[i];
      while(tempdof)
	{
	  cout<<tempdof->value<<"\n"<<flush;
	  tempdof=tempdof->next;
	}
      printf("linked chain finished: %d\n", myid);
    }
  MPI_Barrier(MPI_COMM_WORLD); */    
  *bpntr = new int[*num_blocks+1];
  *(*bpntr) = 0;
  counter = 0;
  for(i=0;i<*num_blocks;i++)
    {
      dofPtr currentPtr =row_block_dofPtr[i]; 
      while(currentPtr)
	{
	  counter++;
	  currentPtr = currentPtr->next;
	}

      *(*bpntr+i+1) = counter;
    }
  //printf("rpntr and bpntr are size: %d\n",*num_blocks+1);
  //printf("bindx and indx are size: %d\n",*(*bpntr+*num_blocks)+1);
  *bindx = new int[*(*bpntr+*num_blocks)+1];
  counter=0;
  for(i=0;i<*num_blocks;i++)
    {

      dofPtr currentPtr =row_block_dofPtr[i]; 
      while(currentPtr)
	{
	  *(*bindx+counter) = currentPtr->value;
	  counter++;
	  currentPtr = currentPtr->next;
	}

    }
  
  counter = 0; //find out size of val & indx arrays
  *indx = new int[*(*bpntr+*num_blocks)+1];
  *(*indx) = 0;
  int counter2 = 0;
  for(i=0;i<*num_blocks;i++)
    {
      dofPtr currentPtr =row_block_dofPtr[i]; 
      while(currentPtr)
	{
	  j=currentPtr->value;
	  counter2 += (*(*rpntr+i+1)-*(*rpntr+i))*(*(*rpntr+j+1)-*(*rpntr+j));
	  *(*indx+counter+1) = counter2;
	  counter++;
	  currentPtr = currentPtr->next;
	}
    }
  for(i=0;i<(*num_blocks);i++)
    deletingdof2(row_block_dofPtr[i]);
  delete []row_block_dofPtr;
  delete []block_size;
  delete []block_loc;
  //printf("stiffness is of size: %d\n",val_size);
  return counter2;
}


//insert new object into linked list
void insert_blockdof(dofPtr* row_block_dofPtr, int blockrow, int blockcol)
{

  dofPtr new_dof, prev_dof;
  
  if(row_block_dofPtr[blockrow] == NULL)
    row_block_dofPtr[blockrow]= new dof(blockcol);
      
  else
    {
      dofPtr curr_dof;
      curr_dof = row_block_dofPtr[blockrow];
      if(curr_dof->value > blockcol)
	{
	  new_dof = new dof(blockcol);
	  new_dof->next = curr_dof;
	  row_block_dofPtr[blockrow] = new_dof;
	}
      else
	{
	  while(curr_dof->value < blockcol)
	    {
	      prev_dof = curr_dof;
	      curr_dof = curr_dof->next;
	      if(curr_dof == NULL)  //new dof is last in linked list
		{
		  new_dof = new dof(blockcol);
		  prev_dof->next = new_dof;
		  return;
		}
	    }
	  if(curr_dof->value != blockcol )
	    {
	      new_dof = new dof(blockcol);
	      new_dof->next = curr_dof;
	      prev_dof->next = new_dof;
	    }
	}
    }
  return;
}
//delete linked list 
void deletingdof2(dofPtr headPtr)
{
  dofPtr currentPtr=headPtr;
  dofPtr tempPtr;

  while(currentPtr)
    {
      tempPtr=currentPtr;
      currentPtr=currentPtr->next;
      delete tempPtr;
    }


}
