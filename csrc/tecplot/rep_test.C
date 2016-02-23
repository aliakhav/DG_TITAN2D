#include "../header/hpfem.h"

void rep_test(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, int which)
{


  int myid;
  int numprocs;
  int material;
  int done = 1; 
  int TECTAG = 123;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  int element_counter = 0;
  Element* EmTemp;
  Node* NodeTemp;
  HashEntry* entryp;
  unsigned* nodes;
  char filename[13] = "reptestx.plt";
  filename[7] = which + 48;
  int order;
  int e_buckets=HT_Elem_Ptr->get_no_of_buckets();

  unsigned key;


  FILE*    fp;
  
  ofstream fout;
 
  if ( myid == 0 )
    {
      fout.open(filename, ios::out);
      fout<<"TITLE= \"MESH OUTPUT\"\n";
      fout<<"VARIABLES = \"X\", \"Y\", \"P\", \"KEY\""<<endl;
    }
  else 
    {
      MPI_Recv(&done, 1, MPI_INT, myid-1, TECTAG, MPI_COMM_WORLD, &status);
      fout.open(filename, ios::app);
    }


  for(int i=0; i<e_buckets; i++)
    {

      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
	{	
	  
	  EmTemp = (Element*)entryp->value;
	  assert(EmTemp);
	  if(!(EmTemp->get_refined_flag()))
	    element_counter++;
	  entryp = entryp->next;
	  
	}
    }


  fout<<'\n';

  fout<<"ZONE N="<<element_counter*4<<", E="<<element_counter<<", F=FEPOINT, ET=QUADRILATERAL\n";

  int elements = HT_Elem_Ptr->get_no_of_buckets();

  for(i=0; i<elements; i++)
    {

      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
	{	
	  
	  EmTemp = (Element*)entryp->value;
	  assert(EmTemp);
	  if(!(EmTemp->get_refined_flag()))
	    {
	      order = 2;
	      for ( int k = 0; k < 5; k++ )
		{
		  int help_order=*(EmTemp->get_order()+k);
		  if ( help_order > order ) order = help_order;
		}

	      nodes = EmTemp->getNode();
	      material=EmTemp->get_material();
	      key=*(EmTemp->pass_key());
	      for(int j=0; j<4; j++)
		{
		  NodeTemp = (Node*) HT_Node_Ptr->lookup(nodes+j*KEYLENGTH);
		  int* dof = NodeTemp->getdof();
		  if ( NodeTemp->getinfo() != S_C_CON )
		    fout<<*(NodeTemp->get_coord())<<" "<<*(NodeTemp->get_coord()+1)<<" "<<order<<" "<<key<<endl;

		  else 
		    fout<<*(NodeTemp->get_coord())<<" "<<*(NodeTemp->get_coord()+1)<<" "<<order<<" "<<key<<endl;
		}  
	    } 
	  
	  entryp = entryp->next; 
	  
	} 
    } 

  fout<<endl;

  for(i=0; i<element_counter;i++) 
    { 
      for(int j=0; j<4; j++) 
	fout<<i*4+j+1<<' '; 

      fout<<'\n'; 
    }
 

  if(myid != numprocs-1) MPI_Send(&done, 1, MPI_INT, myid+1, TECTAG, MPI_COMM_WORLD); 

} 

