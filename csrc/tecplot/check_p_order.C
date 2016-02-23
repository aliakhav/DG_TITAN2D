#include "../header/hpfem.h"

void check_p_order(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr)
{


  int myid;
  int numprocs;
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
  char filename[13] = "porplott.plt";
  //filename[7] = which + 48;
  int order;


  //ofstream outData;
  FILE*    fp;
 
  if ( myid == 0 )
    {
      //outData.open(filename, ios::out); 
      fp = fopen ( filename, "w" );
      //outData<<"TITLE= \"MESH OUTPUT\""<<'\n';
      //outData<<"VARIABLES = \"X\", \"Y\", \"SOL\", \"P\"";
      fprintf ( fp, "TITLE= \"MESH OUTPUT\"\n" );
      fprintf ( fp, "VARIABLES = \"X\", \"Y\", \"SOL\", \"P\"" );
    }
  else 
    {
      MPI_Recv(&done, 1, MPI_INT, myid-1, TECTAG, MPI_COMM_WORLD, &status);
      //outData.open(filename, ios::app);
      fp = fopen ( filename, "a+" );
    }


  for(int i=0; i<500; i++)
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


  //outData<<'\n';
  fprintf ( fp, "\n" );

  //outData<<"ZONE N="<<element_counter*4<<", E="<<element_counter<<", F=FEPOINT, ET=QUADRILATERAL"<<'\n';
  fprintf ( fp, "ZONE N=%d, E=%d, F=FEPOINT, ET=QUADRILATERAL\n", element_counter*4, element_counter );

  for(i=0; i<500; i++)
    {

      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
	{	
	  
	  EmTemp = (Element*)entryp->value;
	  assert(EmTemp);
	  if(!(EmTemp->get_refined_flag()))
	    {
	      
	      nodes = EmTemp->pass_nodes();
	      for(int j=0; j<4; j++)
		{
		  order=*(EmTemp->get_order()+j);
		  NodeTemp = (Node*) HT_Node_Ptr->lookup(nodes+j*KEYLENGTH);
		  int* dof = NodeTemp->getdof();		  
		  
		  fprintf ( fp, "%f %f %d %d\n", *(NodeTemp->get_coord()), *(NodeTemp->get_coord()+1), 0, order);

		}  
	    } 
	  
	  entryp = entryp->next; 
	  
	} 
    } 

  //outData<<'\n'; 
  fprintf ( fp, "\n" );

  for(i=0; i<element_counter;i++) 
    { 
      for(int j=0; j<4; j++) 
	//outData<<i*4+j+1<<' '; 
	fprintf (fp, "%d ", i*4+j+1);

      //outData<<'\n'; 
      fprintf ( fp, "\n" );
    }
 
  fclose ( fp );

  if(myid != numprocs-1) MPI_Send(&done, 1, MPI_INT, myid+1, TECTAG, MPI_COMM_WORLD); 

} 

