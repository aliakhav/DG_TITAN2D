#include "../header/hpfem.h"


/*************************************
**************************************
**************************************
*  output mesh only
**************************************
**************************************
*************************************/
void meshplotter(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, int which, 
		 MatProps* matprops)
{
  int myid, i, k;
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
  char filename[20] = "mshplotxxxxxxxx.plt";
  filename[14] = (which % 10) + 48;
  filename[13] = (which % 100)/10 + 48;
  filename[12] = (which % 1000)/100 + 48;
  filename[11] = (which % 10000)/1000 + 48;
  filename[10] = (which % 100000)/10000 + 48;
  filename[9] = (which % 1000000)/100000 + 48;
  filename[8] = (which % 10000000)/1000000 + 48;
  filename[7] = (which % 100000000)/10000000 + 48;
  filename[6] = (myid%10) +48;
  filename[5] = (myid%100)/10 +48;
  int order;
  int e_buckets=HT_Elem_Ptr->get_no_of_buckets();

  double momentum_scale = matprops->HEIGHT_SCALE * sqrt(matprops->LENGTH_SCALE * (matprops->GRAVITY_SCALE)); // scaling factor for the momentums

  FILE*    fp;
 
    {
      fp = fopen ( filename, "w" );
      fprintf ( fp, "TITLE= \"MESH OUTPUT\"\n" );
      fprintf ( fp, "VARIABLES = \"X\" \"Y\" \"NODAL_ELEVATION\" \"PROC\" \"PILE_HEIGHT\" \"XMOMENTUM\" \"YMOMENTUM\" \"KEY0\" \"KEY1\" \"GENERATION\" \"SON\" \"ORDER\" \"ELMLOC1\" \"ELMLOC2\" \"Flux\" \"Error1\" \"Error2\" \"Error3\" " );
    }


  for(i=0; i<e_buckets; i++)
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

  int elements = HT_Elem_Ptr->get_no_of_buckets();

  for(i=0; i<elements; i++)
    {//begin loop elements

      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
	{	
	  
	  EmTemp = (Element*)entryp->value;
	  assert(EmTemp);
	  if(!(EmTemp->get_refined_flag()))
	    {
	      /*if(*(EmTemp->get_elm_loc()) == 16 && *(EmTemp->get_elm_loc()+1) == 8)
		k = 0;*/
	      nodes = EmTemp->getNode();
	      	      for(int j=0; j<4; j++)
	      		{//begin corner loop 
	     		  double sol[3] = {0,0,0};
		  int ndofe = EmTemp->get_ndof()/EQUATIONS;
		  double* elsol = EmTemp->get_el_solution();
		  int l, order = EmTemp->get_order();
		  double* f = new double[ndofe];
		  double x[2];
		  		  if(j==0 || j==3)
		    x[0] = -1;
		  else
		    x[0] = 1;
		  if(j==3 || j == 2)
		    x[1] = 1;
		  else 
		  x[1] = -1;
		  //		  x[0] = 0;x[1]=0;

		  shape2dg_(&order, x, (x+1), f);
		  if(elsol[0] > (double) .00001)
		    l = 1;
		  for(l=0;l<ndofe;l++) 
		    for(k=0;k<EQUATIONS;k++)
		      sol[k] += f[l]*elsol[l*EQUATIONS+k];
		  delete []f;
		  NodeTemp = (Node*) HT_Node_Ptr->lookup(nodes+j*KEYLENGTH);
		  int* dof = NodeTemp->getdof();
		  if(NodeTemp->getinfo() != S_C_CON) {
		    fprintf ( fp, "%e %e %e %d %e %e %e %u %u %d %d %d %d %d %e %e %e %e \n", 
			      (*(NodeTemp->get_coord()))*(matprops)->LENGTH_SCALE,
			      (*(NodeTemp->get_coord()+1))*(matprops)->LENGTH_SCALE, 
			      NodeTemp->get_elevation() * (matprops->LENGTH_SCALE), 
			      myid, sol[0], sol[1], sol[2], *(EmTemp->pass_key()), 
			      *(EmTemp->pass_key()+1), EmTemp->get_gen(), 
			      EmTemp->get_which_son(), EmTemp->get_order(), 
			      *(EmTemp->get_elm_loc()), *(EmTemp->get_elm_loc()+1),
			      EmTemp->get_fluxIntegral(),*(EmTemp->get_el_error()),
			      *(EmTemp->get_el_error()+1),*(EmTemp->get_el_error()+2));

		  }
		  else { // S_C_CON will have a discontinuity in the elevation so fix that by interpolation
		    double elev;
		    int neighside, mynode;
		    if(j == EmTemp->get_which_son() + 1) {
		      mynode = j - 1;
		      neighside = j;
		    }
		    else if(j == EmTemp->get_which_son() - 1) {
		      mynode = j+1;
		      neighside = j-1;
		      if(neighside < 0)
			neighside = 3;
		    }
		    else if(EmTemp->get_which_son() == 0) {
		      mynode = 0;
		      neighside = 2;
		    }
		    else if(EmTemp->get_which_son() == 3) {
		      mynode = 3;
		      neighside = 0;
		    }
		    else {
			mynode = 1;
			neighside = 1;
		      //assert(0);
		    }
		    Node* NodeTemp2 = (Node*) HT_Node_Ptr->lookup(nodes+mynode*KEYLENGTH);
		    elev = .5 * NodeTemp2->get_elevation();
		    Element* EmTemp2 = (Element*) HT_Elem_Ptr->lookup((EmTemp->get_neighbors()+KEYLENGTH*neighside));
		    NodeTemp2 = (Node*) HT_Node_Ptr->lookup(EmTemp2->getNode()+j*KEYLENGTH);
		    elev += .5 * NodeTemp2->get_elevation();
		    fprintf ( fp, "%e %e %e %d %e %e %e %u %u %d %d %d %d %d %e %e %e %e \n", 
			      (*(NodeTemp->get_coord()))*(matprops)->LENGTH_SCALE,
			      (*(NodeTemp->get_coord()+1))*(matprops)->LENGTH_SCALE, 
			      elev * (matprops->LENGTH_SCALE), 
			      myid, sol[0], sol[1], sol[2], *(EmTemp->pass_key()), 
			      *(EmTemp->pass_key()+1), EmTemp->get_gen(), 
			      EmTemp->get_which_son(), EmTemp->get_order(), 
			      *(EmTemp->get_elm_loc()), *(EmTemp->get_elm_loc()+1),
			      EmTemp->get_fluxIntegral(),*(EmTemp->get_el_error()),
			      *(EmTemp->get_el_error()+1),*(EmTemp->get_el_error()+2));
		  }
		  
		  //  if( *(EmTemp->get_elm_loc())== 11 && *(EmTemp->get_elm_loc()+1)==17 ||  *(EmTemp->get_elm_loc())== 11 && *(EmTemp->get_elm_loc()+1)==12){
		  //    printf("(in tecplot.C) elmloc %d %d \n",*(EmTemp->get_elm_loc()),*(EmTemp->get_elm_loc()+1));
		  //    for(l=0;l<3;l++)
		  //     printf("corner %d solution[%d]= %e \n ",j,l,sol[l]);
		    // printf("(in tecplot.C)\n");
		  // }
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

  
} 
