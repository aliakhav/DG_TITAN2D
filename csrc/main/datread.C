#include "../header/hpfem.h"

extern "C" void cnst_();
extern "C" void setshape_();

//#ifdef DEBUG

void Read_data(HashTable** NodeTable, int myid, 
	       HashTable** ElemTable, MatProps* matprops_ptr,  
	       int* timesteps, double* maxtime, int* numoutput, 
	       int* adaptflag_ptr, int* viz_flag_ptr, int* order_flag_ptr)
{
  ifstream inscale("scale.data",ios::in);
  if(inscale.fail() == 0) {  
    inscale>>matprops_ptr->LENGTH_SCALE;
    inscale>>matprops_ptr->HEIGHT_SCALE;
    inscale>>matprops_ptr->GRAVITY_SCALE;
  }
  else{  // if this file doesn't exist, assume no scaling
    printf("Can't find scale.data for processor %d\n",myid);
    matprops_ptr->LENGTH_SCALE = 1;
    matprops_ptr->HEIGHT_SCALE = 1;
    matprops_ptr->GRAVITY_SCALE = 1;
  }
  matprops_ptr->epsilon = matprops_ptr->HEIGHT_SCALE/matprops_ptr->LENGTH_SCALE;

  inscale.close();
  double junk = 0;
  int NODE_TABLE_SIZE=400000;
  
  //char  filename[14] = "lsh4800xx.inp";
  char  filename[14] = "funkyxxxx.inp";
  unsigned min_key[KEYLENGTH];
  unsigned max_key[KEYLENGTH];
  double XRange[2];
  double YRange[2];
  unsigned key[KEYLENGTH];
  double coord[DIMENSION];
  Node* NodeP;
  int i, j, k;
  
  //set up for parabaloid pile
  ifstream inD2("simulation.data", ios::in);
  if(inD2.fail()) {
    printf("can't find simulation.data file\n");
    exit(0);
  }
  int numpiles = 0;  // the number of separate files to use
  inD2>>numpiles;
  double* parabaloid_pile_info = new double[numpiles*5];
  for(i=0;i<numpiles;i++) {
    inD2>>parabaloid_pile_info[5*i]; // pile height
    inD2>>parabaloid_pile_info[5*i+1]; // pile x-center
    inD2>>parabaloid_pile_info[5*i+2]; // pile y-center
    inD2>>parabaloid_pile_info[5*i+3]; // pile x-radius
    inD2>>parabaloid_pile_info[5*i+4]; // pile y-radius
  }
  inD2>>*timesteps;
  inD2>>*maxtime;
  inD2>>*numoutput;
  inD2>>*(adaptflag_ptr);
  inD2>>*viz_flag_ptr;
  inD2>>*order_flag_ptr;
  
  // read in GIS information
  char gis_main[300];
  char gis_sub[100];
  char gis_mapset[100];
  char gis_map[100];
  inD2>>gis_main;
  inD2>>gis_sub;
  inD2>>gis_mapset;
  inD2>>gis_map;
  inD2.close();
  i = Initialize_GIS_data(gis_main, gis_sub, gis_mapset, gis_map);
  if(i!= 0) {
    printf("Problem with GIS on processor %d\n",myid);
    exit(0);
  }
  
  // read in nodal data
  filename[5] = 48 + myid/1000;  
  filename[6] = 48 + (myid%1000)/100;
  filename[7] = 48 + (myid%100)/10;  
  filename[8] = 48 + myid%10;
  
  ifstream inDatafile(filename, ios::in);
  
  if(inDatafile.fail())
    {
      printf("Can't open file for %d \n", myid);
      //cerr << "Can't open file "
      //   << filename << endl;
      exit(0);
    }
  int Node_Num;
  inDatafile>>Node_Num;
  for(i=0;i<KEYLENGTH;i++) 
    inDatafile>>min_key[i];
  for(i=0;i<KEYLENGTH;i++) 
    inDatafile>>max_key[i];
  
  for(i=0; i < 2; i++)
    inDatafile>>XRange[i];
  
  for(i=0; i < 2; i++)
    inDatafile>>YRange[i];
  
  for(i=0;i<2;i++)
    XRange[i] = XRange[i]/matprops_ptr->LENGTH_SCALE;
  for(i=0;i<2;i++)
    YRange[i] = YRange[i]/matprops_ptr->LENGTH_SCALE;
  
  
  /*.....to tune the hashtable!!!.....*/
  
  // max_key[0]=4267889329;
  // min_key[0]=0;//43652745;//107232333;
  
  *NodeTable  = new HashTable(min_key, max_key, NODE_TABLE_SIZE, 2017, XRange, YRange); 
  
  for(i=0;i<Node_Num;i++)
    {
      double height;
      for(j=0;j<KEYLENGTH;j++) inDatafile >> key[j];
      for(j=0;j<DIMENSION;j++) inDatafile >> coord[j];
      for(j=0;j<2;j++)
	coord[j] = coord[j]/matprops_ptr->LENGTH_SCALE;
      NodeP = new Node(key, coord, matprops_ptr);
      (*NodeTable)->add(key, NodeP);     
    }   
  
  
  //done reading in node data
  //start reading in element data
  
  
  int EL_TABLE_SIZE=100000;
  
  //char  filename[14] = "lsh4800xx.inp";
  int material;
  int elm_loc[2];
  unsigned opposite_brother[KEYLENGTH];
  
  for(i=0;i<2;i++)
    {
      XRange[i]=0;
      YRange[i]=0;
    }
  
  Element* Quad9P;
  void*    p;
  int*     assocp;/*--*/
  unsigned*    keyP;
  
  unsigned nodes[9][2];
  unsigned neigh[4][2];
  int      neighbor_proc[4];
  
  int temp2;
  int interflag;
  
  int Elem_Num;
  inDatafile>>Elem_Num;  //--number of the elements assigned to the proc
  
  for(i=0;i<KEYLENGTH;i++)
    inDatafile>>min_key[i];
  for(i=0;i<KEYLENGTH;i++)
    inDatafile>>max_key[i];
  /*.....to tune the hashtable!!!.....*/
  
  // max_key[0]=4267122543;
  // min_key[0]=0;//43672387;//110456257;
  
  *ElemTable  = new HashTable(min_key, max_key, EL_TABLE_SIZE, 503, XRange, YRange); 
  
  for(i=0;i<Elem_Num;i++)
    {      
      for(j=0; j<9; j++) 
	for(k=0; k<KEYLENGTH; k++)
	  inDatafile>>nodes[j][k];
      
      interflag = 0;//---switch for interface
      for(j=0; j<4; j++)
	{
	  inDatafile>>temp2;//--read the neighbor info
	  if(temp2!=-1)//--if there is neighbor(-1 means the edge is bound)
	    {
	      neighbor_proc[j] = temp2;
	      if(neighbor_proc[j]!=myid) //--the neighbor belongs to other proc
		interflag = 1;//--switch is used for avoiding nominating neighbor twice
	      
	      for(k=0; k<KEYLENGTH; k++)
		inDatafile>>neigh[j][k];//--read the left parts of the key
	      
	    }
	  
	  else 
	    {
	      neighbor_proc[j]=-1;//--there is no neighbor 
	      for(k=0; k<KEYLENGTH; k++) neigh[j][k]=0;
	    }
	}
      
      //.....the essential boundary conditions....
      
      for(j=0; j < 4; j++)
	{
	  inDatafile>>temp2;
	  if(temp2!=-1)//--there is bound constraint
	    for(k=0;k<2;k++)
	      inDatafile>>junk;
	}
      
      //.....the natural boundary conditions.....
      for(j=0; j < 4; j++)
	{
	  inDatafile>>temp2;
	  if(temp2!=-1)//--there is bound constraint
	    {
	      for(k=0; k<2; k++)		     
		inDatafile>>junk;
	    }
	}
      inDatafile>>elm_loc[0];
      inDatafile>>elm_loc[1];
      inDatafile>>opposite_brother[0];
      inDatafile>>opposite_brother[1];
      
      inDatafile>>material;
      material++; //needs to be incremented 1 because of fortran

      //put in the pile height right here...
      // initial pile_height is linearly interpolated -- might want to change this if ORDER gets changed
      double pile_height[4] = {0,0,0,0};
      Node* ndtemp;
      double* ndcoord;
      // for multiple piles which may overlap, the highest value is used..
      for(k=0;k<4;k++) {
	ndtemp = (Node*) (*NodeTable)->lookup(&nodes[k][0]);
	ndcoord = ndtemp->get_coord();
	double height=0;
	for(j=0;j<numpiles;j++) {

	  //	if (((ndcoord[0]*matprops_ptr->LENGTH_SCALE-
	  //	     parabaloid_pile_info[5*j+1])/parabaloid_pile_info[5*j+3])*
	  //	    ((ndcoord[0]*matprops_ptr->LENGTH_SCALE-
	  //	     parabaloid_pile_info[5*j+1])/parabaloid_pile_info[5*j+3])+
	  //	    ((ndcoord[1]*matprops_ptr->LENGTH_SCALE-
	  //      parabaloid_pile_info[5*j+2])/parabaloid_pile_info[5*j+4])*
	  //	    ((ndcoord[1]*matprops_ptr->LENGTH_SCALE-
	  //     parabaloid_pile_info[5*j+2])/parabaloid_pile_info[5*j+4])
	  //   <1) {
	  //  height = parabaloid_pile_info[5*j]/matprops_ptr->HEIGHT_SCALE;
	  //}
	  height = parabaloid_pile_info[5*j]*
	    (1.-pow((ndcoord[0]*matprops_ptr->LENGTH_SCALE-
		     parabaloid_pile_info[5*j+1])/
		    parabaloid_pile_info[5*j+3],2)-
	     pow((ndcoord[1]*matprops_ptr->LENGTH_SCALE-
		  parabaloid_pile_info[5*j+2])/
		 parabaloid_pile_info[5*j+4],2))/matprops_ptr->HEIGHT_SCALE;

	  if(pile_height[k] < height)
	    pile_height[k] = height;
	}
      }
      
      //      if (pile_height[0]>0)
      //       printf("nonzero pile height %e\n", pile_height[0]); 

      double pile_height2[3] = {0,0,0};
      pile_height2[0] = .25*(pile_height[0]+pile_height[1]+pile_height[2]+pile_height[3]);
      pile_height2[1] = .25*(pile_height[1]+pile_height[2]-pile_height[0]-pile_height[3]);
      pile_height2[2] = .25*(pile_height[2]+pile_height[3]-pile_height[1]-pile_height[0]);
      if((dabs(pile_height2[1])+dabs(pile_height2[2])) > pile_height2[0]) {
	remove_neg_pile_height(pile_height, pile_height2[0], (pile_height2+1), (pile_height2+2));
      }
      
      Quad9P=new Element(nodes, neigh, neighbor_proc, material, pile_height2, myid,
			 elm_loc, opposite_brother);
      (*ElemTable)->add(nodes[8], Quad9P);
    }
  
  //read in material properties
  int material_count = 0;
  inDatafile>>material_count;
  double* lambda = new double[material_count];
  double* mi = new double[material_count];
  for(i=0;i<material_count;i++)
    inDatafile>>lambda[i]>>mi[i];

  matprops_ptr->intfrict = lambda[0]*PI/180.;
  matprops_ptr->bedfrict = mi[0]*PI/180.;
  inDatafile.close();
  
//initialize element stiffness routines info
  cnst_();
  setshape_();
  delete []parabaloid_pile_info;
  delete []lambda;
  delete []mi;
  
  
}
