#ifndef NODE_C
#define NODE_C
#include "../header/node.h"
#include "../header/GisApi.h"
#include "../header/properties.h"
#include <mpi.h>


Node:: Node() {
  printf("creating a node without setting its values\n");
}

Node:: Node(unsigned* keyi, double* coordi, MatProps* matprops_ptr)
{
  int      i;
  nextptr =0;
  preptr = 0;
  
  info   = INIT;
	 
  for(i=0; i<DIMENSION; i++)
    coord[i] = *(coordi+i);
  
  for(i=0; i<KEYLENGTH; i++)
    key[i] = *(keyi+i);
  dof[0] = INIT;
  dof[1] = INIT;
  // find the max resolution of the GIS info and then get the elevation at this node
  double resolution = 0;
  double xcoord = coord[0]*(matprops_ptr->LENGTH_SCALE);
  double ycoord = coord[1]*(matprops_ptr->LENGTH_SCALE);
  i = Get_max_resolution(&resolution);
  if(i != 0) {
    printf("error in Get_max_resolution\n");
    exit(1);
  }    
  i = Get_elevation(resolution, xcoord, ycoord, &elevation);
  if(i != 0) {
    printf("error in Get_elevation\n");
    exit(1);
  }    
  elevation = elevation/matprops_ptr->LENGTH_SCALE;
}

Node::Node(unsigned* keyi, double* coordi, int inf, MatProps* matprops_ptr)  //for refined
{
  int      i;
  
  nextptr =0;
  preptr = 0;
    info   = INIT;
  
  for(i=0; i<DIMENSION; i++)
    coord[i] = *(coordi+i);
  
  for(i=0; i<KEYLENGTH; i++)
    key[i] = *(keyi+i);
  dof[0] = INIT;
  dof[1] = INIT;
  
  info=inf;

  // find the max resolution of the GIS info and then get the elevation at this node
  double resolution = 0;
  i = Get_max_resolution(&resolution);
  if(i != 0) {
    printf("error in Get_max_resolution\n");
    exit(1);
  }    
  double xcoord = coord[0]*(matprops_ptr->LENGTH_SCALE);
  double ycoord = coord[1]*(matprops_ptr->LENGTH_SCALE);
  i = Get_elevation(resolution, xcoord, ycoord, &elevation);
  if(i != 0) {
    printf("error in Get_elevation\n");
    exit(1);
  }  
  elevation = elevation/matprops_ptr->LENGTH_SCALE;

  return;
}

Node:: Node(unsigned* keyi, double* coordi, int inf, 
	    double elev)
{
  int      i;
  
  nextptr =0;
  preptr = 0;
  info   = INIT;
  
  for(i=0; i<DIMENSION; i++)
    coord[i] = *(coordi+i);
  
  for(i=0; i<KEYLENGTH; i++)
    key[i] = *(keyi+i);
  dof[0] = INIT;
  dof[1] = INIT;
  
  info=inf;
  elevation = elev;

  return;
}

void Node::putinfo(int in)
{
  info = in;
/*  if(key[0] == (unsigned) 2962355296) {
    int mmmyid;
    MPI_Comm_rank(MPI_COMM_WORLD, &mmmyid);
    printf("?????????????????????????????????????????????????????? \n");
    printf("?????????????????????????????????????????????????????? \n");
    printf("changing info %u %u on %d $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n",key[0], key[1], mmmyid);
    printf("?????????????????????????????????????????????????????? \n");
    printf("?????????????????????????????????????????????????????? \n");
    
  }*/
}


Node:: ~Node() { }

#endif





