#include "../../csrc/header/GisApi.h"
#include <stdio.h>
#include <stdlib.h>

/* command to run:
   ./preproc.x 5 /volcano/elchichon1/acbauer/grass.data/grass5 Colima ColimaR ColimaR  
   ./preproc.x 5 /home/acbauer/grass.data/grass5 Colima ColimaR ColimaR  
*/


/***************************************************************
 * main code
 **************************************************************/
extern "C" {
	void* createfunky_(double*, double*, double*, double*, int*); }

int main( int argc, char ** argv )
{
  double xmin, xmax, ymin, ymax;
  double res;
  int nrows, ncols;
  int row, col, i;
  double *grid = 0;
  int number_y_grid_points = 1;
  
/*  printf("%d\n",argc);
  for (i=0; i<argc; i++)
    printf("%s\n",argv[i]);*/
  
  if(argc != 6 && argc != 10) {
    printf(" argc is %d\n", argc);
    printf("Need either 5 or 9 command line arguments!\n");
    exit(0);
  }
  
  
  number_y_grid_points = atoi(argv[1]);
    
  
  if ( ! Initialize_GIS_data( argv[2], argv[3], argv[4], argv[5]) )
    {
      Get_max_resolution (&res);
      Get_xmax (res, &xmax);
      Get_xmin (res, &xmin);
      Get_ymax (res, &ymax);
      Get_ymin (res, &ymin);
      
      printf("The GIS region is:  x = %e to %e   y = %e to %e.\n",xmin, xmax, ymin, ymax);

      if(argc == 10) {
	double min_max_coords[4];
	for(i=0;i<4;i++)
	  min_max_coords[i] = atof(argv[i+6]);
	/* check the values for user input min and max coordinates ... */
	if(min_max_coords[0] > xmin && min_max_coords[0] < xmax)
	  xmin = min_max_coords[0];
	if(min_max_coords[1] > ymin && min_max_coords[1] < ymax)
	  ymin = min_max_coords[1];
	if(min_max_coords[2] > xmin && min_max_coords[2] < xmax)
	  xmax = min_max_coords[2];
	if(min_max_coords[3] > ymin && min_max_coords[3] < ymax)
	  ymax = min_max_coords[3];
      }
      printf("The simulation region is:  x = %e to %e   y = %e to %e.\n",xmin, xmax, ymin, ymax);

      createfunky_(&xmax, &xmin, &ymax, &ymin, &number_y_grid_points);

      Delete_GIS_data();
    }
  else
    printf("Couldn't initialize the GIS information.\n");

  exit(0);
}

