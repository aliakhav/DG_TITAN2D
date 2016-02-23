#include "GisApi.h"
#include <stdlib.h>

/* command to run:
   ./preproc.x 5 /volcano/elchichon1/acbauer/grass.data/grass5 Colima ColimaR ColimaR  
   ./preproc.x 5 /home/acbauer/grass.data/grass5 Colima ColimaR ColimaR  
*/


/***************************************************************
 * main code
 **************************************************************/

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
  
  if(argc != 6) {
    printf("Need 5 command line arguments!\n");
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

      createfunky_(&xmax, &xmin, &ymax, &ymin, &number_y_grid_points);

      Delete_GIS_data();
    }
  else
    printf("Couldn't initialize the GIS information.\n");

  exit(0);
}

/***************************************************************
 * GIS stuff
 **************************************************************/
int write_SPR_header( int ncols, int nrows, double ulx, double uly, double resx, double resy, double nodata)
{
  fprintf (stdout, "GRIDREG\n" );
  fprintf (stdout, "INFO\n" );
  fprintf (stdout, "//Rectangular Grid SPRING format ASCII file\n" );
  fprintf (stdout, "//Format GRIDDEF  <ncols>  <nlins>  <X1>  <Y2>  <resX>  <resY>  <nodatavalue>\n" );
  fprintf (stdout, "GRIDDEF" );
  fprintf (stdout, " %d %d", ncols, nrows);
  fprintf (stdout, " %8.4f %8.4f", ulx, uly+10000000);
  fprintf (stdout, " %8.4f %8.4f", resx, resy);
  fprintf (stdout, " %2.1e\n", nodata);
  fprintf (stdout, "INFO_END\n" );
}

int write_SPR_end ()
{
  fprintf (stdout, "END\n" );
}
