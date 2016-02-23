#ifndef GIS_API_H
#define GIS_API_H

//#include "GisRasterHdr.h"

typedef struct
{
	double zmin;
	double zmax;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double resolution;
	int nrows;
	int ncols;
	double wxmin;
	double wxmax;
	double wymin;
	double wymax;
	double wresolution;
	int wnrows;
	int wncols;
	float** elev;
	float** xslope;
	float** yslope;
	float** slope;
	float** xcurv;
	float** ycurv;
} Gis_Grid;

#define G_API_MAXFLOAT 1.0e31
#define G_API_BIGFLOAT 1.0e32

/***************************************************************/
/* Error codes:
	 0 -> everything okay
	-1 -> resolution finer than available information
	-2 -> requested (x,y) location is not in the data set
	-3 -> memory problem
	-4 -> something else wrong
*/
/***************************************************************/

/***************************************************************/
/* SELECTION OF DATA */
/***************************************************************/
int Initialize_GIS_data( char* GISDbase, char* location, char* mapset, char* raster_file);
/*
Select mapset at location in database
Input:
	GISDbase - Full path of database
	location - Location
	mapset - Mapset
	raster_file - Raster map name
Output:
Return:
	0 if error, 1 otherwise
*/

int Delete_GIS_data();
/*
Clears all GIS_data
Input:
Output:
Return:
	0 if error, 1 otherwise
*/

/*BASIC INFORMATION RECOVERY */
int Get_xmax(double resolution, double* xmax);
/*
Return extents of original grid
Input:
	resolution - resolution
Output:
	xmax - maximum X coordinate of original grid
Return:
	0 if OK, see table otherwise
*/

int Get_xmin(double resolution, double* xmin);
/*
Return extents of original grid
Input:
	resolution - resolution
Output:
	xmin - minimum X coordinate of original grid
Return:
	0 if OK, see table otherwise
*/

int Get_ymax(double resolution, double* ymax);
/*
Return extents of original grid
Input:
	resolution - resolution
Output:
	ymax - maximum Y coordinate of original grid
Return:
	0 if OK, see table otherwise
*/

int Get_ymin(double resolution, double* ymin);
/*
Return extents of original grid
Input:
	resolution - resolution
Output:
	ymin - minimum Y coordinate of original grid
Return:
	0 if OK, see table otherwise
*/

int Get_elev_min(double resolution, double* elevmin);
/*
Return minimum elevation of original grid
Input:
	resolution - resolution
Output:
	elevmin - minimum elevation of original grid
Return:
	0 if OK, see table otherwise
*/

int Get_elev_max(double resolution, double* elevmax);
/*
Return maximum elevation of original grid
Input:
	resolution - resolution
Output:
	elevmax - maximum elevation of original grid
Return:
	0 if OK, see table otherwise
*/

int Get_window(double *xmin, double *xmax, double *ymin, double *ymax);
/*
Return extents of active window
Input:
Output:
	xmin - minimum X coordinate of active region
	xmax - maximum X coordinate of active region
	ymin - minimum Y coordinate of active region
	ymax - maximum Y coordinate of active region
Return:
	0 if error, 1 otherwise
*/

int Get_max_resolution(double* resolution);
/*
Return original resolution of grid
Input:
Output:
	resolution - resolution
Return:
	0 if OK, see table otherwise
*/

int Get_number_of_rows(int *rows);
/*
Return number of rows of original grid
Input:
Output:
	rows - number of rows of original grid
Return:
	0 if error, 1 otherwise
*/

int Get_number_of_columns(int *cols);
/*
Return number of columns of original grid
Input:
Output:
	cols - number of columns of original grid
Return:
	0 if error, 1 otherwise
*/
/***************************************************************/

/***************************************************************/
/*GETTING VALUES FOR SINGLE POINTS*/
int Get_elevation(double resolution, double x, double y, double* elev);
/*
Return elevation at point XY of original grid
Input:
	resolution - resolution
	x - point X coordinate
	y - Point Y coordinate
Output:
	elev - elevation at point XY of original grid
Return:
	0 if OK, see table otherwise
*/

int Get_slope(double resolution, double x, double y, double* xslope, double* yslope);
/*
Return slope at point XY of original grid
Input:
	resolution - resolution
	x - point X coordinate
	y - Point Y coordinate
Output:
	xslope - slope at point XY of original grid in X direction
	yslope - slope at point XY of original grid in Y direction
Return:
	0 if OK, see table otherwise
*/

int Get_curvature(double resolution, double x, double y, double* xcurv, double* ycurv);
/*
Return curvature at point XY of original grid
Input:
	resolution - resolution
	x - point X coordinate
	y - Point Y coordinate
Output:
	xcurv - curvature at point XY of original grid in X direction
	ycurv - curvature at point XY of original grid in Y direction
Return:
	0 if OK, see table otherwise
*/

/***************************************************************/
/* GETTING VALUES FOR MULTIPLE POINTS */
/***************************************************************/
int Get_elevation_array(double* resolution, double* x, double* y, double* elev, int number_of_locations);
/*
Return elevation at points XY of original grid
Input:
	resolution - resolution
	x - array with points X coordinate
	y - array with points Y coordinate
	number_of_locations - number of points
Output:
	elev - elevation at points XY of original grid
Return:
	0 if OK, see table otherwise
*/

int Get_slope_array(double* resolution, double* x, double* y, double* xslope, double* yslope, int number_of_locations);
/*
Return slope at points XY of original grid
Input:
	resolution - resolution
	x - array with points X coordinate
	y - array with points Y coordinate
	number_of_locations - number of points
Output:
	xslope - slopes at points XY of original grid in X direction
	yslope - slopes at points XY of original grid in Y direction
Return:
	0 if OK, see table otherwise
*/

int Get_curvature_array(double* resolution, double* x, double* y, double* xcurv, double* ycurv, int number_of_locations);
/*
Return curvature at points XY of original grid
Input:
	resolution - resolution
	x - array with points X coordinate
	y - array with points Y coordinate
	number_of_locations - number of points
Output:
	xcurv - curvature at point XY of original grid in X direction
	ycurv - curvature at point XY of original grid in Y direction
Return:
	0 if OK, see table otherwise
*/

int Get_elevation_grid(double resolution, double xmin, double xmax, double ymin, double ymax, double* elev);
/*
Return elevation at points XY of original grid
Input:
	resolution - resolution
	xmin - minimum X coordinate of points window
	xmax - maximum X coordinate of points window
	ymin - minimum Y coordinate of points window
	ymax - maximum Y coordinate of points window
Output:
	elev - elevation at points XY of selected window
Return:
	0 if OK, see table otherwise
*/

int Get_slope_grid(double resolution, double xmin, double xmax, double ymin, double ymax, double* slope);
/*
Return elevation at points XY of original grid
Input:
	resolution - resolution
	xmin - minimum X coordinate of points window
	xmax - maximum X coordinate of points window
	ymin - minimum Y coordinate of points window
	ymax - maximum Y coordinate of points window
Output:
	elev - elevation at points XY of selected window
Return:
	0 if OK, see table otherwise
*/

int Get_curvature_grid(double resolution, double xmin, double xmax, double ymin, double ymax, double* xcurv, double* ycurv);
/*
Return elevation at points XY of original grid
Input:
	resolution - resolution
	xmin - minimum X coordinate of points window
	xmax - maximum X coordinate of points window
	ymin - minimum Y coordinate of points window
	ymax - maximum Y coordinate of points window
Output:
	elev - elevation at points XY of selected window
Return:
	0 if OK, see table otherwise
*/

/***************************************************************/
/* INTERNAL USE FUNCTIONS */
/***************************************************************/
//int set_from_header(GisRasterHdr& gisHeader);

float **alloc_float_matrix( int rows, int cols );
/* Allocate matrix of floats
Input:
	nrows - matrix number of rows
	ncols - matrix number of columns
Output:
Return:
	pointer to allocated space, 0 if error.
*/

int free_float_matrix( float **m );
/* Free matrix of floats
Input:
	m - matrix of floats
Output:
Return:
*/
int clear_gis_grid();
/* Clear gis_grid structure and frees matrices of floats */

int print_grid();
/* Utility to print gis grid contents on stdout */

int calculate_slope();
/* Calculate slope from elevation grid */

int calculate_curvature();
/* Calculate curvature from elevation grid */

void find_min_max();
/* Find elevation grid maximum and minimum */

void get_grid(double resolution, double xmin, double xmax, double ymin, double ymax, float** ingrid, double* outgrid);
/* Copy any grid to output vector
	resolution - resolution
	xmin - minimum X coordinate of points window
	xmax - maximum X coordinate of points window
	ymin - minimum Y coordinate of points window
	ymax - maximum Y coordinate of points window
	ingrid - input grid
Output:
	outgrid - output vector
*/

#endif



