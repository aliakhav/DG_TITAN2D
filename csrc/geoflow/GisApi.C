#include "../header/GisApi.h"
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>

#include "../header/GisBinFile.h"
#include "../header/GisRasterHdr.h"

Gis_Grid gis_grid;

int set_from_header(GisRasterHdr& gisHeader);

/***************************************************************/
/* SELECTION OF DATA */
/***************************************************************/
int Initialize_GIS_data( char* GISDbase, char* location, char* mapset, char* raster_file)
{
	int nrows, ncols;
	int row;

	char gisPath[200];
	char gisFullPath[250];
#if defined WIN32
	char* gisSlash = "\\";
#else
	char* gisSlash = "/";
#endif

	clear_gis_grid();

	if ( GISDbase && location && mapset && raster_file )
	{
		strcpy (gisPath, GISDbase);
		sprintf(gisPath,"%s%s%s%s",gisPath,gisSlash,location,gisSlash);
		sprintf(gisPath,"%s%s%s",gisPath,mapset,gisSlash);
		strcpy(gisFullPath,gisPath);
		sprintf(gisFullPath,"%scellhd%s%s",gisFullPath,gisSlash,raster_file);

		GisRasterHdr gisHeader (gisFullPath);
		if ( ( ! gisHeader.good() ) ||
			 ( set_from_header(gisHeader) != 0 ) )
				return -4;

		nrows = gisHeader.Rows();
		ncols = gisHeader.Cols();

		if ( nrows < 1 || ncols < 1 )
			return -4;

		if ( ! ( gis_grid.elev = alloc_float_matrix ( nrows, ncols ) ) )
			return -3;	// memory error

		strcpy(gisFullPath,gisPath);
		sprintf(gisFullPath,"%sfcell%s%s",gisFullPath,gisSlash,raster_file);

		GisBinFile binFile (gisFullPath);
		if ( binFile.good() )
		{
			binFile.setEndian ("big");
			binFile.setWordSize (4);
			binFile.setIsInteger (false);
			binFile.isCompressed (gisHeader.isCompressed());
			binFile.nRows(gisHeader.Rows());
			binFile.nCols(gisHeader.Cols());

			for (row = 0; row < nrows; row++)
			{
				if ( ! binFile.readRow (row, gis_grid.elev[row]) )
				{
					return -4;
				}
			}
			gis_grid.zmin = G_API_BIGFLOAT;
			gis_grid.zmax = -G_API_BIGFLOAT;
			return 0;
		}
	}
	return -4;
}

int Delete_GIS_data()
{
	clear_gis_grid();
	return 0;
}

/***************************************************************/
/*BASIC INFORMATION RECOVERY */
/***************************************************************/
int Get_xmax(double resolution, double* xmax)
{
	if ( gis_grid.resolution > resolution )
		return -1;
	*xmax = gis_grid.xmax;
	return 0;
}

int Get_xmin(double resolution, double* xmin)
{
	if ( gis_grid.resolution > resolution )
		return -1;
	*xmin = gis_grid.xmin;
	return 0;
}

int Get_ymax(double resolution, double* ymax)
{
	if ( gis_grid.resolution > resolution )
		return -1;
	*ymax = gis_grid.ymax;
	return 0;
}

int Get_ymin(double resolution, double* ymin)
{
	if ( gis_grid.resolution > resolution )
		return -1;
	*ymin = gis_grid.ymin;
	return 0;
}

int Get_elev_min(double resolution, double* elevmin)
{
	if ( gis_grid.zmin > gis_grid.zmax )
		find_min_max();

	*elevmin = gis_grid.zmin;
	return 0;
}

int Get_elev_max(double resolution, double* elevmax)
{
	if ( gis_grid.zmin > gis_grid.zmax )
		find_min_max();

	*elevmax = gis_grid.zmax;
	return 0;
}

int Get_window(double* xmin, double* xmax, double* ymin, double* ymax)
{
	if ( gis_grid.xmin > 0. && gis_grid.xmax > gis_grid.xmin &&
		 gis_grid.ymin > 0. && gis_grid.ymax > gis_grid.ymin )
	{
		*xmin = gis_grid.xmin;
		*xmax = gis_grid.xmax;
		*ymin = gis_grid.ymin;
		*ymax = gis_grid.ymax;
		return 0;
	}
	return -4;
}

int Get_max_resolution(double* resolution)
{
	*resolution = gis_grid.resolution;
	return 0;
}

int Get_number_of_rows(int *rows)
{
	*rows = gis_grid.nrows;
	return 0;
}

int Get_number_of_columns(int *cols)
{
	*cols = gis_grid.ncols;
	return 0;
}

/***************************************************************/
/*GETTING VALUES FOR SINGLE POINTS*/
/***************************************************************/
double interpolate_bilinear_at ( double resolution, double x, double y, float** ingrid )
{
  double dx, dy, dx1, dy1, p1, p2;
  int row, col;

  col = (int)x;
  row = (int)y;

  dx  = x - (double)col;
  dy  = y - (double)row;
  dx1 = 1.0 - dx;
  dy1 = 1.0 - dy;
  p1  = ingrid[row  ][col  ] * dy1 +
	ingrid[row+1][col  ] * dy;
  p2  = ingrid[row  ][col+1] * dy1 +
        ingrid[row+1][col+1] * dy;
  return ( p1 * dx1 + p2 * dx );
}

int Get_elevation(double resolution, double x, double y, double* elev)
{
	if ( x >= gis_grid.xmin && x <= gis_grid.xmax &&
		 y >= gis_grid.ymin && y <= gis_grid.ymax )
	{
	  x = ( x - gis_grid.xmin ) / gis_grid.resolution;
	  y = ( gis_grid.ymax - y ) / gis_grid.resolution;
	  if ( (int)y >= gis_grid.nrows-1 ||
	       (int)x >= gis_grid.ncols-1 )
	    *elev = gis_grid.elev[(int)y-1][(int)x-1];
	  else
	    *elev = interpolate_bilinear_at ( gis_grid.resolution, x, y, gis_grid.elev );
		//		*elev = gis_grid.elev[row][col];
	}
	return 0;
}

int Get_slope(double resolution, double x, double y, double* xslope, double* yslope)
{
	int status;

	if ( gis_grid.xslope == 0 )
	{
		status = calculate_slope();
		if ( status != 0 )
			return status;
	}

	if ( x >= gis_grid.xmin && x <= gis_grid.xmax &&
		 y >= gis_grid.ymin && y <= gis_grid.ymax )
	{
	  //		row = ( gis_grid.ymax - y ) / gis_grid.resolution;
	  //		col = ( x - gis_grid.xmin ) / gis_grid.resolution;
	  //		*xslope = gis_grid.xslope[row][col];
	  //		*yslope = gis_grid.yslope[row][col];
	  x = ( x - gis_grid.xmin ) / gis_grid.resolution;
	  y = ( gis_grid.ymax - y ) / gis_grid.resolution;
	  if ( (int)y >= gis_grid.nrows-1 ||
	       (int)x >= gis_grid.ncols-1 )
	    {
	       *xslope = gis_grid.xslope[(int)y-1][(int)x-1];
	       *yslope = gis_grid.yslope[(int)y-1][(int)x-1];
	    }
	  else
	    {
	       *xslope = interpolate_bilinear_at ( gis_grid.resolution, x, y, gis_grid.xslope );
	       *yslope = interpolate_bilinear_at ( gis_grid.resolution, x, y, gis_grid.yslope );
	    }

	}
	return 0;
}

int Get_curvature(double resolution, double x, double y, double* xcurv, double* ycurv)
{
	int row, col, i, j;
	int status;
	float** xxcurv;
	float** yycurv;
	double x1, y1;

	xxcurv = alloc_float_matrix(2,2);
	yycurv = alloc_float_matrix(2,2);
	if ( !xxcurv || !yycurv)
	  return -3;

	if ( gis_grid.xslope == 0 )
	{
		status = calculate_slope();
		if ( status != 0 )
			return status;
	}

	if ( x >= gis_grid.xmin && x <= gis_grid.xmax &&
		 y >= gis_grid.ymin && y <= gis_grid.ymax )
	{
		row = (int)(( gis_grid.ymax - y ) / gis_grid.resolution);
		if ( row == 0 )
			row++;
		if ( row == gis_grid.nrows - 1 )
			row--;
		col = (int)(( x - gis_grid.xmin ) / gis_grid.resolution);
		if ( col == 0 )
			col++;
		if ( col == gis_grid.ncols - 1 )
			col--;
		if ( col >= (gis_grid.ncols - 3) || row >= gis_grid.nrows - 3 )
		{
		    *xcurv = ( ( gis_grid.xslope[row-1][col+1] - gis_grid.xslope[row-1][col-1] ) +
			   2 * ( gis_grid.xslope[row  ][col+1] - gis_grid.xslope[row  ][col-1] ) +
			       ( gis_grid.xslope[row+1][col+1] - gis_grid.xslope[row+1][col-1] ) ) /
			   ( 8 * gis_grid.resolution );

		    *ycurv = ( ( gis_grid.yslope[row-1][col-1] - gis_grid.yslope[row+1][col-1] ) +
			   2 * ( gis_grid.yslope[row-1][col  ] - gis_grid.yslope[row+1][col  ] ) +
			       ( gis_grid.yslope[row-1][col+1] - gis_grid.yslope[row+1][col+1] ) ) /
			   ( 8 * gis_grid.resolution );
		}
		else
		{
		  for (i = 0; i < 2; i++)
		    {
		      for (j = 0; j < 2; j++)
			{
			  xxcurv[i][j] = ( ( gis_grid.xslope[row+i-1][col+j+1] - gis_grid.xslope[row+i-1][col+j-1] ) +
			   2 * ( gis_grid.xslope[row+i][col+j+1] - gis_grid.xslope[row+i][col+j-1] ) +
			       ( gis_grid.xslope[row+i+1][col+j+1] - gis_grid.xslope[row+i+1][col+j-1] ) ) /
			   ( 8 * gis_grid.resolution );
			  yycurv[i][j] = ( ( gis_grid.yslope[row+i-1][col+j-1] - gis_grid.yslope[row+i+1][col+j-1] ) +
			   2 * ( gis_grid.yslope[row+i-1][col+j] - gis_grid.yslope[row+i+1][col+j] ) +
			       ( gis_grid.yslope[row+i-1][col+j+1] - gis_grid.yslope[row+i+1][col+j+1] ) ) /
			   ( 8 * gis_grid.resolution );
			}
		    }
		    x1 = x - gis_grid.xmin - gis_grid.resolution*(double)col;
		    y1 = gis_grid.ymax - gis_grid.resolution*(double)row - y;
		    x1 = x1/gis_grid.resolution;
		    y1 = y1/gis_grid.resolution;
		    *xcurv = interpolate_bilinear_at ( gis_grid.resolution, x1, y1, xxcurv );
	            *ycurv = interpolate_bilinear_at ( gis_grid.resolution, x1, y1, yycurv );

		}
	}
	free_float_matrix (xxcurv);
	free_float_matrix (yycurv);
	return 0;
}

/***************************************************************/
/* GETTING VALUES FOR MULTIPLE POINTS */
/***************************************************************/
int Get_elevation_array(double* resolution, double* x, double* y, double* elev, int number_of_locations)
{
	double x1, y1;
	int i;
	for ( i = 0; i < number_of_locations; i++ )
	{
		if ( x[i] >= gis_grid.xmin && x[i] <= gis_grid.xmax &&
			y[i] >= gis_grid.ymin && y[i] <= gis_grid.ymax )
		{
//			row = ( gis_grid.ymax - y[i] ) / gis_grid.resolution;
//			col = ( x[i] - gis_grid.xmin ) / gis_grid.resolution;
//			elev[i] = gis_grid.elev[row][col];
		  x1 = ( x[i] - gis_grid.xmin ) / gis_grid.resolution;
		  y1 = ( gis_grid.ymax - y[i] ) / gis_grid.resolution;
		  if ( (int)x1 >= gis_grid.nrows-1 ||
		       (int)y1 >= gis_grid.ncols-1 )
		    elev[i] = gis_grid.elev[(int)y1-1][(int)x1-1];
		  else
		    elev[i] = interpolate_bilinear_at ( gis_grid.resolution, x1, y1, gis_grid.elev);
		}
	}
	return 0;
}

int Get_slope_array(double* resolution, double* x, double* y, double* xslope, double* yslope, int number_of_locations)
{
	int row, col;
	int i, status;

	if ( gis_grid.xslope == 0 )
	{
		status = calculate_slope();
		if ( status != 0 )
			return status;
	}

	for ( i = 0; i < number_of_locations; i++ )
	{
		if ( x[i] >= gis_grid.xmin && x[i] <= gis_grid.xmax &&
			 y[i] >= gis_grid.ymin && y[i] <= gis_grid.ymax )
		{
			row = (int)(( gis_grid.ymax - y[i] ) / gis_grid.resolution);
			col = (int)(( x[i] - gis_grid.xmin ) / gis_grid.resolution);
			xslope[i] = gis_grid.xslope[row][col];
			yslope[i] = gis_grid.yslope[row][col];
		}
	}
	return 0;
}

int Get_curvature_array(double* resolution, double* x, double* y, double* xcurv, double* ycurv, int number_of_locations)
{
	int row, col;
	int i, status;

	if ( gis_grid.slope == 0 )
	{
		status = calculate_slope();
		if ( status != 0 )
			return status;
	}
	if ( gis_grid.xcurv == 0 )
	{
		status = calculate_curvature();
		if ( status != 0 )
			return status;
	}
	for ( i = 0; i < number_of_locations; i++ )
	{
		if ( x[i] >= gis_grid.xmin && x[i] <= gis_grid.xmax &&
			 y[i] >= gis_grid.ymin && y[i] <= gis_grid.ymax )
		{
			row = (int)(( gis_grid.ymax - y[i] ) / gis_grid.resolution);
			col = (int)(( x[i] - gis_grid.xmin ) / gis_grid.resolution);
			xcurv[i] = gis_grid.xcurv[row][col];
			ycurv[i] = gis_grid.ycurv[row][col];
		}
	}
	return 0;
}

int Get_elevation_grid(double resolution, double xmin, double xmax, double ymin, double ymax, double* elev)
{
	get_grid( resolution, xmin, xmax, ymin, ymax, gis_grid.elev, elev);
	return 0;
}

int Get_slope_grid(double resolution, double xmin, double xmax, double ymin, double ymax, double* slope)
{
	int status;

	if ( gis_grid.slope == 0 )
	{
		status = calculate_slope();
		if ( status != 0 )
			return status;
	}

	get_grid( resolution, xmin, xmax, ymin, ymax, gis_grid.slope, slope);
	return 0;
}

int Get_curvature_grid(double resolution, double xmin, double xmax, double ymin, double ymax, double* xcurv, double* ycurv)
{
	int status;

	if ( gis_grid.slope == 0 )
	{
		status = calculate_slope();
		if ( status != 0 )
			return status;
	}

	if ( gis_grid.xcurv == 0 )
	{
		status = calculate_curvature();
		if ( status != 0 )
			return status;
	}

	get_grid( resolution, xmin, xmax, ymin, ymax, gis_grid.xcurv, xcurv);
	get_grid( resolution, xmin, xmax, ymin, ymax, gis_grid.ycurv, ycurv);
	
	return 0;
}

/***************************************************************/
/* INTERNAL USE FUNCTIONS */
/***************************************************************/

void get_grid(double resolution, double xmin, double xmax, double ymin, double ymax, float** ingrid, double* outgrid)
{
	int row, col, i, j, li;
	int irow, icol;			/* initial row and column */
	int frow, fcol;			/* final row and column */
	double icold, irowd;
	double resr;

	icold = (xmin - gis_grid.xmin) / gis_grid.resolution;
	icol  = (int) icold;

	irowd = (gis_grid.ymax - ymax) / gis_grid.resolution;
	irow  = (int) irowd;

	fcol  = (int) ( (xmax - gis_grid.xmin) / gis_grid.resolution);
	frow  = (int) ( (gis_grid.ymax - ymin) / gis_grid.resolution);

	resr = resolution/gis_grid.resolution;

	li = 0;
	j = 0;
	row = irow;
	if (fabs(resr - 1.) < 0.0000001)
	{
		while ( row < frow )
		{
			col = icol;
			i = 0;
			while ( col < fcol )
			{
				outgrid[li++] = ingrid[row][col++];
			}
			row++;
		}
	}
	else
	{
		while ( row < frow )
		{
			col = icol;
			i = 0;
			while ( col < fcol )
			{
				outgrid[li++] = ingrid[row][col];
				col = (int) (icold + (double)(++i) * resr);
			}
			row = (int) (irowd  + (double)(++j) * resr);
		}
	}
}

int clear_gis_grid()
{
	gis_grid.zmin = 0;
	gis_grid.zmax = 0;
	gis_grid.xmin = 0;
	gis_grid.xmax = 0;
	gis_grid.ymin = 0;
	gis_grid.ymax = 0;
	gis_grid.resolution = 0;
	gis_grid.wxmin = 0;
	gis_grid.wxmax = 0;
	gis_grid.wymin = 0;
	gis_grid.wymax = 0;
	gis_grid.wresolution = 0;

	free_float_matrix ( gis_grid.elev );
	free_float_matrix ( gis_grid.xslope );
	free_float_matrix ( gis_grid.yslope );
	free_float_matrix ( gis_grid.slope );
	free_float_matrix ( gis_grid.xcurv );
	free_float_matrix ( gis_grid.ycurv );
    return 0;
}

float **alloc_float_matrix( int nrows, int ncols )
{
    float **m=0;
    int i;

    if ( m = (float **) calloc (nrows, sizeof(float *)) )
	{
		if ( m[0] = (float *) calloc (nrows*ncols, sizeof(float)) )
		{
			for (i = 1; i < nrows; i++)
				m[i] = m[i-1] + ncols;
		}
		else
			return 0;
	    return m;
	}
	return 0;
}

int free_float_matrix( float **m )
{
 	if ( m )
	{
		if ( m[0] )
			free ( m[0] );
		free (m);
	}
    return 0;
}

int calculate_slope()
{
	int row, col;

	free_float_matrix ( gis_grid.xslope );
	free_float_matrix ( gis_grid.yslope );
	free_float_matrix ( gis_grid.slope );

	free_float_matrix ( gis_grid.xcurv );
	free_float_matrix ( gis_grid.ycurv );

	if ( ! ( gis_grid.xslope = alloc_float_matrix ( gis_grid.nrows, gis_grid.ncols ) ) )
		return -3;	/*memory error*/

	if ( ! ( gis_grid.yslope = alloc_float_matrix ( gis_grid.nrows, gis_grid.ncols ) ) )
		return -3;	/*memory error*/

	if ( ! ( gis_grid.slope = alloc_float_matrix ( gis_grid.nrows, gis_grid.ncols ) ) )
		return -3;	/*memory error*/

	for (row = 1; row < gis_grid.nrows - 1; row++)
	{
		for (col = 1; col < gis_grid.ncols - 1; col++)
		{
			gis_grid.xslope[row][col] = ( ( gis_grid.elev[row-1][col+1] - gis_grid.elev[row-1][col-1] ) +
					2 * ( gis_grid.elev[row  ][col+1] - gis_grid.elev[row  ][col-1] ) +
						( gis_grid.elev[row+1][col+1] - gis_grid.elev[row+1][col-1] ) ) /
						( 8 * gis_grid.resolution );
			gis_grid.yslope[row][col] = ( ( gis_grid.elev[row-1][col-1] - gis_grid.elev[row+1][col-1] ) +
					2 * ( gis_grid.elev[row-1][col  ] - gis_grid.elev[row+1][col  ] ) +
						( gis_grid.elev[row-1][col+1] - gis_grid.elev[row+1][col+1] ) ) /
						( 8 * gis_grid.resolution );
			gis_grid.slope[row][col] = sqrt ( gis_grid.xslope[row][col]*gis_grid.xslope[row][col] +
												gis_grid.yslope[row][col]*gis_grid.yslope[row][col] );
		}
	}
	for (col = 1; col < gis_grid.ncols - 1; col++)
	{
		gis_grid.xslope[0][col] = gis_grid.xslope[1][col];
		gis_grid.xslope[gis_grid.nrows-1][col] = gis_grid.xslope[gis_grid.nrows-2][col];
		gis_grid.yslope[0][col] = gis_grid.yslope[1][col];
		gis_grid.yslope[gis_grid.nrows-1][col] = gis_grid.yslope[gis_grid.nrows-2][col];
		gis_grid.slope[0][col] = gis_grid.slope[1][col];
		gis_grid.slope[gis_grid.nrows-1][col] = gis_grid.slope[gis_grid.nrows-2][col];
	}
	for (row = 0; row < gis_grid.nrows; row++)
	{
		gis_grid.xslope[row][0] = gis_grid.xslope[row][1];
		gis_grid.xslope[row][gis_grid.ncols-1] = gis_grid.xslope[row][gis_grid.nrows-2];
		gis_grid.yslope[row][0] = gis_grid.yslope[row][1];
		gis_grid.yslope[row][gis_grid.ncols-1] = gis_grid.yslope[row][gis_grid.nrows-2];
		gis_grid.slope[row][0] = gis_grid.slope[row][1];
		gis_grid.slope[row][gis_grid.ncols-1] = gis_grid.slope[row][gis_grid.nrows-2];
	}
	return 0;
}

int calculate_curvature()
{
	int row, col;
	int status;

	if ( gis_grid.xslope == 0 )
	{
		status = calculate_slope();
		if ( status != 0 )
			return status;
	}

	free_float_matrix ( gis_grid.xcurv );
	free_float_matrix ( gis_grid.ycurv );

	if ( ! ( gis_grid.xcurv = alloc_float_matrix ( gis_grid.nrows, gis_grid.ncols ) ) )
		return -3;	/*memory error*/

	if ( ! ( gis_grid.ycurv = alloc_float_matrix ( gis_grid.nrows, gis_grid.ncols ) ) )
		return -3;	/*memory error*/

	for (row = 1; row < gis_grid.nrows - 1; row++)
	{
		for (col = 1; col < gis_grid.ncols - 1; col++)
		{
			gis_grid.xcurv[row][col] = ( ( gis_grid.slope[row-1][col+1] - gis_grid.slope[row-1][col-1] ) +
					2 * ( gis_grid.slope[row  ][col+1] - gis_grid.slope[row  ][col-1] ) +
						( gis_grid.slope[row+1][col+1] - gis_grid.slope[row+1][col-1] ) ) /
						( 8 * gis_grid.resolution );
			gis_grid.ycurv[row][col] = ( ( gis_grid.slope[row-1][col-1] - gis_grid.slope[row+1][col-1] ) +
					2 * ( gis_grid.slope[row-1][col  ] - gis_grid.slope[row+1][col  ] ) +
						( gis_grid.slope[row-1][col+1] - gis_grid.slope[row+1][col+1] ) ) /
						( 8 * gis_grid.resolution );
		}
	}
	for (col = 1; col < gis_grid.ncols - 1; col++)
	{
		gis_grid.xcurv[0][col] = gis_grid.xcurv[1][col];
		gis_grid.xcurv[gis_grid.nrows-1][col] = gis_grid.xcurv[gis_grid.nrows-2][col];
		gis_grid.ycurv[0][col] = gis_grid.ycurv[1][col];
		gis_grid.ycurv[gis_grid.nrows-1][col] = gis_grid.ycurv[gis_grid.nrows-2][col];
	}
	for (row = 0; row < gis_grid.nrows; row++)
	{
		gis_grid.xcurv[row][0] = gis_grid.xcurv[row][1];
		gis_grid.xcurv[row][gis_grid.ncols-1] = gis_grid.xcurv[row][gis_grid.nrows-2];
		gis_grid.ycurv[row][0] = gis_grid.ycurv[row][1];
		gis_grid.ycurv[row][gis_grid.ncols-1] = gis_grid.ycurv[row][gis_grid.nrows-2];
	}
	return 0;
}

void find_min_max()
{
	int nrows, ncols;
	int row, col;
	nrows = gis_grid.nrows;
	ncols = gis_grid.ncols;

	for (row = 0; row < nrows; row++)
	{
		for (col = 0; col < ncols; col++)
		{
			if ( gis_grid.elev[row][col] > gis_grid.zmax )
				gis_grid.zmax = gis_grid.elev[row][col];
			if ( gis_grid.elev[row][col] < gis_grid.zmin )
				gis_grid.zmin = gis_grid.elev[row][col];
		}
	}
}

int print_grid()
{
	int row, col;

	if ( ! gis_grid.elev )
		return -4;

	for (row = 0; row < gis_grid.nrows; row++)
	{
		fprintf (stdout, "row %d\n", row );
		for (col = 0; col < gis_grid.ncols; col++)
			fprintf (stdout, "%12.8f ", gis_grid.elev[row][col] );
		fprintf (stdout, "\n" );
	}
	fprintf (stdout, "\n" );
	fprintf (stdout, "x=%12.8f y=%12.8f \n", gis_grid.wxmin, gis_grid.wymin);

	return 0;
}

int set_from_header(GisRasterHdr& gisHeader)
{
//	Copy the resolutions
	double xres = gisHeader.XRes();
	double yres = gisHeader.YRes();
	if ( fabs(xres- yres)>=xres/10000.0 )
		return -4;	// resolution must be the same
	gis_grid.wresolution = gis_grid.resolution = xres;

//	Copy the edges of the region
	//	Set rows and cols
	gis_grid.nrows = gis_grid.wnrows = gisHeader.Rows();
	gis_grid.ncols = gis_grid.wncols = gisHeader.Cols();

	gis_grid.wxmin = gis_grid.xmin = gisHeader.West();
//	Use this instead of
//	gis_grid.wxmax = gis_grid.xmax = gisHeader.east
//	to ensure correct ncols
	gis_grid.wxmax = gis_grid.xmax = gisHeader.West() + gisHeader.Cols() * xres;
	gis_grid.wymin = gis_grid.ymin = gisHeader.South();
//	Use this instead of
//	gis_grid.wymax = gis_grid.ymax = gisHeader.north
//	to ensure correct nrows
	gis_grid.wymax = gis_grid.ymax = gisHeader.South() + gisHeader.Rows() * yres;
	
	return 0;
}
