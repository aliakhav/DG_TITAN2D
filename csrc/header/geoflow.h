#ifndef __GEOFLOW
#define __GEOFLOW

/* geoflow header file */
#define GEOFLOW_TINY 0.0001
#define WEIGHT_ADJUSTER 1
#define REFINE_LEVEL 3 
#define REFINE_THRESHOLD 0.006  //for refinement ONLY--if current elemnt height is less than 
                               //this threshol and my neighbor is greater than this 
                               //then refine current element

void checknodesol(HashTable*);

double element_weight(HashTable* El_Table, HashTable*, int myid, int nump);

void calc_volume(HashTable* El_Table, HashTable* NodeTable, 
		 int myid, int nump, MatProps* matprops_ptr);

void setup_geoflow(HashTable* El_Table, HashTable* NodeTable, int myid, int nump, MatProps*);

double get_coef_and_eigen(HashTable* El_Table, HashTable*, MatProps*, int);

void move_data(int nump, int myid, HashTable* El_Table, HashTable* NodeTable);

void delete_ghost_elms(HashTable* El_Table, int myid) ;

void calc_edge_states(HashTable* El_Table, HashTable* NodeTable,
		      int myid);

/* c++ sgn function */

inline double c_sgn(double zz)
{
  double sgn;
  if(zz > GEOFLOW_TINY)
    sgn = 1;
  else if(zz < -GEOFLOW_TINY)
    sgn = -1;
  else 
    sgn = 0;

  return sgn;
}

/* c++ dmin1 functions */
inline double c_dmin1(double d1, double d2) {
  
  if(d1> d2)
    d1 = d2;

  return d1;
}

inline double c_dmin1(double d1, double d2, double d3) {
  
  if(d1> d2)
    d1 = d2;
  if(d1 > d3)
    d1 = d3;

  return d1;
}

/* c++ dmax1 functions */
inline double c_dmax1(double d1, double d2) {
  
  if(d1 < d2)
    d1 = d2;

  return d1;
}

inline double c_dmax1(double d1, double d2, double d3) {
  
  if(d1 < d2)
    d1 = d2;
  if(d1 < d3)
    d1 = d3;

  return d1;
}

/* c++ dabs function */
inline double dabs(double dd) {
  if(dd < 0)
    dd = -dd;

  return dd;
}
  
/* fortran calls */
#ifdef SUNOS 
extern "C" void gmfggetcoef_(double*, double*, double*, double*, double*,
			     double*, double*, double*, double*, double*);
extern "C" void eigen_(double*, double*, double*, double*, double*, double*, double*);
extern "C" void predict_(double*, double*, double*, double*,
			 double*, double*, double*, double*, double*,
			 double*, double*, double*, double*);
extern "C" void correct_(double*, double*, double*, double*, double*,
			 double*, double*, double*, double*, double*,
			 double*, double*, double*, double*, double*, 
			 double*, double*, double*, double*);
extern "C" void shape2dg_(int*,double*,double*,double*);
extern "C" void dshap2dg_(int*,double*,double*,double*);
extern "C" void shape2_(int*,double*,double*,double*);
extern "C" void gshape_(double*,double*,double*);
extern "C" void getkactxy_(int*, int*, double*, double*, double*, double*, double*, double*, double*, double*);

#endif
#ifdef IBMSP
extern "C" void gmfggetcoef(double*, double*, double*, double*,
				 double*, double*, double*, double*, 
				 double*, double*, double*);
extern "C" void eigen(double*, double*, double*, double*, double*, double*);

#endif


#endif
