#define  SUNOS  //definition for gmake architecture
//define HAVE_LIBHDF5
#define TOPO_DATA

#ifdef CRAY
#include <fortran.h>
#endif 

#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <fstream>
#include "blas.h"
#include "properties.h"
#include "node.h"
#include "element2.h"
#include "interface.h"
#include "hashtab.h"
#include <assert.h>
#include "recv.h"
#include "extfun.h"
//#include "glvar.h"
#include "post.h"
#include "geoflow.h"
#include "scale.h"
#include "GisApi.h"



#undef CRAY

