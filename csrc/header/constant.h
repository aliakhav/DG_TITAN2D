#ifndef CONSTANTS
#define CONSTANTS
const int  KEYLENGTH = 2;

//const   NODE_TABLE_SIZE=6000;
//const   EL_TABLE_SIZE=2000;


const int  MAX_PROCS = 2056;

const int  ZERO      = '0';
const int  NONZERO   = 111;
//const   NBUCKETS  = 2000;        //every hashtable posseses 5000 entries.
const int   PRIME     = 2017;        //used for creating hash key
const int   DIMENSION = 2;
const int   EQUATIONS = 3;
const int   ORDER     = 2;
const int   MAX_ORDER = 3;  // 6 shape functions when order is 3
const int  ELM_DOF = EQUATIONS*9;  // the number of dof associated with an element

const double PI   = 3.1415926;

const int   BCTYPE    = 3;

const float  C = 1.5;

const int   NODEINIT    = 0x0000;
const int   CORNER      = 0x0002;//--corner dof
const int   BUBBLE      = 0x0006;
const int   SIDE        = 0x0004;//--side dof
const int   CONSTRAINED = 0x0001;
const int   S_C_CON     = 0x0007;//--side dof
const int   S_S_CON     = 0x0005;//--no dof
const int   ASSIGNED    = 1;
const int   UNASSIGNED  = 0;

const double  UN_CONSTRAINED =-999999.0;/*for the unconstrained dof*/

const int   ON = 1;
const int   OFF = 0;

const int   NEW = 1;
const int   OLD = 0;

const int   INIT = -1;

const float  MAX_X = 1.0;
const float  MAX_Y = 1.0;
const float  MIN_X = -1.0;
const float  MIN_Y = -1.0;

const float  SMALL = 1.0e-5;
const float  XBC = -1.0;
const float  YBC = -1.0;
const float  LOAD_BALANCE_TOLERANCE = 1.0001;

extern void fhsfc3d(double*, unsigned*, unsigned*);
extern void fhsfc2d_(double*, unsigned*, unsigned*);

/* geoflow data */
const int GHOST = -9876;  // indicates a "ghost" element for data storage

#endif
