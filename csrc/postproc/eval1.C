//----------------------------------------------------------------------
//   routine name       - eval1
//   latest revision    - May, 1993
//   purpose            - routine evaluates value of a 1-D solution
//                        at a point with master coodinate Eta
//   arguments :
//     in:     Norder   - order of approximation
//             Eta      - the master coordinate 
//             Nvalues  - node values(solution) in hierarchical shape function
//     out:    Val      - values of the solution
//----------------------------------------------------------------------

#include "../header/hpfem.h"

#ifdef SUNOS
extern "C" void shape1_(int* Norder, double* Eta, double* vshap);
#endif

#ifdef IBMSP
extern "C" void shape1(int* Norder, double* Eta, double* vshap);
#endif

#ifdef CRAY
extern "C" void SHAPE(int* Norder, double* Eta, double* vshap);
#endif


void eval1(int Norder, double Eta, double Nvalues[2][11], double Val[2])
{      
 
  double  vshap[11];

  //  ...explanation of local variables:
  //     vshap   - 1-D hierarchical shape functions
  //  ...evaluate shape functions at the point
#ifdef SUNOS
  shape1_(&Norder, &Eta, vshap);
#endif

#ifdef IBMSP
  shape1(&Norder, &Eta, vshap);
#endif

#ifdef CRAY
  SHAPE(&Norder, &Eta, vshap);
#endif
  
  //  ...accumulate for the final values
  for ( int ivar=0; ivar < EQUATIONS; ivar++ ) 
    {
      double s = 0.0;

      for ( int k = 0; k < Norder+1; k++ ) //need to add values from linears at both corners
	s = s + Nvalues[ivar][k]*vshap[k];

      Val[ivar] = s;
    }

}
