#include "../header/hpfem.h"

/*...calculates the solution at a given master element coordinate...*/
#ifdef SUNOS
extern "C" void shape2_(int* order, double* X, double* Y, double* vshape);
#endif

#ifdef IBMSP
extern "C" void shape2(int* order, double* X, double* Y, double* vshape);
#endif

#ifdef CRAY
extern "C" void SHAPE2(int* order, double* X, double* Y, double* vshape);
#endif

void calculate_solution(double coeff[], int order[], int ndof, double X, double Y, double solution[])
{


  double* vshape;
  
  int i, j;

  //double vshape[9];
  
  int functions=ndof/EQUATIONS;/*...the number of shape functions...*/
  //cout<<"the size now is: "<<size<<endl<<flush;
  vshape=new double[functions];
  
  for(i=0; i<functions; i++)
    vshape[i]=0;


  solution[0]=0;
  solution[1]=0;

  /*...call shape2 to get the values of shape functions...*/
#ifdef SUNOS
  shape2_(order, &X, &Y, vshape);
#endif

#ifdef IBMSP
  shape2(order, &X, &Y, vshape);
#endif

#ifdef CRAY
   SHAPE2(order, &X, &Y, vshape);
#endif
  

  for(i=0; i<EQUATIONS; i++)
    for(j=0; j<functions; j++)
      {
	solution[i]+=coeff[j*EQUATIONS+i] * vshape[j];
      }


  //delete[] vshape;


}
