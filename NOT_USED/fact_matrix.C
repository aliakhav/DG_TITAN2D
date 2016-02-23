#include "../header/hpfem.h"
#include "../header/blas.h"

void fact_matrix(double stiff[], double load[],int dof, int int_dof)
{
/* stiff matrix is:
| EE EB | |u_e|   |f_e|  EE, EB, BE and BB are in stiff
| BE BB | |u_b| = |f_b|  dof is dimension of stiff, int_dof is interior unknowns

*/

  int i=0, j=0, k=0;  
  int ext_dof = dof-int_dof;
  //  mat_zero_(stiff, &dof, &dof);
  //mat_view_(stiff, &dof, &dof);

  
  // schur complement for BB matrix
  double* rhs2 = new double[dof];
  double* u2 = new double[dof];

#ifdef SUNOS
  schur5_(stiff, (stiff+ext_dof*dof), (stiff+ext_dof), (stiff+(dof+1)*ext_dof),
	  &ext_dof, &dof, &int_dof, &dof, load, (load+ext_dof), u2, rhs2);
#endif
#ifdef IBMSP
  schur5(stiff, (stiff+ext_dof*dof), (stiff+ext_dof), (stiff+(dof+1)*ext_dof),
	  &ext_dof, &dof, &int_dof, &dof, load, (load+ext_dof), u2, rhs2);
#endif
#ifdef CRAY
  SCHUR5(stiff, (stiff+ext_dof*dof), (stiff+ext_dof), (stiff+(dof+1)*ext_dof),
	  &ext_dof, &dof, &int_dof, &dof, load, (load+ext_dof), u2, rhs2);
#endif
  delete []rhs2;
  delete []u2;
  
  return;
}
