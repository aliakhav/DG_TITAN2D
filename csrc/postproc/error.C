#include "../header/hpfem.h"

/*--  u = 3/2 r^(2/3) sin (2theta/3) ---*/
void error(HashTable* ht_elem_ptr, HashTable* ht_node_ptr, Element* emtemp,
	   double* errsq, double* solsq)
{
  int          i, j, k;
  Node*        ndtemp;
  Node*        ndtemp1;
  Node*        ndtemp2;
  unsigned*    keyp;
  void*        p;
  Element*     father;
  double       exact, x, y, r, theta, sol, err, ramp;
  *solsq = 0.0; *errsq = 0.0;
  err = 0.0; sol = 0.0;
  
  for(i=0;i<4;i++)
    {
      ndtemp = (Node*)(ht_node_ptr->lookup(emtemp->pass_nodes()+i*KEYLENGTH));
      x = *ndtemp->get_coord();  
      y = *(ndtemp->get_coord()+1);
      r = sqrt(x*x+y*y);
      if(r<1.0e-10) ramp = .0;
      else ramp = y/r;
      theta = asin(ramp);
      if(x<0.0) theta = PI-theta;
      exact = (3.0/2.0)*pow(r, 0.6667)*sin(0.6667*theta);
      if(ndtemp->getinfo()==CORNER) 
	{
	  sol = *ndtemp->getsol();
	  if(sol<1.0e-3) err = 0.0;
	  else err = (sol - exact)/sol;
	  *errsq = *errsq+err*err;
	  *solsq = *solsq+sol*sol;
	}
      else if(ndtemp->getinfo()==S_C_CON)
	*solsq = *solsq+exact*exact;
    }
  *errsq = *errsq/4.0;
  *solsq = *solsq/4.0;
  //cout<<*errsq<<"  "<<*solsq<<'\n'<<flush;

}
