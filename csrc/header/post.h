#ifndef P_H
#define P_H


class Neigh_Sol{

public:
  unsigned*  key;
  double     solu[2][121];
  int        nside;
  int*       norder;
  double     Xnod[18];
  Neigh_Sol* next;

  Neigh_Sol(unsigned* in_key)
    {
      key = in_key;
      nside =-1;
      
      for(int i=0;i<2;i++)
	for(int j=0;j<121;j++)
	  solu[i][j] =0.0;
      
      next=NULL;
    }
  
  
  Neigh_Sol()
    {
      key = NULL;
      nside =-1;
      
      for(int i=0;i<2;i++)
	for(int j=0;j<121;j++)
	  solu[i][j] =0.0;
      
      next = NULL;
    }

};//end of class


#endif
