#ifndef STRUCT_H
#define STRUCT_H


struct ElemPack{
  int        myprocess;
  int        generation;
  int        material;/*flag added by andrew*/
  int        neigh_proc[8];
  int        order;
  int        neigh_gen[8];
  int        ndof;
  int        refined;
  int        which_son;
  int        new_old;
  int        n_info[9];
  int        elm_loc[2];
  int        opposite_brother_flag;

  unsigned   key[KEYLENGTH];/*contains the 9th node key*/
  unsigned   node_key[8][KEYLENGTH];
  unsigned   neighbor[8][KEYLENGTH];
  unsigned   son[4][KEYLENGTH]; 
  unsigned   brothers[4][KEYLENGTH];

  double node_elevation[9]; 
  double n_coord[9][2];
  double el_error[EQUATIONS];
  double el_solution[ELM_DOF];
  double prev_el_solution[ELM_DOF];
  double lb_weight;

};

//                             \|||/    
//                             (o o)   
//---------Elementlink------oo0-(_)-0oo----------STARTS HERE-------


struct ElementLink{
  int             target_proc;
  int             new_proc;
  unsigned        elkey[KEYLENGTH]; 
  unsigned        targetkey[KEYLENGTH];

  ElementLink*    pre;
  ElementLink*    next;
  
  ElementLink(unsigned* keyi,unsigned* key2, int tp, int np)
    {
      target_proc = tp;
      new_proc    = np;
      int i;
      for(i=0;i<KEYLENGTH; i++)
	{
	  elkey[i] = *(keyi+i);
	  targetkey[i]= *(key2+i);
	}
      next   = NULL;
      pre    = NULL;
    }
  ElementLink(){
    next   = NULL;
    pre    = NULL;
  }
  ~ElementLink(){
    if(next)
      next->pre = pre;
    if(pre)
      pre->next = next;
  }


};

typedef ElementLink* ELinkPtr;


//                              \|||/    
//                              (o o)   
//--------MPI datatype-------oo0-(_)-0oo----------STARTS HERE-------


struct NeighborPack{
 
  int             target_proc;
  int             new_proc;
  unsigned        elkey[KEYLENGTH]; 
  unsigned        targetkey[KEYLENGTH];


};

struct Neigh_Sol_Pack{


int      nside;
int      norder[5];
unsigned key[KEYLENGTH];
double   solu[2][121];
double   Xnod[18];

};

typedef NeighborPack* NePtr;

#endif
