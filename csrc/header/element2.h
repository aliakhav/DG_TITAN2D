#ifndef ELEMENT_H
#define ELEMENT_H
#include <math.h>
#include "hashtab.h"
#include "geoflow.h"
#include "node.h"
#include "struct.h"
#include <fstream>
#include <iostream>
#include <assert.h>

class Element {

  friend class HashTable;

  friend void BSFC_combine_elements(int, Element* EmTemp, HashTable*, 
				    HashTable*, int);
    
  friend void Pack_element(Element* sendel, ElemPack** elemptr, 
			   HashTable* HT_Node_Ptr, int);
  friend void Pack_element(Element* sendel, ElemPack* elem, HashTable* HT_Node_Ptr, int);

  friend void destroy_element(Element* r_element, HashTable* HT_Elem_Ptr, 
			      HashTable* HT_Node_Ptr, int target_pro, ELinkPtr* EL_head);
  
  friend void create_element(ElemPack* elem2, HashTable* HT_Elem_Ptr, 
			     HashTable* HT_Node_Ptr, int myid, double* e_error);

  friend void construct_el(Element* newelement, ElemPack* elem2,
			   HashTable* HT_Node_Ptr, int myid);

  friend void unrefine_elements(Element**, HashTable*, HashTable*);
  
 
 public:

  Element() {};  

  Element(unsigned[][KEYLENGTH], unsigned[][KEYLENGTH], 
	  int[], int mat, double pile_height[], int myid, int*, unsigned*);//original element

  /*for the refinement:::*/
  Element(unsigned nodekeys[][2], unsigned neigh[][2], 
	  int n_pro[], int gen, int ord, int gen_neigh[], int mat, 
	  Element* fthTemp, HashTable* El_Table, 
	  HashTable* NodeTable, int myid, MatProps* matprops_ptr, int elm_loc_in[]);
  /* for the unrefinement */
  Element(Element* sons[], HashTable* NodeTable, HashTable* El_Table, MatProps* matprops_ptr);

  ~Element();
  
  unsigned* pass_key(){return key;};
  int       get_material(){return material;};
  void      get_stiffness(HashTable*, HashTable*, double*, double*, Element*);
  void      project_sol(HashTable* NodeTable, HashTable* HT_Elem_Ptr, int parent_son_flag);
  unsigned* getNode();
  int*      getassoc();
  int       get_no_of_dof();
  void      put_gen(int);
  void      put_fluxIntegral(double);
  void      putson(unsigned*);
  void      putbrothers(unsigned*);
  unsigned* get_brothers() {return &brothers[0][0];};
  void      putfather(unsigned*);
  void      putassoc(int, int);
  void      putneighbor(unsigned*, int);
  void      put_neigh_proc(int, int);
  void      put_order(int);
  int       get_order();     
  unsigned* getfather();
  unsigned* getson();
  void      putel_sq(double, double);
  double*   get_el_solution();
  double*   get_prev_el_solution();
  double*   get_temp_el_solution();
  double*   get_el_error();
  unsigned* get_neighbors();
  int*      get_neigh_proc();
  int       get_gen();
  double    get_fluxIntegral();
  int       which_neighbor(unsigned*);  //used to find out which neighbor for a given key
  int       get_refined_flag();
  void      change_neighbor(unsigned*, int, int, int);
  void      put_refined_flag(int);
  int*      get_neigh_gen();
  void      put_neigh_gen(int, int);
  void      put_which_son(int);
  int       get_which_son();
  void      put_new_old(int);
  int       get_new_old();
  void      update_ndof();
  int       get_ndof();
  void      change_neighbor_process(int, int);
  double*   get_el_err();
  void      get_nelb_icon(HashTable*, HashTable*,int*,int*);
  double    get_lb_weight() {return lb_weight;};
  void      put_lb_weight(double dd_in) {lb_weight = dd_in;};
  unsigned* get_lb_key() {return lb_key;};
  void      put_lb_key(unsigned* in_key);
  void      copy_key_to_lb_key();
  void      put_myprocess(int in_proc) {myprocess = in_proc;};
  int       get_myprocess() {return myprocess;};
  int       BSFC_check_send_elements(); // send all of the parents
  int       BSFC_check_send_elements(HashTable*); // to make sure a father gets sent if an elm has a constrained node
  // abani adds
  void      store_el_solution(double *,int); //
  void      store_prev_el_solution(double *,int); //
  void      store_temp_el_solution(double *, int);
  void      store_temp_prev_solution(double *, int);
  void      update_el_solution(); // move solution values to previous and update with new solution
  void      store_prev_el_solution();// move prev solution values to new solution
  void      store_temp_prev_solution();// move prev solution values to new solution 
  void      slopelimit(HashTable*, HashTable*, MatProps* );

  //void      get_icon(HashTable*, HashTable*, int[4]);
  //void      get_boundary(int[4], double[4]);
  /* geoflow functions */
  void      correct(HashTable* NodeTable, HashTable* El_Table,
		    double dt, MatProps*);
  int       determine_refinement(double);
  double    calc_volume(HashTable*);
  double    calc_time_step(HashTable*, MatProps*);
  
  //unrefinement
  int find_brothers(HashTable* El_Table, HashTable*, double, int, MatProps*);
  int check_unrefinement(double);
  void change_neigh_info(unsigned*, unsigned*, int);
  void check_fathers_neighbors(HashTable* El_Table, HashEntryPtr* parentlistptr, int myid);
  void      find_opposite_brother(HashTable*);
  int       get_opposite_brother_flag() {return opposite_brother_flag;};
  // unrefinement beyond the coarse grid
  int* get_elm_loc() {return elm_loc;};
  void put_elm_loc(int* int_in) {elm_loc[0] = int_in[0]; elm_loc[1] = int_in[1];};
  void      calc_which_son();

 private:  
  int myprocess;
  int generation;
  double fluxIntegral;
  int material;/*! ! ! THE MAT. FLAG ! ! !*/
  double lb_weight; //load-balancing weight
  unsigned lb_key[KEYLENGTH];  //key for load-balancing (if there is no constrained node, it is the element key, otherwise it is a construct of the element "bunch"
  unsigned key[KEYLENGTH];
  unsigned node_key[8][KEYLENGTH];
  unsigned neighbor[8][KEYLENGTH];//--considering 1-irr
  unsigned son[4][KEYLENGTH];//--garantee ccw
  int neigh_proc[8];//--considering 1-irr, neigh_proc[4:7] != -2 only if it has 2 neighbors on that side, -1 means bc
  int order;
  int neigh_gen[8];
  int ndof;
  // el_solution is stored as linear h, hu, hv at node 0, linear h, hu, hv at node 1, ... 
  // quadratic h, hu, hv at node 0, quadratic h, hu, hv at node 1,..., quadratic h, hu, hv at bubble node
  double el_error[EQUATIONS];  
  double el_solution[ELM_DOF];
  double prev_el_solution[ELM_DOF];
  double temp_el_solution[ELM_DOF];
  int    refined;
  int    which_son;
  int    new_old;
  unsigned brothers[4][KEYLENGTH];
  int elm_loc[2];
  int opposite_brother_flag;  //flag to indicate if we have the correct key for the non-neighbor brother (0:= don't have info, 1:= have info)
};

inline int* Element:: getassoc(){ return neigh_proc;}

inline unsigned* Element:: getNode(){ return &(node_key[0][0]);}

inline int Element::get_no_of_dof(){return ndof;}

inline void Element::put_gen(int g){generation = g;}

inline void Element::put_fluxIntegral(double flux){fluxIntegral = flux;}

inline void Element::put_neigh_proc(int i, int proc) { neigh_proc[i] = proc;}

inline void Element::put_order(int ord) {order = ord;}

inline unsigned* Element::getson(){return &(son[0][0]);}

inline int Element::get_order(){return order;}

inline void Element::putson(unsigned* s)
{

  for(int i=0; i<4; i++)
    for(int j=0; j<KEYLENGTH;j++)  
      son[i][j] = *(s+i*KEYLENGTH+j);

  refined=1;
}

inline void Element::putbrothers(unsigned* s)
{

  for(int i=0; i<4; i++)
    for(int j=0; j<KEYLENGTH;j++)  
      brothers[i][j] = *(s+i*KEYLENGTH+j);

  return;
}

inline void Element::putneighbor(unsigned* n, int i)
{
  int j;
  for(j=0; j<KEYLENGTH;j++)  neighbor[i][j] = *(n+j);
}

inline void Element::putassoc(int a, int i){neigh_proc[i] = a;}

inline double* Element::get_el_solution()
{return el_solution;}

inline double* Element::get_prev_el_solution()
{return prev_el_solution;}

inline double* Element::get_temp_el_solution()
{return temp_el_solution;}

inline double* Element::get_el_error()
{return el_error;}

inline unsigned* Element::get_neighbors()
{return &neighbor[0][0];}

inline int* Element::get_neigh_proc()
{return neigh_proc;}

inline int Element::get_gen()
{return generation;}

inline double Element::get_fluxIntegral()
{return fluxIntegral;}

inline int Element::get_refined_flag()
{return refined;} 

inline void Element::put_refined_flag(int i)
{refined = i;}

inline int* Element::get_neigh_gen() {return neigh_gen;}

inline void Element::put_neigh_gen(int i, int gen) {neigh_gen[i] = gen;}

inline void Element::put_which_son(int i){ which_son = i;}

inline int  Element::get_which_son(){return which_son;}

inline int  Element::get_new_old(){return new_old;}

inline void Element::put_new_old(int i) {new_old = i;}

inline void Element::change_neighbor_process(int which, int newp){neigh_proc[which]=newp;}

//abani adds
inline void Element::store_el_solution(double * newsol,int Nc)
{
  int i;
  for(i=0;i<Nc;i++){
    el_solution[i]=newsol[i];
    if (el_solution[0]<=0){
      el_solution[i]=0;
    }
  }  

return;
}

inline void Element::store_prev_el_solution()
{
  int i;
  for(i=0;i<ndof;i++)
    *(get_el_solution()+i)=*(get_prev_el_solution()+i);
  
  
  for(i=ndof;i<ELM_DOF;i++)
    el_solution[i]=0;
  
  
  for(i=0;i<ELM_DOF;i++)
    if (el_solution[0] <= 0)
      el_solution[i]=0;

  for(i=0;i<ELM_DOF;i++)
    {
      //  printf("prev_el_sol[%d]=%e\n",i,prev_el_solution[i]);
      assert(prev_el_solution[i] >= -1.E30);
    }
  
  return;
}

inline void Element::store_temp_prev_solution()
{
  int i;
  for(i=0;i<ndof;i++)
    prev_el_solution[i]=temp_el_solution[i];
  
  for(i=ndof;i<ELM_DOF;i++)
    prev_el_solution[i]=0;  
  
  for(i=0;i<ELM_DOF;i++)
    if (prev_el_solution[0]<=0)
      prev_el_solution[i]=0;
  
  
 for(i=0;i<ELM_DOF;i++)
    {
      // printf("prev_el_sol[%d]=%e\n",i,prev_el_solution[i]);
      assert(prev_el_solution[i] >=-1.E30);
    }
 
  return;
}


inline void Element::store_temp_el_solution(double * newsol,int Nc)
{
  int i;
  for(i=0;i<Nc;i++){
    temp_el_solution[i]=newsol[i];
  }

  for(i=Nc;i<ELM_DOF;i++){
    temp_el_solution[i]=0.;
  }

  for(i=0;i<ELM_DOF;i++)
    if (temp_el_solution[0]<=0.)
      temp_el_solution[i]=0.;
    
  
  return;
}

inline void Element::store_prev_el_solution(double * newsol,int Nc)
{
  int i;
  for(i=0;i<Nc;i++)
    prev_el_solution[i]=newsol[i];
  return;
}
inline void Element::update_el_solution()
{
  int i;
  for(i=0;i<ndof;i++)
    prev_el_solution[i]=el_solution[i];

  for(i=ndof;i<ELM_DOF;i++)
    prev_el_solution[i]=0;
  

  for(i=0;i<ELM_DOF;i++)
    if (el_solution[0]<=0){
      prev_el_solution[i]=0;
      el_solution[i]=0;
    }
  
  //if (el_solution[0] <= 0.5*GEOFLOW_TINY && i!=0){
  //  prev_el_solution[i]=0;
  //  el_solution[i]=0;
  //  }

  return;
}

inline int Element::get_ndof() {
  return ndof;
}


#endif


