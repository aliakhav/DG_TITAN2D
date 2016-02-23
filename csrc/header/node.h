#ifndef NODE_H
#define NODE_H

#include "properties.h"
#include "constant.h"
#include "hashtab.h"
#include "struct.h"

class Element;
class Node { 
  
  friend class Element;

  friend void Pack_element(Element* sendel, ElemPack** elemptr, HashTable* HT_Node_Ptr, int);

  friend void Pack_element(Element* sendel, ElemPack* elem, HashTable* HT_Node_Ptr, int);

  friend void destroy_element(Element* r_element, 
		    HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr);

  friend void create_element(ElemPack* elem2, HashTable* HT_Elem_Ptr, 
		    HashTable* HT_Node_Ptr, double* e_error);
  
  friend void unrefine_elements(Element**, HashTable*, HashTable*);

 public:
  Node(unsigned* keyi, double*, MatProps*);

  Node(unsigned*, double*, int, MatProps*);/*for refined*/

  Node(unsigned* keyi, double* coordi, int inf, 
       double elev);

  
  Node();

  ~Node();
  
  void putdof(int,int);
  int* getdof();
  void putglnum(int);
  int  getglnum();
  int  getinfo() {return info;}
  void putinfo(int in);
  unsigned* pass_key();
  double* get_coord() {return coord;};
  int  get_reconstructed();
  void put_reconstructed(int);
  void put_id(int id_in) {id = id_in;};
  int get_id() {return id;};
  /*  geoflow methods */
  double get_elevation() {return elevation;};
  
  
 protected:
  int       id;/*--used in delete_unused_nodes_and_elements() function --*/
  int       info;
  double    coord[DIMENSION];
  unsigned  key[KEYLENGTH];
  void*     nextptr;
  void*     preptr;
  int       dof[2];/*--dof[1]-dof[0]+1 = dof of the node---*/
  int       glnum;/*--the node occupies the position from glnum to glnum+dof--*/
  int       reconstructed;
  double    elevation; // this elevation should currently be the GIS elevation at the finest "scale"

};


inline void Node:: putdof(int lower, int up){
     dof[0] = lower;
     dof[1] = up;
}

inline int* Node:: getdof(){
     return dof;
}

inline void Node:: putglnum(int numbering){
     glnum =  numbering;
}

inline int Node:: getglnum(){
     return glnum;
}


inline unsigned* Node:: pass_key(){
     return key;
}

inline void Node::put_reconstructed(int i)
{
  reconstructed = i;
}

inline int Node::get_reconstructed()
{
  return reconstructed;
}

#endif








































