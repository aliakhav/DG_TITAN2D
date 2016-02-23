#ifndef ELEMENT_H
#define ELEMENT_H
#include "node.h"
#include "boundary.h"

class Element {


  friend void drawit(Element);

 public:

  Element();
  void setparameters(int, Node*[], int, int*);
  void order_nodes();
  void case1();
  void case2();
  void case3();
  void case4();
  void case5();
  void create_m_node(double*, double*);
  int getid(){return elementid;};
  //void determine_the_key(unsigned, double*, double*);
  unsigned* pass_key(); 
  void determine_neighbors(int, Element*);
  Node** get_element_node(){return element_nodes;};
  void myproc(int, int, int);
  void write_node_data(ofstream*);
  void write_element_data(ofstream*);
  void reset_written_flag();
  void find_boundary(int, int, Boundary*);
  Edge* get_element_edges();
  Element* get_neighbors(int);
  void determine_opposite_brother();

 private:
  int elementid;
  Boundary* boundary[4][2];
  Node* element_nodes[9];
  Edge* element_edges[4];
  Element* neighbor[4];
  int material;
  int myprocess; 
  int elm_loc[2];
  int which_son;
  Element* opposite_brother;

};
#endif
