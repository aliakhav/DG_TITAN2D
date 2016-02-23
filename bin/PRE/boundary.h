#ifndef BOUNDARY_H
#define BOUNDARY_H
#include "node.h"


class Boundary {

friend class Element;

 public:

  Boundary();
  void setparameters(Node*, double, double, int);
  Node* get_boundary_node(){return node;};
  int get_type() {return type;};
  double get_x_value(){return x_value;};
  double get_y_value(){return y_value;};
  void write_b_data(ofstream*);

 private:
  Node* node;
  int type;//-2 essential -3 natural
  double x_value;
  double y_value;

};
#endif
