// node.h for the definition of node class
#ifndef EDGE_H
#define EDGE_H
#include "node.h"

class Edge {

 public:

  int getid();
  int* getowners();
  Node* get_edge_nodes();

 private:

  int edgeid;
  Node* edge_nodes[3];
  int owner[2];


};
#endif
