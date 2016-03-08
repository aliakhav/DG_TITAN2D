#include <iostream>
#include <fstream>
#include <stdio.h>
#include "boundary.h"


Boundary::Boundary(){}


void Boundary::setparameters(Node* n, double xv, double yv, int t)
{

  node=n;
  x_value=xv;
  y_value=yv;
  type=t;

}

void Boundary::write_b_data(ofstream* out)

{
  *out<<type<<' '<<node->key[0]<<' '<<node->key[1]<<' '<<x_value<<' '<<y_value<<' '<<'\n';

}
