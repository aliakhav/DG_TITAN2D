#ifndef SCALE_H
#define SCALE_H
#include <iostream>
#include <fstream>
using namespace std;
struct ScaleValues {
  double lengthscale;
  double heightscale;
  double gravityscale;

  ScaleValues() {
    std::ifstream inDatafile("scale.data", ios::in);
    if(inDatafile.fail())
      {
	// assume no scaling then...
	lengthscale = 1;
	heightscale = 1;
	gravityscale = 1;
      }
    inDatafile>>lengthscale;
    if(lengthscale < GEOFLOW_TINY)
      lengthscale = 1;
    inDatafile>>heightscale;
    if(heightscale < GEOFLOW_TINY)
      heightscale = 1;
    inDatafile>>heightscale;
    if(gravityscale < GEOFLOW_TINY)
      gravityscale = 1;
    
    inDatafile.close();
  }
};

#endif
