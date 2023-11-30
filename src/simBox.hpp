#ifndef DEFINED_SIMBOX
#define DEFINED_SIMBOX
#include "dims3D.hpp"
template<typename T>
class simBox{
public:
  simBox(T xIn, T yIn, T zIn): dimensions(xIn, yIn, zIn){};
  dims3D<T> dimensions;

  void resize(T xin, T yin, T zin){
    dimensions.resize(xin, yin, zin);
  }
};

#endif
