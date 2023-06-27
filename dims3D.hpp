#ifndef DEFINED_DIMS3D
#define DEFINED_DIMS3D
template<typename T>
class dims3D{
public:
  dims3D(T xIn, T yIn, T zIn): x(xIn), y(yIn), z(zIn){};
  T x;
  T y;
  T z;
  void resize(T x_new, T y_new, T z_new){
    x=x_new;
    y=y_new;
    z=z_new;
  }
};
#endif
