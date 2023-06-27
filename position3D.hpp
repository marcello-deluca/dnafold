#ifndef DEFINED_POSITION3D
#define DEFINED_POSITION3D
#include <vector>
template<typename T>
class position3D{
public:
  T x;
  T y;
  T z;

  std::vector<T> a;
  std::vector<T> b;
  std::vector<T> c;
  
  // Constructor: xyz.
  position3D(): x(0), y(0), z(0){
    for (size_t i = 0; i < 3; i++){
      a.push_back(0);
      b.push_back(0);
      c.push_back(0);
    }
  };
  position3D(T xIn, T yIn, T zIn): x(xIn), y(yIn), z(zIn){};

  // Set Coordinates.
  void set(T xIn, T yIn, T zIn){
    x = xIn;
    y = yIn;
    z = zIn;
  }

  void setA(T in1, T in2, T in3){
    a[0]=in1;
    a[1]=in2;
    a[2]=in3;
  }

  void setB(T in1, T in2, T in3){
    b[0]=in1;
    b[1]=in2;
    b[2]=in3;
  }

  void setC(T in1, T in2, T in3){
    c[0]=in1;
    c[1]=in2;
    c[2]=in3;
  }

  // Returns xyz vector of coords.
  std::vector<T> getCoords(){
    std::vector<T> answer;
    answer.push_back(x);
    answer.push_back(y);
    answer.push_back(z);
    return answer;
  }
};

#endif
