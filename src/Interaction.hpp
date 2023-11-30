#ifndef DEFINED_INTERACTION
#define DEFINED_INTERACTION

class Interaction{
public:
  int FirstParticleIndex;
  int SecondParticleIndex;

  //constructor:
  Interaction(): FirstParticleIndex(-1), SecondParticleIndex(-1){
  }

  Interaction(int First, int Second): FirstParticleIndex(First), SecondParticleIndex(Second){};
};

#endif
