#ifndef VERLET_LIST_H
#define VERLET_LIST_H

class VerletList{

public:
  //constructor:
  std::vector<size_t> NeighborList;
  size_t ParticleID;
  VerletList(size_t ID): NeighborList(), ParticleID(ID){};
  void AddToList(size_t ParticleNum){
    NeighborList.push_back(ParticleNum);
  }
  std::vector<size_t> getNeighborList(){
    return NeighborList;
  }
private:


};


#endif
