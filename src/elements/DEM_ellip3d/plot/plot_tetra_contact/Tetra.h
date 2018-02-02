#ifndef TETRA_H
#define TETRA_H

#include "realtypes.h"
#include "Vec.h"
#include <map>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

class Tetra {

 private:
  int i, j, k, l;    // node number of a cell (tetrahedron)
  std::map<int, Vec> ptclMap; // map between node number and particle coordinates

 public:
 Tetra()
   : i(0), j(0), k(0), l(0) {
  }
    
 Tetra(int d1, int d2, int d3, int d4, Vec v1, Vec v2, Vec v3, Vec v4)
   : i(d1), j(d2), k(d3), l(d4) {
    ptclMap[i] = v1;
    ptclMap[j] = v2;
    ptclMap[k] = v3;
    ptclMap[l] = v4;
  }
    
  int getI() const {return i;}
  int getJ() const {return j;}
  int getK() const {return k;}
  int getL() const {return l;}

  Vec getParticleI() {return ptclMap[i];}
  Vec getParticleJ() {return ptclMap[j];}
  Vec getParticleK() {return ptclMap[k];}
  Vec getParticleL() {return ptclMap[l];}

  REAL getVolume();
  void setNodeOrder();
  std::vector<REAL> getAngles();

};

#endif
