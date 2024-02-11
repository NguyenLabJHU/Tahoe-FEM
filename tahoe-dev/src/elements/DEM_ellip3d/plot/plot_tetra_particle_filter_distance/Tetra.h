#ifndef TETRA_H
#define TETRA_H

#include "realtypes.h"
#include "Vec.h"
#include <utility>
#include <map>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

class Particle {
 public:
  Vec coord;
  double radius; // maximum radius
};

class Tetra {

 private:
  int i, j, k, l;    // node number of a cell (tetrahedron)
  std::map<int, Particle> ptclMap; // map between node number and particle coordinates and radius.

 public:
 Tetra()
   : i(0), j(0), k(0), l(0) {
  }
    
 Tetra(int d1, int d2, int d3, int d4, Vec v1, Vec v2, Vec v3, Vec v4, double a1, double a2, double a3, double a4)
   : i(d1), j(d2), k(d3), l(d4) {
    ptclMap[i].coord = v1;
    ptclMap[j].coord = v2;
    ptclMap[k].coord = v3;
    ptclMap[l].coord = v4;
    ptclMap[i].radius = a1;
    ptclMap[j].radius = a2;
    ptclMap[k].radius = a3;
    ptclMap[l].radius = a4;
  }
    
  int getI() const {return i;}
  int getJ() const {return j;}
  int getK() const {return k;}
  int getL() const {return l;}

  Particle getParticleI() {return ptclMap[i];}
  Particle getParticleJ() {return ptclMap[j];}
  Particle getParticleK() {return ptclMap[k];}
  Particle getParticleL() {return ptclMap[l];}

  REAL getVolume();
  void setNodeOrder();
  bool isValid(double gapRatio);
  std::vector<REAL> getInfo();

};

#endif
