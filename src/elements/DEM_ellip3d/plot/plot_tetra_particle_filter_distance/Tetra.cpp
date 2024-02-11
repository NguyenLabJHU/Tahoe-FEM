#include "Tetra.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <algorithm>

const double PI = 3.1415927;

REAL Tetra::getVolume() {
  REAL x1, y1, z1;
  REAL x2, y2, z2;
  REAL x3, y3, z3;
  REAL x4, y4, z4;
  x1 = ptclMap[i].coord.getX();
  y1 = ptclMap[i].coord.getY();
  z1 = ptclMap[i].coord.getZ();
  x2 = ptclMap[j].coord.getX();
  y2 = ptclMap[j].coord.getY();
  z2 = ptclMap[j].coord.getZ();
  x3 = ptclMap[k].coord.getX();
  y3 = ptclMap[k].coord.getY();
  z3 = ptclMap[k].coord.getZ();
  x4 = ptclMap[l].coord.getX();
  y4 = ptclMap[l].coord.getY();
  z4 = ptclMap[l].coord.getZ();

  return (x2*y3*z4-x2*y4*z3-x3*y2*z4+x3*y4*z2+x4*y2*z3-x4*y3*z2
	  -x1*y3*z4+x1*y4*z3+x3*y1*z4-x3*y4*z1-x4*y1*z3+x4*y3*z1
	  +x1*y2*z4-x1*y4*z2-x2*y1*z4+x2*y4*z1+x4*y1*z2-x4*y2*z1
	  -x1*y2*z3+x1*y3*z2+x2*y1*z3-x2*y3*z1-x3*y1*z2+x3*y2*z1)/6.0;
}

void Tetra::setNodeOrder() {
  if(getVolume() < 0) { // swap node 2 and node 3
    /*
      std::cout << "tetra nodes are not numbered counter-clockwise, volume=" << getVolume() << " iteration=" << iteration << std::endl
      << std::setw(OWID) << i << std::setw(OWID) << j << std::setw(OWID) << k << std::setw(OWID) << l << std::endl
      << std::setw(OWID) << ptclMap[i].coord.getX() << std::setw(OWID) << ptclMap[i].coord.getY() << std::setw(OWID) << ptclMap[i].coord.getZ() << std::endl
      << std::setw(OWID) << ptclMap[j].coord.getX() << std::setw(OWID) << ptclMap[j].coord.getY() << std::setw(OWID) << ptclMap[j].coord.getZ() << std::endl
      << std::setw(OWID) << ptclMap[k].coord.getX() << std::setw(OWID) << ptclMap[k].coord.getY() << std::setw(OWID) << ptclMap[k].coord.getZ() << std::endl
      << std::setw(OWID) << ptclMap[l].coord.getX() << std::setw(OWID) << ptclMap[l].coord.getY() << std::setw(OWID) << ptclMap[l].coord.getZ() << std::endl;
    */
    int jTmp = j;
    j = k;
    k = jTmp;
  }
}

typedef std::map<int, Particle> M;

bool value_comparer(M::value_type &i1, M::value_type &i2) {
  return i1.second.coord.getZ() < i2.second.coord.getZ();
}

std::vector<REAL> Tetra::getInfo() {
  M::iterator iter = std::max_element(ptclMap.begin(), ptclMap.end(), value_comparer);
  //std::cout << i << " " << j << " " << k << " " << l << ", " << iter->first << ":" << iter->second.getZ() << std::endl << std::endl;
  // now iter->first is the top vertex, then ensure i is the top vertex to use the equation of top solid angle
  if (iter->first == i) {
    ;
  } else if (iter->first == j) {
    int tmp = i;
    i = j;
    j = tmp;
  } else if (iter->first == k) {
    int tmp = i;
    i = k;
    k = tmp;
  } else if (iter->first == l) {
    int tmp = i;
    i = l;
    l = tmp;
  }
  setNodeOrder(); // this does not affect node i, which is always the top vertex.

  REAL D = 288.0 * pow(getVolume(), 2);
  REAL eij = vfabs(ptclMap[i].coord - ptclMap[j].coord);
  REAL eik = vfabs(ptclMap[i].coord - ptclMap[k].coord);
  REAL eil = vfabs(ptclMap[i].coord - ptclMap[l].coord);
  REAL ejk = vfabs(ptclMap[j].coord - ptclMap[k].coord);
  REAL ekl = vfabs(ptclMap[k].coord - ptclMap[l].coord);
  REAL ejl = vfabs(ptclMap[j].coord - ptclMap[l].coord);

  REAL Ni = (eij + eik) * (eik + eil) * (eil + eij) - (eij*ekl*ekl + eik*ejl*ejl + eil*ejk*ejk);
  //if (Ni <= 0 ) std::cout << eij << " " << eik << " " << eil << " " << ejk << " " << ejl << " " << ekl << " Ni= " << Ni << " D=" << D << " V=" << getVolume() << std::endl;
  REAL phi, theta;
  if (Ni <= 0) {
    phi = 2 * PI;
    theta = 180;
  } else {
    phi = 2 * atan(sqrt(D/2) / Ni); // solid angle in sr
    theta = acos(1 - phi/(2*PI)) *2 *180/PI; // theta x2 (not theta) or total angle of a cone in degree
  }
  std::vector<REAL> info;
  info.push_back((eij + eik + eil + ejk + ekl + ejl)/6);
  info.push_back((eij + eik + eil)/3);
  info.push_back((ejk + ekl + ejl)/3);
  info.push_back((eij + eik + eil)/(ejk + ekl + ejl));
  info.push_back(phi);
  info.push_back(theta);
  info.push_back(getVolume());
 
  return info;
}

bool Tetra::isValid(double gapRatio) {
  REAL eij = vfabs(ptclMap[i].coord - ptclMap[j].coord);
  REAL eik = vfabs(ptclMap[i].coord - ptclMap[k].coord);
  REAL eil = vfabs(ptclMap[i].coord - ptclMap[l].coord);
  REAL ejk = vfabs(ptclMap[j].coord - ptclMap[k].coord);
  REAL ekl = vfabs(ptclMap[k].coord - ptclMap[l].coord);
  REAL ejl = vfabs(ptclMap[j].coord - ptclMap[l].coord);

  REAL dij = (ptclMap[i].radius + ptclMap[j].radius);
  REAL dik = (ptclMap[i].radius + ptclMap[k].radius);
  REAL dil = (ptclMap[i].radius + ptclMap[l].radius);
  REAL djk = (ptclMap[j].radius + ptclMap[k].radius);
  REAL dkl = (ptclMap[k].radius + ptclMap[l].radius);
  REAL djl = (ptclMap[j].radius + ptclMap[l].radius);
 
  //gapRatio is defined as: d - (a1+a2) < gapRatio * d
  if (eij/dij > 1/(1-gapRatio)) return false;
  if (eik/dik > 1/(1-gapRatio)) return false;
  if (eil/dil > 1/(1-gapRatio)) return false;
  if (ejk/djk > 1/(1-gapRatio)) return false;
  if (ekl/dkl > 1/(1-gapRatio)) return false;
  if (ejl/djl > 1/(1-gapRatio)) return false;

  return true;
}
