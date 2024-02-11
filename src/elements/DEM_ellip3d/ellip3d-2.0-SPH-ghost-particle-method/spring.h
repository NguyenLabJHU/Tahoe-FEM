#ifndef SPRING_H
#define SPRING_H

#include "realtypes.h"
#include "particle.h"
#include "vec.h"
#include <cassert>

namespace dem{

class spring {

 public:
  spring(particle &p1, particle &p2, REAL young);
  spring(std::vector<particle*> &ParticleVec, int id1, int id2, REAL young);

  REAL getLength0() const {return length0;}
  REAL getLength() const {return vfabs( p2.getCurrPosition() - p1.getCurrPosition() );}
  vec  getDeformation(); 
  void applyForce();
  int getParticleId1() const {return p1.getID();}
  int getParticleId2() const {return p2.getID();}

 private:
  particle &p1;
  particle &p2;
  REAL Young;   // Young's modulus
  REAL ks;      // stiffness
  REAL length0; // equilibrium length

  void init(particle &p1, particle &p2) {
    length0 = vfabs( p2.getCurrPosition() - p1.getCurrPosition() );
    REAL radius = p1.getA();
    assert (radius == p2.getA() );
    ks = Young * 4 * radius * radius / length0;
  } 
};



}

#endif
