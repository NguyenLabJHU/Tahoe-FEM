#ifndef SPRING_H
#define SPRING_H

#include "realtypes.h"
#include "Particle.h"
#include "Vec.h"
#include <cassert>

namespace dem{
  
  class Spring {
    
  public:
    Spring(Particle &p1, Particle &p2, REAL young);
    Spring(std::vector<Particle*> &ParticleVec, int id1, int id2, REAL young);
    
    REAL getLength0() const {return length0;}
    REAL getLength() const {return vfabs( p2.getCurrPos() - p1.getCurrPos() );}
    Vec  getDeformation(); 
    void applyForce();
    int getParticleId1() const {return p1.getId();}
    int getParticleId2() const {return p2.getId();}
    
  private:
    Particle &p1;
    Particle &p2;
    REAL Young;   // Young's modulus
    REAL ks;      // stiffness
    REAL length0; // equilibrium length
    
    void init(Particle &p1, Particle &p2) {
      length0 = vfabs( p2.getCurrPos() - p1.getCurrPos() );
      REAL radius = p1.getA();
      assert (radius == p2.getA() );
      ks = Young * 4 * radius * radius / length0;
    } 

  };
  
}

#endif
