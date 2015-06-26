#ifndef PARTICLE_H
#define PARTICLE_H

#include "realtypes.h"
#include "vec.h"

class particle{
  public:
    particle(int n, REAL* posi_x, REAL* posi_y, REAL* posi_z);
    
    int  getID() const {return ID;}

    REAL getVolume() const {return volume;}
    REAL getMass() const {return mass;}
    REAL getDensity() const {return density;}
    vec  getCurrPosition() const {return curr_position;}
    vec  getPrevPosition() const {return prev_position;}
    vec  getInitPosition() const {return init_position;}	// initial position for granular strain
    vec  getStartPosition() const {return start_position;}

    vec  getCurrVelocity() const {return curr_velocity;}
    vec  getPrevVelocity() const {return prev_velocity;}
    vec  getCurrOmga() const {return curr_omga;}
    vec  getPrevOmga() const {return prev_omga;}
    vec  getCurrAcceleration() const {return curr_acceleration;}
    vec  getPrevAcceleration() const {return prev_acceleration;}

    
    void setID(int n){ID=n;}

    void setCurrPosition(vec vv){curr_position=vv;}
    void setPrevPosition(vec vv){prev_position=vv;}
    void setInitPosition(vec vv){init_position=vv;}	// initial center for granular strain
    void setStartPosition(vec vv){start_position=vv;}	// position at starting step when tessellation is regenerated
    void setCurrVelocity(vec vv){curr_velocity=vv;}
    void setPrevVelocity(vec vv){prev_velocity=vv;}
    void setCurrOmga(vec vv){curr_omga=vv;}
    void setPrevOmga(vec vv){prev_omga=vv;}
    void setCurrAcceleration(vec vv){curr_acceleration=vv;}
    void setPrevAcceleration(vec vv){prev_acceleration=vv;}

    void update(int step, REAL time_interval);
    
    
  private:

    int ID;          
    REAL* position_x;
    REAL* position_y;
    REAL* position_z;

    vec curr_position;   // particle center
    vec prev_position;
    vec init_position;	// initial center
    vec start_position;	// position at starting step when tessellation is regenerated

    vec curr_velocity;   // the velocity of the mass center
    vec prev_velocity;
    vec curr_omga;       // angular velocity in global frame!
    vec prev_omga;
    vec curr_acceleration;
    vec prev_acceleration;

    REAL density; // specific gravity
    REAL mass;
    REAL volume;


};
  


#endif
