#ifndef RECTANGLE_H
#define RECTANGLE_H

#include "realtypes.h"
#include "vec.h"
#include "shape.h"

namespace dem {
  
  class rectangle: public shape {
    
  public:
    rectangle()
      :dimx(0),dimy(0),dimz(0),center(0),v1(0),v2(0)
      {}
    
    rectangle(REAL dx, REAL dy, REAL dz, vec c)
      :dimx(dx),dimy(dy),dimz(dz),center(c) {
      v1.setx( c.getx() - dx/2.0);
      v1.sety( c.gety() - dy/2.0);
      v1.setz( c.getz() - dz/2.0);
      v2.setx( c.getx() + dx/2.0);
      v2.sety( c.gety() + dy/2.0);
      v2.setz( c.getz() + dz/2.0);
    }
		     
    rectangle(REAL dx, REAL dy, REAL dz, vec ref, int i)
      :dimx(dx),dimy(dy),dimz(dz),v1(ref){
      center.setx(v1.getx() + dx/2.0);
      center.sety(v1.gety() + dy/2.0);
      center.setz(v1.getz() + dz/2.0);
      v2.setx(v1.getx() + dx);
      v2.sety(v1.gety() + dy);
      v2.setz(v1.getz() + dz);
    }
  
    rectangle(rectangle &rec) {
      dimx = rec.getDimx();
      dimy = rec.getDimy();
      dimz = rec.getDimz();
      center.set(rec.getCenter());
      v1.set(rec.getMinCorner());
      v2.set(rec.getMaxCorner());
    }
    
    REAL getDimx() const {return dimx;}
    REAL getDimy() const {return dimy;}
    REAL getDimz() const {return dimz;}
    vec  getCenter() const {return center;}
    vec  getMinCorner() const {return v1;}
    vec  getMaxCorner() const {return v2;}
    REAL getVolume() const {return dimx*dimy*dimz;}
    
    void setDimx(REAL dx) {dimx=dx;}
    void setDimy(REAL dy) {dimy=dy;}
    void setDimz(REAL dz) {dimz=dz;}
    void setCenter(vec v) {center=v;}
    void setV1(vec v) {v1=v;}
    void setV2(vec v) {v2=v;}
    vec  randomPoint() const;
    void print() const;
    void set(REAL dx, REAL dy, REAL dz, vec c) {
      dimx = dx;
      dimy = dy;
      dimz = dz;
      center = c;
      v1.setx( c.getx() - dx/2.0);
      v1.sety( c.gety() - dy/2.0);
      v1.setz( c.getz() - dz/2.0);
      v2.setx( c.getx() + dx/2.0);
      v2.sety( c.gety() + dy/2.0);
      v2.setz( c.getz() + dz/2.0);
    }
    
  private:
    REAL dimx;
    REAL dimy;
    REAL dimz;
    vec  center;
    vec  v1; // lower corner of the rectangle, minimum x,y,z value
    vec  v2; // lower corner of the rectangle, maximum x,y,z value

  };
  
} // namespace dem

#endif
