#ifndef RECTANGLE_H
#define RECTANGLE_H

#include "realtypes.h"
#include "Vec.h"

namespace dem {
  
  class Rectangle {
    
  public:
    Rectangle()
      :dimx(0), dimy(0), dimz(0), center(0), v1(0), v2(0)
      {}
    
    Rectangle(REAL dx, REAL dy, REAL dz, Vec c)
      :dimx(dx), dimy(dy), dimz(dz), center(c) {
      v1.setx( c.getx() - dx/2.0);
      v1.sety( c.gety() - dy/2.0);
      v1.setz( c.getz() - dz/2.0);
      v2.setx( c.getx() + dx/2.0);
      v2.sety( c.gety() + dy/2.0);
      v2.setz( c.getz() + dz/2.0);
    }
		     
    Rectangle(REAL dx, REAL dy, REAL dz, Vec ref, int i)
      :dimx(dx), dimy(dy), dimz(dz), v1(ref){
      center.setx(v1.getx() + dx/2.0);
      center.sety(v1.gety() + dy/2.0);
      center.setz(v1.getz() + dz/2.0);
      v2.setx(v1.getx() + dx);
      v2.sety(v1.gety() + dy);
      v2.setz(v1.getz() + dz);
    }
  
    Rectangle(REAL x1, REAL y1, REAL z1, REAL x2, REAL y2, REAL z2)
      :dimx(x2-x1), dimy(y2-y1), dimz(z2-z1), v1(x1, y1, z1), v2(x2, y2, z2) {
      center = (v1 + v2) / 2;
    }

    Rectangle(Vec _v1, Vec _v2)
      :v1(_v1), v2(_v2) {
      center = (v1 + v2) / 2;
      Vec vt = v2 - v1;
      dimx = vt.getx();
      dimy = vt.gety();
      dimz = vt.getz();
    }

    REAL getDimx() const {return dimx;}
    REAL getDimy() const {return dimy;}
    REAL getDimz() const {return dimz;}
    Vec  getCenter() const {return center;}
    Vec  getMinCorner() const {return v1;}
    Vec  getMaxCorner() const {return v2;}
    REAL getVolume() const {return dimx*dimy*dimz;}
    
    void setDimx(REAL dx) {dimx=dx;}
    void setDimy(REAL dy) {dimy=dy;}
    void setDimz(REAL dz) {dimz=dz;}
    void setCenter(Vec v) {center=v;}
    void setV1(Vec v) {v1=v;}
    void setV2(Vec v) {v2=v;}
    Vec  randomPoint() const;
    void print() const;
    void set(REAL dx, REAL dy, REAL dz, Vec c) {
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
    Vec  center;
    Vec  v1; // lower corner of the rectangle, minimum x,y,z value
    Vec  v2; // lower corner of the rectangle, maximum x,y,z value  

    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & dimx;
      ar & dimy;
      ar & dimz;
      ar & center;
      ar & v1;
      ar & v2;
    }

  };
  
} // namespace dem

#endif
