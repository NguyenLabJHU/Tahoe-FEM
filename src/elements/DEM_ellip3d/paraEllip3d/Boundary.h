#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "realtypes.h"
#include "Vec.h"
#include "Cylinder.h"
#include <map>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

namespace dem {
  
  class Particle; // forward declaration, only use pointer to class Particle
  
  /////////////////////////////////////
  class BoundaryTgt {  
  public:
    int  particleId;
    Vec  tgtForce;
    Vec  tgtDisp;
    bool tgtLoading;
    Vec  tgtDispStart;
    REAL tgtPeak;
    
  public:
  BoundaryTgt()
    :particleId(0), tgtForce(0), tgtDisp(0), tgtLoading(false), tgtDispStart(0), tgtPeak(0)
      {}
    
  BoundaryTgt(int _particleId, Vec _v1, Vec _v2, bool _b, Vec _v3, REAL _tp)
    :particleId(_particleId), tgtForce(_v1), tgtDisp(_v2), tgtLoading(_b), tgtDispStart(_v3), tgtPeak(_tp)
    {}
    
  };
  
  /////////////////////////////////////
  class BdryContact {
  public:
    Particle *ptcl;
    Vec point;
    Vec normal;
    Vec tangt;
    REAL penetr;

  public:
  BdryContact()
    :ptcl(NULL), point(0), normal(0), tangt(0), penetr(0) 
    {}
    
  BdryContact(Particle *p, Vec pt, Vec nm, Vec tg, REAL pntr)
    :ptcl(p), point(pt), normal(nm), tangt(tg), penetr(pntr) 
    {}
    
  void print(std::ostream &os) {
    os << std::setw(OWID) << point.getX()
       << std::setw(OWID) << point.getY()
       << std::setw(OWID) << point.getZ()
       << std::setw(OWID) << normal.getX()
       << std::setw(OWID) << normal.getY()
       << std::setw(OWID) << normal.getZ() 
       << std::setw(OWID) << tangt.getX()
       << std::setw(OWID) << tangt.getY()
       << std::setw(OWID) << tangt.getZ() 
       << std::setw(OWID) << penetr << std::endl;

  }
  private:
    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & ptcl;
      ar & point;
      ar & normal;
      ar & tangt;
      ar & penetr;
    }
  };

  ///////////////////////////////////////
  class Boundary { // abstract base class
  protected:
    int id;
    int type;

    std::vector<Particle *> possParticle;
    std::vector<BdryContact> contactInfo;
    int  contactNum;
    Vec  normal;
    Vec  tangt;
    REAL penetr;
    
  private:
    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & id;
      ar & type;
      ar & possParticle;
      ar & contactInfo;
      ar & contactNum;
      ar & normal;
      ar & tangt;
      ar & penetr;
    }
    
  public:
    Boundary(int i = 0, int tp = 0)
      :id(i), type(tp), contactNum(0), normal(0), tangt(0), penetr(0)  
      {}

    Boundary(int type, std::ifstream &ifs);
    virtual ~Boundary() {} // polymorphic base class requires a virtual destructor
    
    int getId() { return id; }
    int getType() { return type; }
    std::vector<Particle *> &getPossParticle () {return possParticle;}
    std::vector<BdryContact> &getContactInfo () {return contactInfo;}
    int getContactNum() const { return contactNum; }
    Vec getNormalForce() const { return normal; }
    Vec getTangtForce() const { return tangt; }
    REAL getAvgPenetr() const { return penetr; }

    virtual void print(std::ostream &os);
    virtual void printContactInfo(std::ostream &os);
    virtual void findBdryContact(std::vector<Particle *> &ptcls) = 0;
    virtual void boundaryForce(std::map<int,std::vector<BoundaryTgt> > &boundaryTgtMap) = 0;
    virtual void updateLocation() = 0;
    virtual void updateStatForce();
    void clearStatForce();
  };

  ///////////////////////////////////////
  class planeBoundary : public Boundary {
  private:
    Vec  direc;
    Vec  point;
    
  private:
    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & boost::serialization::base_object<Boundary >(*this);
      ar & direc;
      ar & point;
    }
    
  public:
    planeBoundary(int i = 0, int tp = 0)
      :Boundary(i, tp), direc(0), point(0) 
      {}

    planeBoundary(int type, std::ifstream &ifs);

    Vec getDirec() const { return direc; }
    Vec getPoint() const { return point; }

    REAL distanceToBdry(Vec pos) const { return (pos - point) % normalize(direc); }
    void print(std::ostream &os);
    void printContactInfo(std::ostream &os);

    void updateLocation() {}
    void findBdryContact(std::vector<Particle *> &ptcls);
    void boundaryForce(std::map<int,std::vector<BoundaryTgt> > &boundaryTgtMap);

  };
  
  ///////////////////////////////////////
  class cylinderBoundary : public Boundary {
  private:
    Vec  direc;
    Vec  point;
    REAL radius;
    
  private:
    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & boost::serialization::base_object<Boundary >(*this);
      ar & direc;
      ar & point;
      ar & radius;
    }
    
  public:
    cylinderBoundary()
      :Boundary(), direc(0), point(0), radius(0)
      {}

    cylinderBoundary(int type, std::ifstream &ifs);

    Vec getDirec() const { return direc; }
    Vec getPoint() const { return point; }
    REAL getRadius() const { return radius; }
 
    void print(std::ostream &os) {
      Boundary::print(os);
      os << std::setw(OWID) << direc.getX()
	 << std::setw(OWID) << direc.getY()
	 << std::setw(OWID) << direc.getZ()
	 << std::setw(OWID) << point.getX()
	 << std::setw(OWID) << point.getY()
	 << std::setw(OWID) << point.getZ()
	 << std::setw(OWID) << radius
	 << std::endl << std::endl;
    }

    void updateLocation() {}
    void findBdryContact(std::vector<Particle *> &ptcls);
    void boundaryForce(std::map<int,std::vector<BoundaryTgt> > &boundaryTgtMap);

  };

} // namespace dem ends

#endif
