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
#include <cstddef>
#include <cstdlib>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

namespace dem {
  
  class Particle; // forward declaration, only use pointer to class Particle
  
  /////////////////////////////////////
  class BoundaryTgt {  
  public:
    std::size_t  particleId;
    Vec  tgtForce;
    Vec  tgtDisp;
    bool tgtLoading;
    Vec  tgtDispStart;
    REAL tgtPeak;
    bool tgtSlide;
    bool tgtRoll;
    
  public:
    BoundaryTgt();
    BoundaryTgt(std::size_t _particleId, Vec _v1, Vec _v2, bool _b, Vec _v3, REAL _tp, bool _s, bool _r);
    
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
    BdryContact();
    BdryContact(Particle *p, Vec pt, Vec nm, Vec tg, REAL pntr);
    ~BdryContact();
    void print(std::ostream &os);

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

  /////////////////////////////////////
  class Plane {
  public:
    Vec direc;
    Vec point;

  public:
    Plane();
    Plane(Vec dir, Vec pt);
    void print(std::ostream &os);
    Vec getDirec() const;
    Vec getPoint() const;

  private:
    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & direc;
      ar & point;
    }
  };

  ///////////////////////////////////////
  class Boundary { // abstract base class
  protected:
    std::size_t id;
    std::size_t type;

    // extra edges that are necessary to define a finite plane
    // e.g., side wall of a top-open container
    std::size_t extraNum;
    std::vector<Plane> extraEdge;

    std::vector<Particle *> possParticle;
    std::vector<BdryContact> contactInfo;
    std::size_t  contactNum;
    Vec  normal;
    Vec  tangt;
    REAL penetr;
    
  private:
    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & id;
      ar & type;
      ar & extraNum;
      ar & extraEdge;
      ar & possParticle;
      ar & contactInfo;
      ar & contactNum;
      ar & normal;
      ar & tangt;
      ar & penetr;
    }
    
  public:
    Boundary(std::size_t i = 0, std::size_t tp = 0, std::size_t en = 0);
    Boundary(std::size_t type, std::ifstream &ifs);
    virtual ~Boundary() {} // polymorphic base class requires a virtual destructor
    
    std::size_t getId();
    std::size_t getType();
    std::vector<Particle *> &getPossParticle();
    std::vector<BdryContact> &getContactInfo();
    std::size_t getContactNum() const;
    Vec  getNormalForce() const;
    Vec  getTangtForce() const;
    REAL getAvgPenetr() const;

    virtual void print(std::ostream &os);
    virtual void printContactInfo(std::ostream &os);
    virtual void findBdryContact(std::vector<Particle *> &ptcls) = 0;
    virtual void boundaryForce(std::map<std::size_t,std::vector<BoundaryTgt> > &boundaryTgtMap) = 0;
    virtual void updateStatForce();
    void clearStatForce();
    void clearContactInfo();

    virtual void updateIsotropic(REAL simga, REAL areaX, REAL areaY, REAL areaZ) {}
    virtual void updateOedometer(REAL simga, REAL areaX, REAL areaY, REAL areaZ) {}
    virtual void updateTriaxial(REAL simga, REAL areaX, REAL areaY, REAL areaZ) {}
    virtual void updatePlaneStrain(REAL simga, REAL areaX, REAL areaY, REAL areaZ) {}
    virtual void updateTrueTriaxial(REAL simga, REAL areaX, REAL areaY, REAL areaZ, REAL sigmaX, REAL sigmaY) {}

    virtual Vec  getPoint() const = 0;
    virtual Vec  getVeloc() const = 0;
    virtual Vec  getPrevPoint() const = 0;
    virtual Vec  getPrevVeloc() const = 0;
    virtual void setPoint(Vec pnt) = 0;
    virtual void setVeloc(Vec vel) = 0;
  };

  ///////////////////////////////////////
  class planeBoundary : public Boundary {
  private:
    Vec  direc;
    Vec  point;
    Vec  prevPoint;
    Vec  veloc;
    Vec  prevVeloc;
    
  private:
    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & boost::serialization::base_object<Boundary >(*this);
      ar & direc;
      ar & point;
      ar & prevPoint;
      ar & veloc;
      ar & prevVeloc;
    }
    
  public:
    planeBoundary(std::size_t i = 0, std::size_t tp = 0, std::size_t en = 0);
    planeBoundary(std::size_t type, std::ifstream &ifs);

    Vec getDirec() const;
    Vec getPoint() const;
    Vec getVeloc() const;
    Vec getPrevPoint() const;
    Vec getPrevVeloc() const;

    void setDirec(Vec dir);
    void setPoint(Vec pnt);
    void setVeloc(Vec vel);

    REAL distanceToBdry(Vec pos) const;
    REAL distanceToBdry(Vec pos, Plane pn) const;

    void print(std::ostream &os);
    void printContactInfo(std::ostream &os);

    void updateIsotropic(REAL sigma, REAL areaX, REAL areaY, REAL areaZ);
    void updateOedometer(REAL simga, REAL areaX, REAL areaY, REAL areaZ);
    void updateTriaxial(REAL simga, REAL areaX, REAL areaY, REAL areaZ);
    void updatePlaneStrain(REAL simga, REAL areaX, REAL areaY, REAL areaZ);
    void updateTrueTriaxial(REAL simga, REAL areaX, REAL areaY, REAL areaZ, REAL sigmaX, REAL sigmaY);
    void findBdryContact(std::vector<Particle *> &ptcls);
    void boundaryForce(std::map<std::size_t,std::vector<BoundaryTgt> > &boundaryTgtMap);

  };
  
  ///////////////////////////////////////
  class cylinderBoundary : public Boundary {
  private:
    Vec  direc;
    Vec  point;
    Vec  prevPoint;
    Vec  veloc;
    Vec  prevVeloc;    
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
    cylinderBoundary();
    cylinderBoundary(std::size_t type, std::ifstream &ifs);

    Vec getDirec() const;
    Vec getPoint() const;
    Vec getVeloc() const;
    Vec getPrevPoint() const;
    Vec getPrevVeloc() const;
    REAL getRadius() const;

    void setDirec(Vec dir);
    void setPoint(Vec pnt);
    void setVeloc(Vec vel);
 
    void print(std::ostream &os);

    void findBdryContact(std::vector<Particle *> &ptcls);
    void boundaryForce(std::map<std::size_t,std::vector<BoundaryTgt> > &boundaryTgtMap);

  };

} // namespace dem ends

#endif
