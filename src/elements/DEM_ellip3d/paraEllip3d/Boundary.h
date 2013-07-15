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
  
  
  typedef struct bdryfunc {
    int  order; // 1-linear; 2-quadratic
    Vec  dirc;  // normal vector if plane, mother line vector if cylinder,it points out of the particles		
    Vec  apt;   // a point on the plane or a point on the axis of the cylinder
    REAL rad;   // zero if plane
    
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
      ar & order;
      ar & dirc;
      ar & apt;
      ar & rad;
    }
    
  public:
    void display() const{
      std::cout << "order: " << order << std::endl;
      std::cout << "dirc: " << dirc.getX() << " " << dirc.getY() << " " << dirc.getZ() << std::endl;
      std::cout << "apt: " << apt.getX() << " " << apt.getY() << " " << apt.getZ() << std::endl;
      std::cout << "radius: " << rad << std::endl;
    }
    
    void display(std::ofstream &ofs) const{
      ofs << std::setw(OWID) << order
	  << std::setw(OWID) << dirc.getX()
	  << std::setw(OWID) << dirc.getY()
	  << std::setw(OWID) << dirc.getZ()
	  << std::setw(OWID) << apt.getX()
	  << std::setw(OWID) << apt.getY()
	  << std::setw(OWID) << apt.getZ()
	  << std::setw(OWID) << rad << std::endl;
    }
  } BdryCoef;
  
  
  typedef struct updatectl{
    Vec  tran;  // tranlation second
    Vec  rote;  // rotate first
    Vec  fixpt; // before any update is made
    REAL expnd; // expand last
    
    updatectl() {tran=0;rote=0;fixpt=0;expnd=1;}
    
    void display() const {
      std::cout << "tran: " << tran.getX() << " " << tran.getY() << " " << tran.getZ() << std::endl;
      std::cout << "rote: " << rote.getX() << " " << rote.getY() << " " << rote.getZ() << std::endl;
      std::cout << "fixpt: " << fixpt.getX() << " " << fixpt.getY() << " " << fixpt.getZ() << std::endl;
      std::cout << "expand:" << expnd << std::endl;
    }
    
  } UPDATECTL;
  
  
  class Boundary {
  public:
    int  boundaryId; // the first record defines the bdry itself, the other 
    std::vector<BdryCoef> coefOfLimits; // limitnum records define the other lines on the bdry 
    REAL avgNormal;  // that give limits to the first boundary.
    REAL avgPenetr;  // average penetration by particles onto this boundary
    int  contactNum; // contact numbers by particles onto this boundary
    REAL area;       // bounary's area
    int  limitNum;   // how many lines the boundary has
    
  private:
    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & boundaryId;
      ar & coefOfLimits;
      ar & avgNormal;
      ar & avgPenetr;
      ar & contactNum;
      ar & area;
      ar & limitNum;
    }
    
  public:
    Boundary() {}
    Boundary(std::ifstream &ifs);
    int getBdryID() { return boundaryId; }
    virtual ~Boundary() {} // base class needs a virtual destructor.
    virtual void display() const {
      std::cout << "area: " << area << " limitNum: " << limitNum << std::endl;
      std::vector<BdryCoef>::const_iterator it;
      for(it = coefOfLimits.begin(); it != coefOfLimits.end(); ++it)
	(*it).display();
    }
    virtual void display(std::ofstream &ofs) const {
      std::vector<BdryCoef>::const_iterator it;
      ofs << std::endl
	  << std::setw(OWID) << (*coefOfLimits.begin()).order << std::endl;
      ofs << std::setw(OWID) << boundaryId
	  << std::setw(OWID) << limitNum << std::endl;
      for(it = coefOfLimits.begin(); it != coefOfLimits.end(); ++it)
	(*it).display(ofs);
    }
    virtual void findBdryContact(std::vector<Particle *> &ptcls) {}
    
    // calculate for each boundary particles the rigid boundary force
    virtual void boundaryForce(std::map<int,std::vector<BoundaryTgt> > &boundaryTgtMap) {}
    
    virtual Vec getNormalForce() const { return 0; }
    virtual REAL getAvgNormal() const { return 0; }
    virtual REAL getAvgPenetr() const { return 0; }
    virtual int getCntnum() const { return 0; }
    virtual Vec getShearForce() const { return 0; }
    virtual Vec getApt() const { return 0; }
    virtual Vec getDirc() const { return 0; }
    virtual void setArea(REAL a) { area=a; }
    virtual REAL getArea() { return area; }
    virtual void update(UPDATECTL &ctl); //the boundary is translating with tran and rotating with rote around axis
  };
  
  class plnBoundary : public Boundary {
  public:
    Vec normal;  // normal force acting on the boundary by all contacting particles 
    Vec tangt;   // tangential force acting on the boundary
    Vec moment;  // moment on the boundary
    std::vector<Particle *> possBdryParticle; // possible boundary particles of this specific boundary
    
  private:
    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & boost::serialization::base_object<Boundary >(*this);
      ar & normal;
      ar & tangt;
      ar & moment;
      ar & possBdryParticle;
    }
    
  public:
    plnBoundary() {}
    plnBoundary(std::ifstream &ifs):Boundary(ifs) {
      normal = 0;
      tangt = 0;
      moment = 0;
      this->avgNormal = 0;
      this->avgPenetr = 0;
      this->contactNum = 0;
    }
    
    int getBdryID() {return this->boundaryId;}
    void display() const;
    REAL distToBdry(Vec posi) const;
    void findBdryContact(std::vector<Particle *> &ptcls);
    Vec getApt() const;
    Vec getDirc() const;
    plnBoundary* getBdry(int bdryid) { return this; }
    const plnBoundary* getBdry(int bdryid) const { return this; }
    void boundaryForce(std::map<int,std::vector<BoundaryTgt> > &boundaryTgtMap);
    Vec getNormalForce() const { return normal; }
    REAL getAvgNormal() const { return this->avgNormal; }
    REAL getAvgPenetr() const { return this->avgPenetr; }
    int getCntnum() const { return this->contactNum; }
    Vec getShearForce() const { return tangt; }
  };
  
  class cylBoundary : public Boundary {
  public:
    Vec normal; 
    std::vector<Particle *> possBdryParticle;

  public:
  cylBoundary(std::ifstream &ifs):Boundary(ifs){normal = 0;}
    void display() const;
    REAL distToBdry(Vec posi) const;
    void findBdryContact(std::vector<Particle *> &ptcls);
    void boundaryForce();
    Vec getNormalForce() const{ return normal; }
  };
  
} // namespace dem ends

#endif
