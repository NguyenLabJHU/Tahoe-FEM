// 1. This is a template class, for which we have to include the implementation in the header file.
//    As we cannot put using statement in a header file, we have to use std::xxx wherever we need
//     to reference anything from standard namespace.
//
// 2. A base class needs a virtual destructor, otherwise it may cause undetermined errors.
//
// 3. When inheritating a template class, it is important to know the idea "Name lookup, templates,
//    and accessing members of base classes". 
//    Reference: http://gcc.gnu.org/onlinedocs/gcc-4.0.2/gcc/Name-lookup.html#Name-lookup

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

namespace dem {
  
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

  
  typedef struct bdryfunc{
    int  order; // 1-linear; 2-quadratic
    Vec  dirc;  // normal vector if plane, mother line vector if cylinder,it points out of the particles		
    Vec  apt;   // a point on the plane or a point on the axis of the cylinder
    REAL rad;   // zero if plane
    int  side;  // zero if plane; side=1, particles are outside the cylinder; =-1, inside the cylinder
    
    void display() const{
      std::cout << "order: " << order << std::endl;
      std::cout << "dirc: " << dirc.getX() << " " << dirc.getY() << " " << dirc.getZ() << std::endl;
      std::cout << "apt: " << apt.getX() << " " << apt.getY() << " " << apt.getZ() << std::endl;
      std::cout << "radius: " << rad << " side: " << side << std::endl;
    }
    
    void display(std::ofstream &ofs) const{
      ofs << std::setw(OWID) << order
	  << std::setw(OWID) << dirc.getX()
	  << std::setw(OWID) << dirc.getY()
	  << std::setw(OWID) << dirc.getZ()
	  << std::setw(OWID) << apt.getX()
	  << std::setw(OWID) << apt.getY()
	  << std::setw(OWID) << apt.getZ()
	  << std::setw(OWID) << rad
	  << std::setw(OWID) << side << std::endl;
    }
  } BdryCoef;
  
  
  typedef struct updatectl{
    Vec  tran;  // tranlation second
    Vec  rote;  // rotate first
    Vec  fixpt; // before any update is made
    REAL expnd; // expand last
    
    updatectl() {tran=0;rote=0;fixpt=0;expnd=1;}
    
    void display() const{
      std::cout << "tran: " << tran.getX() << " " << tran.getY() << " " << tran.getZ() << std::endl;
      std::cout << "rote: " << rote.getX() << " " << rote.getY() << " " << rote.getZ() << std::endl;
      std::cout << "fixpt: " << fixpt.getX() << " " << fixpt.getY() << " " << fixpt.getZ() << std::endl;
      std::cout << "expand:" << expnd << std::endl;
    };
    
  } UPDATECTL;
  
  
  template<class T> class Boundary {
    
 public:
  int  boundaryId; // the first record defines the bdry itself, the other 
  std::vector<BdryCoef> CoefOfLimits; // limitnum records define the other lines on the bdry 
  REAL avgNormal;  // that give limits to the first boundary.
  REAL avgPenetr;  // average penetration by particles onto this boundary
  int  contactNum; // contact numbers by particles onto this boundary
  REAL area;       // bounary's area
  int  limitNum;   // how many lines the boundary have

 public:
  Boundary(std::ifstream &ifs);
  int getBdryID() {return boundaryId;}
  virtual ~Boundary() {} // base class needs a virtual destructor.
  virtual void display() const{
    std::cout << "area: " << area << " limitNum: " << limitNum << std::endl;
    std::vector<BdryCoef>::const_iterator it;
    for(it=CoefOfLimits.begin();it!=CoefOfLimits.end();++it)
      (*it).display();
  }
  virtual void display(std::ofstream &ofs) const{
    std::vector<BdryCoef>::const_iterator it;
    ofs << std::endl
	<< std::setw(OWID) << (*CoefOfLimits.begin()).order << std::endl;
    ofs << std::setw(OWID) << boundaryId
	<< std::setw(OWID) << limitNum << std::endl;
    for(it=CoefOfLimits.begin();it!=CoefOfLimits.end();++it)
      (*it).display(ofs);
  }
  virtual void findParticleOnBoundary(std::vector<T*>& ptcls){};

  // calculate for each boundary particles the rigid boundary force
  virtual void boundaryForce(std::map<int,std::vector<BoundaryTgt> >& boundaryTgtMap) {}

  virtual Vec getNormalForce() const{return 0;}
  virtual REAL getAvgNormal() const{return 0;}
  virtual REAL getAvgPenetr() const{return 0;}
  virtual int getCntnum() const{return 0;}
  virtual Vec getShearForce() const{return 0;}
  virtual Vec getApt() const{return 0;}
  virtual Vec getDirc() const{return 0;}
  virtual void setArea(REAL a){area=a;}
  virtual REAL getArea(){return area;}
  virtual void update(UPDATECTL& ctl); //the boundary is translating with tran and rotating with rote around axis
};
 
template <class T>
Boundary<T>::Boundary(std::ifstream &ifs){
  BdryCoef tmp;
  REAL x,y,z;
  CoefOfLimits.clear();
  ifs >> boundaryId >> limitNum;
  for (int k = 0; k < limitNum; k++){
    ifs >> tmp.order >> x >> y >> z;
    tmp.dirc=Vec(x,y,z);
    ifs >> x >> y >> z;
    tmp.apt=Vec(x,y,z);
    ifs >> tmp.rad >> tmp.side;
    CoefOfLimits.push_back(tmp);
  }
};

template<class T>
void Boundary<T>::update(UPDATECTL& ctl){
  BdryCoef tmp;
  std::vector<BdryCoef>::iterator it;
  Vec nv, napt;
  for (it=CoefOfLimits.begin();it!=CoefOfLimits.end();++it){
    tmp=*it;
    nv=rotateVec(tmp.dirc,ctl.rote);
    napt=ctl.tran+ctl.fixpt+rotateVec(tmp.apt-ctl.fixpt,ctl.rote);
    (*it).dirc=nv;
    (*it).apt=napt;
  }
};
 
template<class T> class plnBoundary:public Boundary<T>{
 public:
  Vec normal;  // normal force acting on the boundary by all contacting particles 
  Vec tangt;   // tangential force acting on the boundary
  Vec moment;  // moment on the boundary
  std::vector<T*> possBdryParticle; // possible boundary particles of this specific boundary
 public:
 plnBoundary(std::ifstream &ifs):Boundary<T>(ifs){
    normal=0;
    tangt=0;
    moment=0;
    this->avgNormal=0;
    this->avgPenetr=0;
    this->contactNum=0;
  };
  int getBdryID() {return this->boundaryId;}
  void display() const;
  REAL distToBdry(Vec posi) const;
  void findParticleOnBoundary(std::vector<T*>& ptcls);
  Vec getApt() const;
  Vec getDirc() const;
  plnBoundary<T>* getBdry(int bdryid) const{
    return this;
  }
  void boundaryForce(std::map<int,std::vector<BoundaryTgt> >& boundaryTgtMap);
  Vec getNormalForce() const{return normal;}
  REAL getAvgNormal() const{return this->avgNormal;}
  REAL getAvgPenetr() const{return this->avgPenetr;}
  int getCntnum() const{return this->contactNum;}
  Vec getShearForce() const{return tangt;}
};

template<class T>
Vec plnBoundary<T>::getApt() const{
  return (*this->CoefOfLimits.begin()).apt;
};

template<class T>
Vec plnBoundary<T>::getDirc() const{
  return (*this->CoefOfLimits.begin()).dirc;
};

template<class T>
void plnBoundary<T>::display() const{
  Boundary<T>::display();
  std::cout << "normal: " << normal.getX() << " " << normal.getY() << " " << normal.getZ() << std::endl;
  typename std::vector<T*>::const_iterator it;
  int i=0;
  for(it=possBdryParticle.begin();it!=possBdryParticle.end();++it){
    if(i++<10)
      std::cout << (*it)->getId();
    else{
      i=0;
      std::cout << (*it)->getId() << std::endl;
    }
  }
};
 
template<class T>
REAL plnBoundary<T>::distToBdry(Vec posi) const{
  Vec dv=(*this->CoefOfLimits.begin()).dirc;
  Vec pt=(*this->CoefOfLimits.begin()).apt;
  Vec ndv=normalize(dv);
  return (posi-pt)%ndv;
};

template<class T>
void plnBoundary<T>::findParticleOnBoundary(std::vector<T*>& ptcls){
  typename std::vector<T*>::iterator it;
  std::vector<BdryCoef>::iterator bt;
  bool next;
  possBdryParticle.clear();
  REAL dist, r;
  Vec posi, ndirc;
  for (it=ptcls.begin();it!=ptcls.end();++it){
    if ( (*it)->getType() == 0 ) { // only process free particles, excluding type 5
      posi=(*it)->getCurrPos();
      dist=distToBdry(posi);
      if(dist>=0 || fabs(dist) > (*it)->getA()) // outside to CoefOfLimits[0] or inside too much
	continue;
      /*
      debugInf << "boundary.h: iter=" << iteration
	       << " bdryId=" << getBdryID()
	       << " ptclId=" << (*it)->getId() << std::endl;
      */
      next=true;
      for (bt=++this->CoefOfLimits.begin();bt!=this->CoefOfLimits.end();++bt){ // CoefOfLimits[1,2,...]
	ndirc=normalize((*bt).dirc);
	r=vfabs((posi-(*bt).apt)-(posi-(*bt).apt)%ndirc*ndirc);
	if( ( (*bt).order==1 && (posi-(*bt).apt)%(*bt).dirc >= 0 ) ||
	    ( (*bt).order==2 && (r-(*bt).rad)*(*bt).side<0 ) ){
	  next=false; // the particle is out of boundary, process next particle
	  break;
	}
      }
      if(next)
	possBdryParticle.push_back(*it);
    }
  }
  
  /*
  if (getBdryID() > 6) { // cavity boundaries
    debugInf << "boundary.h: iter=" << iteration
	     << " bdryId=" << getBdryID()
	     << " ptcl#=" << possBdryParticle.size();
    for(it = possBdryParticle.begin(); it != possBdryParticle.end(); ++it)
      debugInf << " " << (*it)->getId();
    debugInf << std::endl;
  }
  */
};

template<class T>
void plnBoundary<T>::boundaryForce(std::map<int,std::vector<BoundaryTgt> >& boundaryTgtMap){
  typename std::vector<T*>::iterator it;
  this->avgNormal=0;
  this->avgPenetr=0;
  this->contactNum=0;
  normal=0;
  tangt=0;
  moment=0;
  
  // for each plane boundary, define a temparory variable vtmp to use,
  // better than define a member variable which needs to be cleared.
  // and vtmp is initialized as empty in each iteration.
  std::vector<BoundaryTgt> vtmp;
  
  // for each possible boundary particle
  REAL penetr=0;
  int count=0;
  for (it=possBdryParticle.begin();it!=possBdryParticle.end();++it){
    penetr=0;
    (*it)->planeRBForce(this,boundaryTgtMap,vtmp,penetr);
    this->avgPenetr += penetr;
    count++;
  }
  if (count>0) this->avgPenetr /= count;
  this->contactNum=count;
  
  // checkout tangential forces and displacements after each particle is processed
  boundaryTgtMap[this->boundaryId]=vtmp;
};

template<class T> class cylBoundary:public Boundary<T>{
 public:
  Vec normal; 
  std::vector<T*> possBdryParticle;
 public:
 cylBoundary(std::ifstream &ifs):Boundary<T>(ifs){normal=0;}
  void display() const;
  REAL distToBdry(Vec posi) const;
  void findParticleOnBoundary(std::vector<T*>& ptcls);
  void boundaryForce();
  Vec getNormalForce() const{return normal;};
};

template<class T>
void cylBoundary<T>::display() const{
  Boundary<T>::display();
  std::cout << "normal: " << normal.getX() << " " << normal.getY() << " " << normal.getZ() << std::endl;
  typename std::vector<T*>::const_iterator it;
  int i=0;
  for(it=possBdryParticle.begin();it!=possBdryParticle.end();++it){
    if(i++<10)
      std::cout << (*it)->getId();
    else{
      i=0;
      std::cout << (*it)->getId() << std::endl;
    }
  }
};

template<class T>
REAL cylBoundary<T>::distToBdry(Vec posi) const{
  Vec ndc=(*this->CoefOfLimits.begin()).dirc;
  Vec napt=(*this->CoefOfLimits.begin()).apt;
  REAL r=(*this->CoefOfLimits.begin()).rad;
  Vec norm=normalize(ndc);
  return fabs(r-vfabs((posi-napt)-(posi-napt)%norm*norm));
};

template<class T>
void cylBoundary<T>::findParticleOnBoundary(std::vector<T*> &ptcls){
  typename std::vector<T*>::iterator it;
  std::vector<BdryCoef>::iterator bt;
  bool next;
  possBdryParticle.clear();
  REAL dist,r;
  Vec posi, ndirc;
  for (it=ptcls.begin();it!=ptcls.end();++it){
    posi=(*it)->getCurrPos();
    dist=distToBdry(posi);
    if (dist>(*it)->getA())
      continue;
    next=false;
    for (bt=++this->CoefOfLimits.begin();bt!=this->CoefOfLimits.end();++bt){
      ndirc=normalize((*bt).dirc);
      r=vfabs((posi-(*bt).apt)-(posi-(*bt).apt)%ndirc*ndirc);
      if( ( (*bt).order==1&&(posi-(*bt).apt)%(*bt).dirc>(*it)->getA() )||
	  ( (*bt).order==2&&(r-(*bt).rad)*(*bt).side<0 ) ){
	next=true;//the particle is outof boundary, process next particle
	break;
      }
    }
    if(!next)
      possBdryParticle.push_back(*it);
  }
};

template<class T>
void cylBoundary<T>::boundaryForce(){
  Cylinder cyl;
  typename std::vector<T*>::iterator it;
  BdryCoef tmp;
  tmp=*this->CoefOfLimits.begin();
  cyl.setRadius(tmp.rad);
  cyl.setCenter(tmp.apt);
  normal=0;
  for (it=possBdryParticle.begin();it!=possBdryParticle.end();++it){
    normal-=(*it)->cylinderRBForce(this->boundaryId,cyl,tmp.side);
  }
};

} // namespace dem ends

#endif
