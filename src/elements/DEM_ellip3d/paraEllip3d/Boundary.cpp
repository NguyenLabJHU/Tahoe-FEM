#include "Boundary.h"
#include "Particle.h"
// use both pointer to and variable of class Particle
// separate interface and implemenation completely to get rid of dependence

namespace dem {

  BoundaryTgt::BoundaryTgt()
    :particleId(0), tgtForce(0), tgtDisp(0), tgtLoading(false), tgtDispStart(0), tgtPeak(0), tgtSlide(false), tgtRoll(false)
  {}
    
  BoundaryTgt::BoundaryTgt(std::size_t _particleId, Vec _v1, Vec _v2, bool _b, Vec _v3, REAL _tp, bool _s, bool _r)
    :particleId(_particleId), tgtForce(_v1), tgtDisp(_v2), tgtLoading(_b), tgtDispStart(_v3), tgtPeak(_tp), tgtSlide(_s), tgtRoll(_r)
  {}
  
  BdryContact::BdryContact()
    :ptcl(NULL), point(0), normal(0), tangt(0), penetr(0) 
  {}
    
  BdryContact::BdryContact(Particle *p, Vec pt, Vec nm, Vec tg, REAL pntr)
    :ptcl(p), point(pt), normal(nm), tangt(tg), penetr(pntr) 
  {}

  BdryContact::~BdryContact() {
    //delete ptcl; // it causes double free or corruption.
  }

  void BdryContact::print(std::ostream &os) {
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

  Plane::Plane()
    :direc(0), point(0) 
  {}
    
  Plane::Plane(Vec dir, Vec pt)
    :direc(dir), point(pt) 
  {}
    
  void Plane::print(std::ostream &os) {
    os << std::setw(OWID) << direc.getX()
       << std::setw(OWID) << direc.getY()
       << std::setw(OWID) << direc.getZ()
       << std::setw(OWID) << point.getX()
       << std::setw(OWID) << point.getY()
       << std::setw(OWID) << point.getZ() << std::endl;
  }

  Vec Plane::getDirec() const { return direc; }
  Vec Plane::getPoint() const { return point; }

  Boundary::Boundary(std::size_t i, std::size_t tp, std::size_t en)
    :id(i), type(tp), extraNum(en), contactNum(0), normal(0), tangt(0), penetr(0)  
  {}

  Boundary::Boundary(std::size_t tp, std::ifstream &ifs) {
    type = tp;
    ifs >> extraNum;
    ifs >> id;
  }
  
  std::size_t Boundary::getId() { return id; }
  std::size_t Boundary::getType() { return type; }
  std::vector<Particle *> & Boundary::getPossParticle() {return possParticle;}
  std::vector<BdryContact> & Boundary::getContactInfo() {return contactInfo;}
  std::size_t Boundary::getContactNum() const { return contactNum; }
  Vec  Boundary::getNormalForce() const { return normal; }
  Vec  Boundary::getTangtForce() const { return tangt; }
  REAL Boundary::getAvgPenetr() const { return penetr; }

  void Boundary::print(std::ostream &os) {
    os << std::endl 
       << std::setw(OWID) << type
       << std::setw(OWID) << extraNum << std::endl
       << std::setw(OWID) << id;
  }
  
  void Boundary::printContactInfo(std::ostream &os) {
    os << std::setw(OWID) << id << std::endl
       << std::setw(OWID) << contactInfo.size() << std::endl
       << std::setw(OWID) << "pos_x"
       << std::setw(OWID) << "pos_y"
       << std::setw(OWID) << "pos_z"
       << std::setw(OWID) << "normal_x"
       << std::setw(OWID) << "normal_y"
       << std::setw(OWID) << "normal_z"
       << std::setw(OWID) << "tangt_x"
       << std::setw(OWID) << "tangt_y"
       << std::setw(OWID) << "tangt_z"
       << std::setw(OWID) << "pentr" << std::endl;

    for (std::vector<BdryContact>::iterator it = contactInfo.begin(); it != contactInfo.end(); ++it)
      it->print(os);
  }
  
  void Boundary::clearStatForce() {
    contactNum = 0;
    normal = 0;
    tangt  = 0;
    penetr = 0;
  }

  void Boundary::updateStatForce() {
    clearStatForce();
    contactNum = contactInfo.size();
    for (std::vector<BdryContact>::iterator it = contactInfo.begin(); it != contactInfo.end(); ++it) {
      normal += it->normal;
      tangt  += it->tangt;
      penetr += it->penetr;
    }
    if (contactNum != 0) 
      penetr /= contactNum;
  }

  void Boundary::clearContactInfo() {
    possParticle.clear();
    contactInfo.clear();
  }

  planeBoundary::planeBoundary(std::size_t i, std::size_t tp, std::size_t en)
    :Boundary(i, tp, en), direc(0), point(0), prevPoint(0), veloc(0), prevVeloc(0)
  {}

  planeBoundary::planeBoundary(std::size_t tp, std::ifstream &ifs)
    :Boundary(tp, ifs) {
    REAL dx, dy, dz, px, py, pz;
    ifs >> dx >> dy >> dz >> px >> py >> pz;
    direc = Vec(dx, dy, dz);
    point = Vec(px, py, pz);
    prevPoint = point;
    prevVeloc = veloc = 0;
    for (std::size_t i = 0; i < extraNum; ++i) {
      ifs >> dx >> dy >> dz >> px >> py >> pz;
      extraEdge.push_back(Plane(Vec(dx, dy, dz), Vec(px, py, pz)));
    }
  }

  Vec planeBoundary::getDirec() const { return direc; }
  Vec planeBoundary::getPoint() const { return point; }
  Vec planeBoundary::getVeloc() const { return veloc; }
  Vec planeBoundary::getPrevPoint() const { return prevPoint; }
  Vec planeBoundary::getPrevVeloc() const { return prevVeloc; }

  void planeBoundary::setDirec(Vec dir) { direc = dir; }
  void planeBoundary::setPoint(Vec pnt) { point = pnt; }
  void planeBoundary::setVeloc(Vec vel) { veloc = vel; }

  REAL planeBoundary::distanceToBdry(Vec pos) const { return (pos - point) * normalize(direc); }
  REAL planeBoundary::distanceToBdry(Vec pos, Plane pn) const { return (pos - pn.getPoint()) * normalize(pn.getDirec()); }

  void planeBoundary::print(std::ostream &os) {
    Boundary::print(os);
    os << std::setw(OWID) << direc.getX()
       << std::setw(OWID) << direc.getY()
       << std::setw(OWID) << direc.getZ()
       << std::setw(OWID) << point.getX()
       << std::setw(OWID) << point.getY()
       << std::setw(OWID) << point.getZ() << std::endl;

    for (std::vector<Plane>::iterator et = extraEdge.begin(); et != extraEdge.end(); ++et)
      os << std::setw(OWID) << " "
	 << std::setw(OWID) << et->getDirec().getX()
	 << std::setw(OWID) << et->getDirec().getY()
	 << std::setw(OWID) << et->getDirec().getZ()
	 << std::setw(OWID) << et->getPoint().getX()
	 << std::setw(OWID) << et->getPoint().getY()
	 << std::setw(OWID) << et->getPoint().getZ()
	 << std::endl;
  }
  
  void planeBoundary::printContactInfo(std::ostream &os) {
    Boundary::printContactInfo(os);
    os << std::setw(OWID) << " "
       << std::setw(OWID) << " "
       << std::setw(OWID) << " "
       << std::setw(OWID) << normal.getX()
       << std::setw(OWID) << normal.getY()
       << std::setw(OWID) << normal.getZ() 
       << std::setw(OWID) << tangt.getX()
       << std::setw(OWID) << tangt.getY()
       << std::setw(OWID) << tangt.getZ() 
       << std::setw(OWID) << penetr << std::endl << std::endl;;
  }

  void planeBoundary::findBdryContact(std::vector<Particle *> &ptcls) {
    possParticle.clear();
    contactInfo.clear();
    clearStatForce();

    for (std::vector<Particle *>::iterator it = ptcls.begin(); it != ptcls.end(); ++it) {
      if ( (*it)->getType() == 0 ) { // only process free particles, excluding type 5
	REAL dist = distanceToBdry((*it)->getCurrPos());
	if (dist < 0 && fabs(dist) <= (*it)->getA()) {
	  bool inside = true;
	  for (std::vector<Plane>::iterator et = extraEdge.begin(); et != extraEdge.end(); ++et) {
	    REAL eDist = distanceToBdry((*it)->getCurrPos(), (*et));
	    if (eDist >= 0 ) {
	      inside = false;
	      break;
	    }
	  }
	  if (inside)
	    possParticle.push_back(*it);
	}
      }
    }
  }
  
  void planeBoundary::boundaryForce(std::map<std::size_t,std::vector<BoundaryTgt> > &boundaryTgtMap) {
    // for each plane boundary, define a temparory variable vtmp to use,
    // better than define a member variable which needs to be cleared.
    // and vtmp is initialized as empty in each iteration.
    std::vector<BoundaryTgt> vtmp;
    
    // for each possible boundary particle
    for (std::vector<Particle *>::iterator it = possParticle.begin(); it != possParticle.end(); ++it)
      (*it)->planeRBForce(this, boundaryTgtMap, vtmp);
    
    // checkout tangential forces and displacements after each particle is processed
    boundaryTgtMap[this->id] = vtmp;

    updateStatForce();
  }

  void planeBoundary::updateIsotropic(REAL sigma, REAL areaX, REAL areaY, REAL areaZ) {

    REAL forceDamp = dem::Parameter::getSingleton().parameter["forceDamp"];
    REAL massScale = dem::Parameter::getSingleton().parameter["massScale"];
    //REAL mass = dem::Parameter::getSingleton().parameter["boundaryMass"];
    REAL boundaryRate = dem::Parameter::getSingleton().parameter["boundaryRate"];
    REAL topSpeedup = dem::Parameter::getSingleton().parameter["topSpeedup"];
    REAL tol = dem::Parameter::getSingleton().parameter["tractionErrorTol"];
    REAL atf = forceDamp*timeStep;

    REAL vel, pos;
    switch (id) {
    case 1: 
      if (fabs(normal.getX()/areaX + sigma)/sigma > tol) {
	vel = ((normal.getX() + sigma*areaX)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getX() * (2-atf) / (2+atf) + (normal.getX() + sigma * areaX) / mass * timeStep * 2 / (2 + atf);
	pos = prevPoint.getX() + vel * timeStep;
	setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ() ));
	setPoint(Vec(pos, getPoint().getY(), getPoint().getZ() )); }
      break;
    case 2:
      if (fabs(normal.getX()/areaX - sigma)/sigma > tol) {
	vel = ((normal.getX() - sigma*areaX)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getX() * (2-atf) / (2+atf) + (normal.getX() - sigma * areaX) / mass * timeStep * 2 / (2 + atf);
	pos = prevPoint.getX() + vel * timeStep;
	setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ() ));
	setPoint(Vec(pos, getPoint().getY(), getPoint().getZ() )); }
      break;
    case 3:
      if (fabs(normal.getY()/areaY + sigma)/sigma > tol) {
	vel = ((normal.getY() + sigma*areaY)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() + sigma * areaY) / mass * timeStep * 2 / (2 + atf);
	pos = prevPoint.getY() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ() ));
	setPoint(Vec(getPoint().getX(), pos, getPoint().getZ() )); }
      break;
    case 4:
      if (fabs(normal.getY()/areaY - sigma)/sigma > tol) {
	vel = ((normal.getY() - sigma*areaY)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() - sigma * areaY) / mass * timeStep * 2 / (2 + atf);
	pos = prevPoint.getY() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ() ));
	setPoint(Vec(getPoint().getX(), pos, getPoint().getZ() )); }
      break;
    case 5:
      if (fabs(normal.getZ()/areaZ + sigma)/sigma > tol) {
	vel = ((normal.getZ() + sigma*areaZ)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() + sigma * areaZ) / mass * timeStep * 2 / (2 + atf);
	pos = prevPoint.getZ() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel ));
	setPoint(Vec(getPoint().getX(), getPoint().getY(), pos )); }
      break;
    case 6:
      if (fabs(normal.getZ()/areaZ - sigma)/sigma > tol) {
	vel = ((normal.getZ() - sigma*areaZ)>0 ? 1:-1) * boundaryRate;
	if (normal.getZ() == 0 ) vel = -boundaryRate*topSpeedup;
	//vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() - sigma * areaZ) / mass * timeStep * 2 / (2 + atf);
	pos = prevPoint.getZ() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel ));
	setPoint(Vec(getPoint().getX(), getPoint().getY(), pos )); }
      break;
    }
    prevPoint = point;
    prevVeloc = veloc;
  }

  void planeBoundary::updateOedometer(REAL sigma, REAL areaX, REAL areaY, REAL areaZ) {

    REAL forceDamp = dem::Parameter::getSingleton().parameter["forceDamp"];
    REAL massScale = dem::Parameter::getSingleton().parameter["massScale"];
    //REAL mass = dem::Parameter::getSingleton().parameter["boundaryMass"];
    REAL boundaryRate = dem::Parameter::getSingleton().parameter["boundaryRate"];
    REAL tol = dem::Parameter::getSingleton().parameter["tractionErrorTol"];
    REAL atf = forceDamp*timeStep;

    REAL vel, pos;
    switch (id) {
    case 5:
      if (fabs(normal.getZ()/areaZ + sigma)/sigma > tol) {
	vel = ((normal.getZ() + sigma*areaZ)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() + sigma * areaZ) / mass * timeStep * 2 / (2 + atf);
	pos = prevPoint.getZ() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel ));
	setPoint(Vec(getPoint().getX(), getPoint().getY(), pos )); }
      break;
    case 6:
      if (fabs(normal.getZ()/areaZ - sigma)/sigma > tol) {
	vel = ((normal.getZ() - sigma*areaZ)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() - sigma * areaZ) / mass * timeStep * 2 / (2 + atf);
	pos = prevPoint.getZ() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel ));
	setPoint(Vec(getPoint().getX(), getPoint().getY(), pos )); }
      break;
    }
    prevPoint = point;
    prevVeloc = veloc;
  }

  void planeBoundary::updateTriaxial(REAL sigma, REAL areaX, REAL areaY, REAL areaZ) {
    std::size_t triaxialType = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["triaxialType"]);
    std::size_t unloadStep = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["unloadStep"]);

    REAL forceDamp = dem::Parameter::getSingleton().parameter["forceDamp"];
    REAL massScale = dem::Parameter::getSingleton().parameter["massScale"];
    //REAL mass = dem::Parameter::getSingleton().parameter["boundaryMass"];
    REAL boundaryRate = dem::Parameter::getSingleton().parameter["boundaryRate"];
    REAL tol = dem::Parameter::getSingleton().parameter["tractionErrorTol"];
    REAL atf = forceDamp*timeStep;

    REAL vel, pos;
    switch (id) {
    case 1: 
      if (fabs(normal.getX()/areaX + sigma)/sigma > tol) {
	vel = ((normal.getX() + sigma*areaX)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getX() * (2-atf) / (2+atf) + (normal.getX() + sigma * areaX) / mass * timeStep * 2 / (2 + atf);
	pos = prevPoint.getX() + vel * timeStep;
	setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ() ));
	setPoint(Vec(pos, getPoint().getY(), getPoint().getZ() )); }
      break;
    case 2:
      if (fabs(normal.getX()/areaX - sigma)/sigma > tol) {
	vel = ((normal.getX() - sigma*areaX)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getX() * (2-atf) / (2+atf) + (normal.getX() - sigma * areaX) / mass * timeStep * 2 / (2 + atf);
	pos = prevPoint.getX() + vel * timeStep;
	setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ() ));
	setPoint(Vec(pos, getPoint().getY(), getPoint().getZ() )); }
      break;
    case 3:
      if (fabs(normal.getY()/areaY + sigma)/sigma > tol) {
	vel = ((normal.getY() + sigma*areaY)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() + sigma * areaY) / mass * timeStep * 2 / (2 + atf);
	pos = prevPoint.getY() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ() ));
	setPoint(Vec(getPoint().getX(), pos, getPoint().getZ() )); }
      break;
    case 4:
      if (fabs(normal.getY()/areaY - sigma)/sigma > tol) {
	vel = ((normal.getY() - sigma*areaY)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() - sigma * areaY) / mass * timeStep * 2 / (2 + atf);
	pos = prevPoint.getY() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ() ));
	setPoint(Vec(getPoint().getX(), pos, getPoint().getZ() )); }
      break;
    case 5:
      if (triaxialType == 1)
	vel = boundaryRate;
      else if (triaxialType == 2) {
	if (iteration <= unloadStep) // loading
	  vel = boundaryRate;
	else if (iteration > unloadStep && fabs(normal.getZ()/areaZ) >= 1.5*sigma && iteration <= 1.5*unloadStep) // unloading
	  vel = -boundaryRate;
	else if (iteration > 1.5 *unloadStep) // reloading. Note there are invalid loops if stress < and iteration <, it is OK.
	  vel = boundaryRate;
      }
      //vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() + sigma * areaZ) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getZ() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel ));
      setPoint(Vec(getPoint().getX(), getPoint().getY(), pos ));
      break;
    case 6:
      if (triaxialType == 1)
	vel = -boundaryRate;
      else if (triaxialType == 2) {
	if (iteration <= unloadStep) // loading
	  vel = -boundaryRate;
	else if (iteration > unloadStep && fabs(normal.getZ()/areaZ) >= 1.5*sigma && iteration <= 1.5*unloadStep) // unloading
	  vel = boundaryRate;
	else if (iteration > 1.5 *unloadStep) // reloading. Note there are invalid loops if stress < and iteration <, it is OK.
	  vel = -boundaryRate;
      }
      //vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() - sigma * areaZ) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getZ() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel ));
      setPoint(Vec(getPoint().getX(), getPoint().getY(), pos ));
      break;
    }
    prevPoint = point;
    prevVeloc = veloc;
  }

  void planeBoundary::updatePlaneStrain(REAL sigma, REAL areaX, REAL areaY, REAL areaZ) {
    std::size_t plnstrnType = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["plnstrnType"]);
    std::size_t unloadStep = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["unloadStep"]);

    REAL forceDamp = dem::Parameter::getSingleton().parameter["forceDamp"];
    REAL massScale = dem::Parameter::getSingleton().parameter["massScale"];
    //REAL mass = dem::Parameter::getSingleton().parameter["boundaryMass"];
    REAL boundaryRate  = dem::Parameter::getSingleton().parameter["boundaryRate"];
    REAL sideRateRatio = dem::Parameter::getSingleton().parameter["sideRateRatio"];
    REAL tol = dem::Parameter::getSingleton().parameter["tractionErrorTol"];
    REAL atf = forceDamp*timeStep;

    REAL vel, pos;
    switch (id) { // boundary x1(1) and boundary x2(2) do not move
    case 3:
      if (fabs(normal.getY()/areaY + sigma)/sigma > tol) {
	vel = ((normal.getY() + sigma*areaY)>0 ? 1:-1) * boundaryRate * sideRateRatio;
	//vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() + sigma * areaY) / mass * timeStep * 2 / (2 + atf);
	pos = prevPoint.getY() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ() ));
	setPoint(Vec(getPoint().getX(), pos, getPoint().getZ() )); }
      break;
    case 4:
      if (fabs(normal.getY()/areaY - sigma)/sigma > tol) {
	vel = ((normal.getY() - sigma*areaY)>0 ? 1:-1) * boundaryRate * sideRateRatio;
	//vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() - sigma * areaY) / mass * timeStep * 2 / (2 + atf);
	pos = prevPoint.getY() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ() ));
	setPoint(Vec(getPoint().getX(), pos, getPoint().getZ() )); }
      break;
    // displacement control, leading to zero volumetric strain
    case 5:
      if (plnstrnType == 1)
	vel = boundaryRate;
      else if (plnstrnType == 2) {
	if (iteration <= unloadStep) // loading
	  vel = boundaryRate;
	else if (iteration > unloadStep && fabs(normal.getZ()/areaZ) >= 1.5*sigma && iteration <= 1.5*unloadStep) // unloading
	  vel = -boundaryRate;
	else if (iteration > 1.5 *unloadStep) // reloading. Note there are invalid loops if stress < and iteration <, it is OK.
	  vel = boundaryRate;
      }
      //vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() + sigma * areaZ) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getZ() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel ));
      setPoint(Vec(getPoint().getX(), getPoint().getY(), pos ));
      break;
    case 6:
      if (plnstrnType == 1)
	vel = -boundaryRate;
      else if (plnstrnType == 2) {
	if (iteration <= unloadStep) // loading
	  vel = -boundaryRate;
	else if (iteration > unloadStep && fabs(normal.getZ()/areaZ) >= 1.5*sigma && iteration <= 1.5*unloadStep) // unloading
	  vel = boundaryRate;
	else if (iteration > 1.5 *unloadStep) // reloading. Note there are invalid loops if stress < and iteration <, it is OK.
	  vel = -boundaryRate;
      }
      //vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() - sigma * areaZ) / mass * timeStep * 2 / (2 + atf);
      pos = prevPoint.getZ() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel ));
      setPoint(Vec(getPoint().getX(), getPoint().getY(), pos ));
      break;
    }
    prevPoint = point;
    prevVeloc = veloc;
  }

  void planeBoundary::updateTrueTriaxial(REAL sigma, REAL areaX, REAL areaY, REAL areaZ, REAL sigmaX, REAL sigmaY) {
    // sigma implies sigmaZ

    REAL forceDamp = dem::Parameter::getSingleton().parameter["forceDamp"];
    REAL massScale = dem::Parameter::getSingleton().parameter["massScale"];
    //REAL mass = dem::Parameter::getSingleton().parameter["boundaryMass"];
    REAL boundaryRate = dem::Parameter::getSingleton().parameter["boundaryRate"];
    REAL tol = dem::Parameter::getSingleton().parameter["tractionErrorTol"];
    REAL atf = forceDamp*timeStep;

    REAL vel, pos;
    switch (id) {
    case 1: 
      if (fabs(normal.getX()/areaX + sigmaX)/sigmaX > tol) {
	vel = ((normal.getX() + sigmaX*areaX)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getX() * (2-atf) / (2+atf) + (normal.getX() + sigma * areaX) / mass * timeStep * 2 / (2 + atf);
	pos = prevPoint.getX() + vel * timeStep;
	setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ() ));
	setPoint(Vec(pos, getPoint().getY(), getPoint().getZ() )); }
      break;
    case 2:
      if (fabs(normal.getX()/areaX - sigmaX)/sigmaX > tol) {
	vel = ((normal.getX() - sigmaX*areaX)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getX() * (2-atf) / (2+atf) + (normal.getX() - sigma * areaX) / mass * timeStep * 2 / (2 + atf);
	pos = prevPoint.getX() + vel * timeStep;
	setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ() ));
	setPoint(Vec(pos, getPoint().getY(), getPoint().getZ() )); }
      break;
    case 3:
      if (fabs(normal.getY()/areaY + sigmaY)/sigmaY > tol) {
	vel = ((normal.getY() + sigmaY*areaY)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() + sigma * areaY) / mass * timeStep * 2 / (2 + atf);
	pos = prevPoint.getY() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ() ));
	setPoint(Vec(getPoint().getX(), pos, getPoint().getZ() )); }
      break;
    case 4:
      if (fabs(normal.getY()/areaY - sigmaY)/sigmaY > tol) {
	vel = ((normal.getY() - sigmaY*areaY)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getY() * (2-atf) / (2+atf) + (normal.getY() - sigma * areaY) / mass * timeStep * 2 / (2 + atf);
	pos = prevPoint.getY() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ() ));
	setPoint(Vec(getPoint().getX(), pos, getPoint().getZ() )); }
      break;
    case 5:
      if (fabs(normal.getZ()/areaZ + sigma)/sigma > tol) {
	vel = ((normal.getZ() + sigma*areaZ)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() + sigma * areaZ) / mass * timeStep * 2 / (2 + atf);
	pos = prevPoint.getZ() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel ));
	setPoint(Vec(getPoint().getX(), getPoint().getY(), pos )); }
      break;
    case 6:
      if (fabs(normal.getZ()/areaZ - sigma)/sigma > tol) {
	vel = ((normal.getZ() - sigma*areaZ)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getZ() * (2-atf) / (2+atf) + (normal.getZ() - sigma * areaZ) / mass * timeStep * 2 / (2 + atf);
	pos = prevPoint.getZ() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel ));
	setPoint(Vec(getPoint().getX(), getPoint().getY(), pos )); }
      break;
    }
    prevPoint = point;
    prevVeloc = veloc;
  }

  cylinderBoundary::cylinderBoundary()
    :Boundary(), direc(0), point(0), prevPoint(0), veloc(0), prevVeloc(0), radius(0)
  {}

  cylinderBoundary::cylinderBoundary(std::size_t tp, std::ifstream &ifs)
    :Boundary(tp, ifs) {
    REAL dx, dy, dz, px, py, pz;
    ifs >> dx >> dy >> dz >> px >> py >> pz >> radius;
    direc = Vec(dx, dy, dz);
    point = Vec(px, py, pz);
  }

  Vec cylinderBoundary::getDirec() const { return direc; }
  Vec cylinderBoundary::getPoint() const { return point; }
  Vec cylinderBoundary::getVeloc() const { return veloc; }
  Vec cylinderBoundary::getPrevPoint() const { return prevPoint; }
  Vec cylinderBoundary::getPrevVeloc() const { return prevVeloc; }
  REAL cylinderBoundary::getRadius() const { return radius; }

  void cylinderBoundary::setDirec(Vec dir) { direc = dir; }
  void cylinderBoundary::setPoint(Vec pnt) { point = pnt; }
  void cylinderBoundary::setVeloc(Vec vel) { veloc = vel; }
 
  void cylinderBoundary::print(std::ostream &os) {
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

  void cylinderBoundary::findBdryContact(std::vector<Particle *> &ptcls) {
    possParticle.clear();
    contactInfo.clear();

    for (std::vector<Particle *>::iterator it = ptcls.begin(); it != ptcls.end(); ++it) {
      if ( (*it)->getType() == 0 ) { // only process free particles, excluding type 5
	;
      }
    }
  }
  
  void cylinderBoundary::boundaryForce(std::map<std::size_t,std::vector<BoundaryTgt> > &boundaryTgtMap) {
    // for each plane boundary, define a temparory variable vtmp to use,
    // better than define a member variable which needs to be cleared.
    // and vtmp is initialized as empty in each iteration.
    std::vector<BoundaryTgt> vtmp;
    
    // for each possible boundary particle
    for (std::vector<Particle *>::iterator it = possParticle.begin(); it != possParticle.end(); ++it)
      ;// (*it)->cylinderRBForce();
    
    // checkout tangential forces and displacements after each particle is processed
    boundaryTgtMap[this->id] = vtmp;

    updateStatForce();
  }

} // namespace dem ends

