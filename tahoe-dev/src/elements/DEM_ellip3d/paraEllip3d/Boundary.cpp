///////////////////////////////////////////////////////////////////////////////////////////////////////
//                                   Code: ParaEllip3d-CFD                                           //
//                                 Author: Dr. Beichuan Yan                                          //
//                                  Email: beichuan.yan@colorado.edu                                 //
//                              Institute: University of Colorado Boulder                            //
///////////////////////////////////////////////////////////////////////////////////////////////////////

#include "Boundary.h"
#include "Particle.h"
// use both pointer to and variable of class Particle

namespace dem {
  
  Boundary::Boundary(std::size_t tp, std::ifstream &ifs) {
    type = tp;
    ifs >> extraNum;
    ifs >> id;
  }
  
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
    moment = 0;
    penetr = 0;
  }

  void Boundary::updateStatForce() {
    clearStatForce();
    contactNum = contactInfo.size();
    for (std::vector<BdryContact>::iterator it = contactInfo.begin(); it != contactInfo.end(); ++it) {
      normal += it->normal;
      tangt  += it->tangt;
      penetr += it->penetr;
      double d1 = it->point.getX();
      double d2 = it->point.getY();
      double d3 = it->point.getZ();
      double f1 = (it->normal + it->tangt).getX();
      double f2 = (it->normal + it->tangt).getY();
      double f3 = (it->normal + it->tangt).getZ();
      moment += Vec(d2*f3 - d3*f2, d3*f1 - d1*f3, d1*f2 - d2*f1);
    }
    if (contactNum != 0) 
      penetr /= contactNum;
  }

  void Boundary::clearContactInfo() {
    possParticle.clear();
    contactInfo.clear();
  }

  void Boundary::clearPossParticle() {
    possParticle.clear();
  }

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
    // for each plane boundary, define a temporary variable vtmp to use,
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

  void planeBoundary::updateMoveWall(REAL sigma, REAL areaX, REAL areaY, REAL areaZ) {
    REAL bdry1Rate = dem::Parameter::getSingleton().parameter["bdry1Rate"];
    REAL bdry2Rate = dem::Parameter::getSingleton().parameter["bdry2Rate"];
    REAL bdry3Rate = dem::Parameter::getSingleton().parameter["bdry3Rate"];
    REAL bdry4Rate = dem::Parameter::getSingleton().parameter["bdry4Rate"];
    REAL bdry6Rate = dem::Parameter::getSingleton().parameter["bdry6Rate"];

    REAL vel, pos;
    switch (id) {
    case 1: 
      vel = - bdry1Rate;
      pos = prevPoint.getX() + vel * timeStep;
      setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ() ));
      setPoint(Vec(pos, getPoint().getY(), getPoint().getZ() ));
      break;
    case 2:
      vel = bdry2Rate;
      pos = prevPoint.getX() + vel * timeStep;
      setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ() ));
      setPoint(Vec(pos, getPoint().getY(), getPoint().getZ() ));
      break;
    case 3:
      vel = - bdry3Rate;
      pos = prevPoint.getY() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ() ));
      setPoint(Vec(getPoint().getX(), pos, getPoint().getZ() ));
      break;
    case 4:
      vel = bdry4Rate;
      pos = prevPoint.getY() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ() ));
      setPoint(Vec(getPoint().getX(), pos, getPoint().getZ() ));
      break;
    case 6:
      vel = bdry6Rate;
      pos = prevPoint.getZ() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel ));
      setPoint(Vec(getPoint().getX(), getPoint().getY(), pos ));
      break;
    }
    prevPoint = point;
    prevVeloc = veloc;
  }

  void planeBoundary::updateIsotropic(REAL sigma, REAL areaX, REAL areaY, REAL areaZ) {
    std::size_t isotropicType = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["isotropicType"]);
    REAL forceDamp = dem::Parameter::getSingleton().parameter["forceDamp"];
    REAL massScale = dem::Parameter::getSingleton().parameter["massScale"];
    //REAL mass = dem::Parameter::getSingleton().parameter["boundaryMass"];
    REAL boundaryRate = dem::Parameter::getSingleton().parameter["boundaryRate"];
    REAL tol = dem::Parameter::getSingleton().parameter["tractionErrorTol"];
    REAL atf = forceDamp;

    // topSpeedup does not exist when isotropicType = 2 or 3
    REAL topSpeedup = dem::Parameter::getSingleton().parameter["topSpeedup"];

    REAL vel, pos;
    switch (id) {
    case 1: 
      if (fabs(normal.getX()/areaX + sigma)/sigma > tol) {
	vel = ((normal.getX() + sigma*areaX)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getX() * (1-atf) / (1+atf) + (normal.getX() + sigma * areaX) / mass * timeStep / (1+atf);
	pos = prevPoint.getX() + vel * timeStep;
	setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ() ));
	setPoint(Vec(pos, getPoint().getY(), getPoint().getZ() )); }
      break;
    case 2:
      if (fabs(normal.getX()/areaX - sigma)/sigma > tol) {
	vel = ((normal.getX() - sigma*areaX)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getX() * (1-atf) / (1+atf) + (normal.getX() - sigma * areaX) / mass * timeStep / (1+atf);
	pos = prevPoint.getX() + vel * timeStep;
	setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ() ));
	setPoint(Vec(pos, getPoint().getY(), getPoint().getZ() )); }
      break;
    case 3:
      if (fabs(normal.getY()/areaY + sigma)/sigma > tol) {
	vel = ((normal.getY() + sigma*areaY)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getY() * (1-atf) / (1+atf) + (normal.getY() + sigma * areaY) / mass * timeStep / (1+atf);
	pos = prevPoint.getY() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ() ));
	setPoint(Vec(getPoint().getX(), pos, getPoint().getZ() )); }
      break;
    case 4:
      if (fabs(normal.getY()/areaY - sigma)/sigma > tol) {
	vel = ((normal.getY() - sigma*areaY)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getY() * (1-atf) / (1+atf) + (normal.getY() - sigma * areaY) / mass * timeStep / (1+atf);
	pos = prevPoint.getY() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ() ));
	setPoint(Vec(getPoint().getX(), pos, getPoint().getZ() )); }
      break;
    case 5:
      if (fabs(normal.getZ()/areaZ + sigma)/sigma > tol) {
	vel = ((normal.getZ() + sigma*areaZ)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getZ() * (1-atf) / (1+atf) + (normal.getZ() + sigma * areaZ) / mass * timeStep / (1+atf);
	pos = prevPoint.getZ() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel ));
	setPoint(Vec(getPoint().getX(), getPoint().getY(), pos )); }
      break;
    case 6:
      if (fabs(normal.getZ()/areaZ - sigma)/sigma > tol) {
	vel = ((normal.getZ() - sigma*areaZ)>0 ? 1:-1) * boundaryRate;
	if (isotropicType == 1 && normal.getZ() == 0 ) vel = -boundaryRate*topSpeedup;
	//vel = prevVeloc.getZ() * (1-atf) / (1+atf) + (normal.getZ() - sigma * areaZ) / mass * timeStep / (1+atf);
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
    REAL atf = forceDamp;

    REAL vel, pos;
    switch (id) {
    case 5:
      if (fabs(normal.getZ()/areaZ + sigma)/sigma > tol) {
	vel = ((normal.getZ() + sigma*areaZ)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getZ() * (1-atf) / (1+atf) + (normal.getZ() + sigma * areaZ) / mass * timeStep / (1+atf);
	pos = prevPoint.getZ() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel ));
	setPoint(Vec(getPoint().getX(), getPoint().getY(), pos )); }
      break;
    case 6:
      if (fabs(normal.getZ()/areaZ - sigma)/sigma > tol) {
	vel = ((normal.getZ() - sigma*areaZ)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getZ() * (1-atf) / (1+atf) + (normal.getZ() - sigma * areaZ) / mass * timeStep / (1+atf);
	pos = prevPoint.getZ() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel ));
	setPoint(Vec(getPoint().getX(), getPoint().getY(), pos )); }
      break;
    }
    prevPoint = point;
    prevVeloc = veloc;
  }

  // Normally it works well for boundary walls to displace outwards or inwards at at constant boundary rate,
  // depending on the difference of internal and external resultant forces. However, 
  // 1. for a smaller number of particles (which implies more inhomogeneous boundary contact), boundary walls 
  //    should follow the same equation of motion as particles, in order to sustain a constant confining pressure 
  //    in triaxial tests.
  // 2. for extremely-shaped particles (which implies an instable particle system), boundary walls should follow 
  //    the same equation of motion as particles, for the same reason.
  // How to use equation of motion for boundary walls?
  // 1. apply the same massScale as particles.
  // 2. let atf = 1 (critial background damping, usually leading to equilibrium) or atf = 0 (no background damping, often leading to zero tractions).
  // 3. virtualMass = 1.0E-7 for inhomogeneous case; = 1.0E-5 for extremely-shaped particles. These values need trial and error tests.
  void planeBoundary::updateTriaxial(REAL sigma, REAL areaX, REAL areaY, REAL areaZ) {
    std::size_t triaxialType = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["triaxialType"]);
    std::size_t unloadStep = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["unloadStep"]);

    REAL forceDamp = dem::Parameter::getSingleton().parameter["forceDamp"];
    REAL massScale = dem::Parameter::getSingleton().parameter["massScale"];
    //REAL mass = dem::Parameter::getSingleton().parameter["boundaryMass"];
    REAL boundaryRate = dem::Parameter::getSingleton().parameter["boundaryRate"];
    REAL tol = dem::Parameter::getSingleton().parameter["tractionErrorTol"];
    REAL atf = forceDamp;

//#define INHOMOGENEOUS 1
#ifdef INHOMOGENEOUS
    REAL virtualMass = 1.0E-7; // virtual mass of boundary walls for triaxial tests.
    REAL mass = virtualMass * massScale;
    atf = 1;
#endif

    REAL vel, pos;
    switch (id) {
    case 1: 
      if (fabs(normal.getX()/areaX + sigma)/sigma > tol) {
#ifdef INHOMOGENEOUS
	vel = prevVeloc.getX() * (1-atf) / (1+atf) + (normal.getX() + sigma * areaX) / mass * timeStep / (1+atf);
#else
	vel = ((normal.getX() + sigma*areaX)>0 ? 1:-1) * boundaryRate;
#endif
	pos = prevPoint.getX() + vel * timeStep;
	setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ() ));
	setPoint(Vec(pos, getPoint().getY(), getPoint().getZ() )); }
      break;
    case 2:
      if (fabs(normal.getX()/areaX - sigma)/sigma > tol) {
#ifdef INHOMOGENEOUS
	vel = prevVeloc.getX() * (1-atf) / (1+atf) + (normal.getX() - sigma * areaX) / mass * timeStep / (1+atf);
#else
	vel = ((normal.getX() - sigma*areaX)>0 ? 1:-1) * boundaryRate;
#endif
	pos = prevPoint.getX() + vel * timeStep;
	setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ() ));
	setPoint(Vec(pos, getPoint().getY(), getPoint().getZ() )); }
      break;
    case 3:
      if (fabs(normal.getY()/areaY + sigma)/sigma > tol) {
#ifdef INHOMOGENEOUS
	vel = prevVeloc.getY() * (1-atf) / (1+atf) + (normal.getY() + sigma * areaY) / mass * timeStep / (1+atf);
#else
	vel = ((normal.getY() + sigma*areaY)>0 ? 1:-1) * boundaryRate;
#endif
	pos = prevPoint.getY() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ() ));
	setPoint(Vec(getPoint().getX(), pos, getPoint().getZ() )); }
      break;
    case 4:
      if (fabs(normal.getY()/areaY - sigma)/sigma > tol) {
#ifdef INHOMOGENEOUS
	vel = prevVeloc.getY() * (1-atf) / (1+atf) + (normal.getY() - sigma * areaY) / mass * timeStep / (1+atf);
#else
	vel = ((normal.getY() - sigma*areaY)>0 ? 1:-1) * boundaryRate;
#endif
	pos = prevPoint.getY() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ() ));
	setPoint(Vec(getPoint().getX(), pos, getPoint().getZ() )); }
      break;

    // displacement control
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
      //vel = prevVeloc.getZ() * (1-atf) / (1+atf) + (normal.getZ() + sigma * areaZ) / mass * timeStep / (1+atf);
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
      //vel = prevVeloc.getZ() * (1-atf) / (1+atf) + (normal.getZ() - sigma * areaZ) / mass * timeStep / (1+atf);
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
    REAL atf = forceDamp;

    REAL vel, pos;
    switch (id) { // boundary x1(1) and boundary x2(2) do not move
    case 3:
      if (fabs(normal.getY()/areaY + sigma)/sigma > tol) {
	vel = ((normal.getY() + sigma*areaY)>0 ? 1:-1) * boundaryRate * sideRateRatio;
	//vel = prevVeloc.getY() * (1-atf) / (1+atf) + (normal.getY() + sigma * areaY) / mass * timeStep / (1+atf);
	pos = prevPoint.getY() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ() ));
	setPoint(Vec(getPoint().getX(), pos, getPoint().getZ() )); }
      break;
    case 4:
      if (fabs(normal.getY()/areaY - sigma)/sigma > tol) {
	vel = ((normal.getY() - sigma*areaY)>0 ? 1:-1) * boundaryRate * sideRateRatio;
	//vel = prevVeloc.getY() * (1-atf) / (1+atf) + (normal.getY() - sigma * areaY) / mass * timeStep / (1+atf);
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
      //vel = prevVeloc.getZ() * (1-atf) / (1+atf) + (normal.getZ() + sigma * areaZ) / mass * timeStep / (1+atf);
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
      //vel = prevVeloc.getZ() * (1-atf) / (1+atf) + (normal.getZ() - sigma * areaZ) / mass * timeStep / (1+atf);
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
    REAL atf = forceDamp;

    REAL vel, pos;
    switch (id) {
    case 1: 
      if (fabs(normal.getX()/areaX + sigmaX)/sigmaX > tol) {
	vel = ((normal.getX() + sigmaX*areaX)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getX() * (1-atf) / (1+atf) + (normal.getX() + sigma * areaX) / mass * timeStep / (1+atf);
	pos = prevPoint.getX() + vel * timeStep;
	setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ() ));
	setPoint(Vec(pos, getPoint().getY(), getPoint().getZ() )); }
      break;
    case 2:
      if (fabs(normal.getX()/areaX - sigmaX)/sigmaX > tol) {
	vel = ((normal.getX() - sigmaX*areaX)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getX() * (1-atf) / (1+atf) + (normal.getX() - sigma * areaX) / mass * timeStep / (1+atf);
	pos = prevPoint.getX() + vel * timeStep;
	setVeloc(Vec(vel, getVeloc().getY(), getVeloc().getZ() ));
	setPoint(Vec(pos, getPoint().getY(), getPoint().getZ() )); }
      break;
    case 3:
      if (fabs(normal.getY()/areaY + sigmaY)/sigmaY > tol) {
	vel = ((normal.getY() + sigmaY*areaY)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getY() * (1-atf) / (1+atf) + (normal.getY() + sigma * areaY) / mass * timeStep / (1+atf);
	pos = prevPoint.getY() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ() ));
	setPoint(Vec(getPoint().getX(), pos, getPoint().getZ() )); }
      break;
    case 4:
      if (fabs(normal.getY()/areaY - sigmaY)/sigmaY > tol) {
	vel = ((normal.getY() - sigmaY*areaY)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getY() * (1-atf) / (1+atf) + (normal.getY() - sigma * areaY) / mass * timeStep / (1+atf);
	pos = prevPoint.getY() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), vel, getVeloc().getZ() ));
	setPoint(Vec(getPoint().getX(), pos, getPoint().getZ() )); }
      break;
    case 5:
      if (fabs(normal.getZ()/areaZ + sigma)/sigma > tol) {
	vel = ((normal.getZ() + sigma*areaZ)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getZ() * (1-atf) / (1+atf) + (normal.getZ() + sigma * areaZ) / mass * timeStep / (1+atf);
	pos = prevPoint.getZ() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel ));
	setPoint(Vec(getPoint().getX(), getPoint().getY(), pos )); }
      break;
    case 6:
      if (fabs(normal.getZ()/areaZ - sigma)/sigma > tol) {
	vel = ((normal.getZ() - sigma*areaZ)>0 ? 1:-1) * boundaryRate;
	//vel = prevVeloc.getZ() * (1-atf) / (1+atf) + (normal.getZ() - sigma * areaZ) / mass * timeStep / (1+atf);
	pos = prevPoint.getZ() + vel * timeStep;
	setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel ));
	setPoint(Vec(getPoint().getX(), getPoint().getY(), pos )); }
      break;
    }
    prevPoint = point;
    prevVeloc = veloc;
  }

  void planeBoundary::updateOedometerImpact(REAL areaX, REAL areaY, REAL areaZ) {

    REAL forceDamp = dem::Parameter::getSingleton().parameter["forceDamp"];
    REAL massScale = dem::Parameter::getSingleton().parameter["massScale"];
    //REAL mass = dem::Parameter::getSingleton().parameter["boundaryMass"];
    REAL boundaryRate = dem::Parameter::getSingleton().parameter["boundaryRate"];
    REAL tol = dem::Parameter::getSingleton().parameter["tractionErrorTol"];
    REAL atf = forceDamp;

    REAL vel, pos;
    switch (id) {
    case 6:
      vel = -boundaryRate;
      //vel = prevVeloc.getZ() * (1-atf) / (1+atf) + (normal.getZ() - sigma * areaZ) / mass * timeStep / (1+atf);
      pos = prevPoint.getZ() + vel * timeStep;
      setVeloc(Vec(getVeloc().getX(), getVeloc().getY(), vel ));
      setPoint(Vec(getPoint().getX(), getPoint().getY(), pos ));
      break;
    }
    prevPoint = point;
    prevVeloc = veloc;
  }

  cylinderBoundary::cylinderBoundary(std::size_t tp, std::ifstream &ifs)
    :Boundary(tp, ifs) {
    REAL dx, dy, dz, px, py, pz;
    ifs >> dx >> dy >> dz >> px >> py >> pz >> radius;
    direc = Vec(dx, dy, dz);
    point = Vec(px, py, pz);
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
    // for each plane boundary, define a temporary variable vtmp to use,
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

