#include "Boundary.h"
#include "Particle.h"
// use both pointer to and variable of class Particle

namespace dem {
  
  Boundary::Boundary(int tp, std::ifstream &ifs) {
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
  
  planeBoundary::planeBoundary(int tp, std::ifstream &ifs)
    :Boundary(tp, ifs) {
    REAL dx, dy, dz, px, py, pz;
    ifs >> dx >> dy >> dz >> px >> py >> pz;
    direc = Vec(dx, dy, dz);
    point = Vec(px, py, pz);
    for (int i = 0; i < extraNum; ++i) {
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
  
  void planeBoundary::boundaryForce(std::map<int,std::vector<BoundaryTgt> > &boundaryTgtMap) {
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

  cylinderBoundary::cylinderBoundary(int tp, std::ifstream &ifs)
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
  
  void cylinderBoundary::boundaryForce(std::map<int,std::vector<BoundaryTgt> > &boundaryTgtMap) {
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

