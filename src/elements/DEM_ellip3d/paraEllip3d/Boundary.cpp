#include "Boundary.h"
#include "Particle.h"
// use both pointer to and variable of class Particle

namespace dem {
  
  Boundary::Boundary(std::ifstream &ifs) {
    BdryCoef tmp;
    REAL x,y,z;
    coefOfLimits.clear();
    ifs >> boundaryId >> limitNum;
    for (int k = 0; k < limitNum; k++){
      ifs >> tmp.order >> x >> y >> z;
      tmp.dirc = Vec(x,y,z);
      ifs >> x >> y >> z;
      tmp.apt = Vec(x,y,z);
      ifs >> tmp.rad >> tmp.side;
      coefOfLimits.push_back(tmp);
    }
  }
  
  void Boundary::update(UPDATECTL &ctl) {
    BdryCoef tmp;
    std::vector<BdryCoef>::iterator it;
    Vec nv, napt;
    for (it = coefOfLimits.begin(); it != coefOfLimits.end(); ++it){
      tmp = *it;
      nv = rotateVec(tmp.dirc,ctl.rote);
      napt = ctl.tran+ctl.fixpt+rotateVec(tmp.apt-ctl.fixpt, ctl.rote);
      (*it).dirc = nv;
      (*it).apt = napt;
    }
  }  
  
  Vec plnBoundary::getApt() const {
    return (*this->coefOfLimits.begin()).apt;
  }
  
  Vec plnBoundary::getDirc() const {
    return (*this->coefOfLimits.begin()).dirc;
  }
  
  void plnBoundary::display() const {
    Boundary::display();
    std::cout << "normal: " << normal.getX() << " " << normal.getY() << " " << normal.getZ() << std::endl;
    typename std::vector<Particle *>::const_iterator it;
    int i = 0;
    for(it = possBdryParticle.begin();it != possBdryParticle.end(); ++it){
      if(i++ < 10)
	std::cout << (*it)->getId();
      else {
	i = 0;
	std::cout << (*it)->getId() << std::endl;
      }
    }
  }
  
  REAL plnBoundary::distToBdry(Vec posi) const {
    Vec dv=(*this->coefOfLimits.begin()).dirc;
    Vec pt=(*this->coefOfLimits.begin()).apt;
    Vec ndv=normalize(dv);
    return (posi-pt) % ndv;
  }
  
  void plnBoundary::findBdryContact(std::vector<Particle *>& ptcls) {
    typename std::vector<Particle *>::iterator it;
    std::vector<BdryCoef>::iterator bt;
    bool next;
    possBdryParticle.clear();
    REAL dist, r;
    Vec posi, ndirc;
    for (it = ptcls.begin(); it != ptcls.end(); ++it){
      if ( (*it)->getType() == 0 ) { // only process free particles, excluding type 5
	posi=(*it)->getCurrPos();
	dist=distToBdry(posi);
	if(dist >= 0 || fabs(dist) > (*it)->getA()) // outside to coefOfLimits[0] or inside too much
	  continue;
	/*
	  debugInf << "boundary.h: iter=" << iteration
	  << " bdryId=" << getBdryID()
	  << " ptclId=" << (*it)->getId() << std::endl;
	*/
	next = true;
	for (bt = ++this->coefOfLimits.begin(); bt != this->coefOfLimits.end(); ++bt) { // coefOfLimits[1,2,...]
	  ndirc = normalize((*bt).dirc);
	  r = vfabs((posi-(*bt).apt) - (posi-(*bt).apt) % ndirc * ndirc);
	  if( ( (*bt).order == 1 && (posi-(*bt).apt)%(*bt).dirc >= 0 ) ||
	      ( (*bt).order == 2 && (r-(*bt).rad)*(*bt).side < 0 ) ){
	    next = false; // the particle is out of boundary, process next particle
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
  }
  
  void plnBoundary::boundaryForce(std::map<int,std::vector<BoundaryTgt> > &boundaryTgtMap) {
    typename std::vector<Particle *>::iterator it;
    this->avgNormal = 0;
    this->avgPenetr = 0;
    this->contactNum = 0;
    normal = 0;
    tangt = 0;
    moment = 0;
    
    // for each plane boundary, define a temparory variable vtmp to use,
    // better than define a member variable which needs to be cleared.
    // and vtmp is initialized as empty in each iteration.
    std::vector<BoundaryTgt> vtmp;
    
    // for each possible boundary particle
    REAL penetr = 0;
    int count = 0;
    for (it = possBdryParticle.begin(); it != possBdryParticle.end(); ++it){
      penetr = 0;
      (*it)->planeRBForce(this, boundaryTgtMap, vtmp, penetr);
      this->avgPenetr += penetr;
      count++;
    }
    if (count>0) this->avgPenetr /= count;
    this->contactNum = count;
    
    // checkout tangential forces and displacements after each particle is processed
    boundaryTgtMap[this->boundaryId] = vtmp;
  }
  
  void cylBoundary::display() const {
    Boundary::display();
    std::cout << "normal: " << normal.getX() << " " << normal.getY() << " " << normal.getZ() << std::endl;
    typename std::vector<Particle *>::const_iterator it;
    int i = 0;
    for(it = possBdryParticle.begin();it!=possBdryParticle.end();++it){
      if(i++<10)
	std::cout << (*it)->getId();
      else{
	i = 0;
	std::cout << (*it)->getId() << std::endl;
      }
    }
  }
  
  REAL cylBoundary::distToBdry(Vec posi) const {
    Vec ndc = (*this->coefOfLimits.begin()).dirc;
    Vec napt = (*this->coefOfLimits.begin()).apt;
    REAL r = (*this->coefOfLimits.begin()).rad;
    Vec norm = normalize(ndc);
    return fabs(r-vfabs((posi-napt)-(posi-napt)%norm*norm));
  }
  
  void cylBoundary::findBdryContact(std::vector<Particle *> &ptcls) {
    typename std::vector<Particle *>::iterator it;
    std::vector<BdryCoef>::iterator bt;
    bool next;
    possBdryParticle.clear();
    REAL dist,r;
    Vec posi, ndirc;
    for (it = ptcls.begin();it!=ptcls.end();++it){
      posi = (*it)->getCurrPos();
      dist = distToBdry(posi);
      if (dist>(*it)->getA())
	continue;
      next = false;
      for (bt = ++this->coefOfLimits.begin();bt!=this->coefOfLimits.end();++bt){
	ndirc = normalize((*bt).dirc);
	r = vfabs((posi-(*bt).apt)-(posi-(*bt).apt)%ndirc*ndirc);
	if( ( (*bt).order==1&&(posi-(*bt).apt)%(*bt).dirc>(*it)->getA() )||
	    ( (*bt).order==2&&(r-(*bt).rad)*(*bt).side<0 ) ){
	  next = true;//the particle is outof boundary, process next particle
	  break;
	}
      }
      if(!next)
	possBdryParticle.push_back(*it);
    }
  }
  
  void cylBoundary::boundaryForce() {
    Cylinder cyl;
    typename std::vector<Particle *>::iterator it;
    BdryCoef tmp;
    tmp = *this->coefOfLimits.begin();
    cyl.setRadius(tmp.rad);
    cyl.setCenter(tmp.apt);
    normal = 0;
    for (it = possBdryParticle.begin();it!=possBdryParticle.end();++it){
      normal-=(*it)->cylinderRBForce(this->boundaryId,cyl,tmp.side);
    }
  }
  
  
} // namespace dem ends

