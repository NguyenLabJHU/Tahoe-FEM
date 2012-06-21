//  This is a template class, for which we have to include the implementation in the header file.
//  As we cannot put using statement in a header file, we have to use std::something wherever we
//  need to refer anything from standard namespace.

#ifndef CONTACT_H
#define CONTACT_H

#include "realtypes.h"
#include "parameter.h"
#include "root6.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

//#define MINDLIN_ASSUMED
//#define MINDLIN_KNOWN

namespace dem {
  
  class ContactTgt {
    
  public:
    int  ptcl1;
    int  ptcl2;
    Vec  tgtForce;
    Vec  tgtDisp;
    bool tgtLoading;
    Vec  tgtDispStart;
    REAL tgtPeak;
    bool tgtSlide;
    
    ContactTgt()
      :ptcl1(0),
      ptcl2(0),
      tgtForce(0),
      tgtDisp(0),
      tgtLoading(0),
      tgtDispStart(0),
      tgtPeak(0),
      tgtSlide(false)
      {}
    
    ContactTgt(int _ptcl1, int _ptcl2, Vec _tf, Vec _td, bool _tl, Vec _tds, REAL _tp, bool _ts)
      :ptcl1(_ptcl1),
      ptcl2(_ptcl2),
      tgtForce(_tf),
      tgtDisp(_td), 
      tgtLoading(_tl),
      tgtDispStart(_tds),
      tgtPeak(_tp),
      tgtSlide(_ts)
      {}
  };

  template <class T> class Contact {
  public:  
    Contact();
    Contact(T* t1, T* t2);
    
    T*   getP1() const;
    T*   getP2() const;
    Vec  getPoint1() const {return point1;}
    Vec  getPoint2() const {return point2;}
    REAL getRadius1() const {return radius1;}
    REAL getRadius2() const {return radius2;}
    REAL getR0() const {return R0;}
    REAL getE0() const {return E0;}
    REAL getVibraTimeStep() const {return vibraTimeStep;}
    REAL getImpactTimeStep() const {return impactTimeStep;}
    
    bool isOverlapped();
    void contactForce(); // calculate normal and tangential force of contact
    REAL getNormalForce() const {return vfabs(normalForce);}
    REAL getTgtForce()  const {return vfabs(tgtForce);}
    REAL getPenetration() const {return penetr;}
    REAL getContactRadius() const {return contactRadius;}
    REAL gettgtDisp() const {return vfabs(tgtDisp);} // total value during a process of contact
    void checkoutTgt(std::vector<ContactTgt>& ContactTgtVec);
    void checkinPreTgt(std::vector<ContactTgt>& ContactTgtVec);
    Vec  normalForceVec() const {return normalForce;}
    Vec  tgtForceVec() const {return tgtForce;}
    bool isRedundant(Contact<T> other) const;
    
  private:
    T*   p1;            // particle 1
    T*   p2;            // particle 2
    REAL penetr;        // penetr
    REAL contactRadius; // radius of contact surface
    Vec  point1;        // point1 on particle 1, innermost to particle 2
    Vec  point2;        // point2 on particle 2, innermost to particle 1
    REAL radius1;       // radius of osculating circles at point1
    REAL radius2;       // radius of osculating circles at point2
    REAL E0;              
    REAL G0;
    REAL R0;
    REAL vibraTimeStep;
    REAL impactTimeStep;
    bool isInContact;
    bool tgtLoading;    // tangential loading or unloading
    Vec  normalForce;   // positive when pointing to paticle 1
    Vec  tgtForce;      // TgtrDirc points along tangential forces exerted on particle 1
    Vec  tgtDisp;       // tangential relative displacment total vector
    Vec  tgtDispStart;  // displacement start value for each loading-unloading loop
    bool tgtSlide;
    Vec  normalDirc;    
    Vec  tgtDirc;
    Vec  cohesionForce;  // cohesion force between particles
    bool prevTgtLoading; // previous loading-unloading status
    Vec  prevNormalForce;
    Vec  prevTgtForce;
    Vec  prevTgtDisp;    // previous tangential relative displacment total vector
    bool prevTgtSlide;
    REAL tgtPeak;       
    Vec  spin_res;
};


  template <class T>
    Contact<T>::Contact(){
    p1 = NULL;
    p2 = NULL;
    isInContact = false;
    tgtLoading = prevTgtLoading = true;
    tgtPeak = 0;
    penetr = 0;
    contactRadius = 0;
    radius1 = radius2 = 0;
    normalForce = prevNormalForce = 0;
    tgtForce = prevTgtForce = 0;
    tgtDisp = prevTgtDisp = 0;
    tgtDispStart = 0;
    normalDirc = 0;
    tgtDirc = 0;
    spin_res = 0;
  }
  
  
  template <class T>
    Contact<T>::Contact(T* t1, T* t2){
    p1 = t1;
    p2 = t2;
    isInContact = false;
    tgtLoading = prevTgtLoading = true;
    tgtPeak = 0;
    penetr = 0;
    contactRadius = 0;
    radius1 = radius2 = 0;
    normalForce = prevNormalForce = 0;
    tgtForce = prevTgtForce = 0;
    tgtDisp = prevTgtDisp = 0;
    tgtDispStart = 0;
    normalDirc = 0;
    tgtDirc = 0;
    spin_res = 0;
  }
  
  template<class T>
    bool Contact<T>::isRedundant(Contact<T> other) const {
    int id1 = getP1() -> getId();
    int id2 = getP2() -> getId();
    int oid1 = ( other.getP1() ) -> getId();
    int oid2 = ( other.getP2() ) -> getId();
    
    if ( (id2 == oid1 && id1 == oid2) || (id1 == oid1 && id2 == oid2 ) ) {
      //std::cout << id1 << " " << id2 << " " << oid1 << " " << oid2 << " " << std::endl; 
      return true;}
    else 
      return false;
  }
  
  template<class T>
    T* Contact<T>::getP1() const {
    return p1;
  }
  

  template<class T>
    T* Contact<T>::getP2() const {
    return p2;
  }
  

  template<class T>
    bool Contact<T>::isOverlapped(){
    REAL coef1[10],coef2[10];
    p1->getGlobalCoef(coef1); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobalCoef(coef2);    
    Vec v[2];
    bool b1 = root6(coef1,coef2,v[0]);
    bool b2 = root6(coef2,coef1,v[1]);
    point1 = v[1];
    point2 = v[0];
    radius1 = p1->getRadius(point1);
    radius2 = p2->getRadius(point2);
    penetr = vfabs(point1-point2);
    REAL minRelOverlap = penetr/(2.0*fmax(radius1,radius2));

    if (b1 && b2 && minRelOverlap > MINOVERLAP) { // a strict detection method
        isInContact = true;
        return true;
    }
    else {
        isInContact = false;
	return false;
    }
  }
  
  
  template<class T>
    void Contact<T>::checkinPreTgt(std::vector<ContactTgt>& ContactTgtVec) {
    if (ContactTgtVec.size()>0) {
      for(std::vector<ContactTgt>::iterator it = ContactTgtVec.begin();it != ContactTgtVec.end(); ++it) {
	if (it->ptcl1 == p1->getId() && it->ptcl2 == p2->getId()) {
	  prevTgtForce   = it->tgtForce;
	  prevTgtDisp    = it->tgtDisp;
	  prevTgtLoading = it->tgtLoading;
	  prevTgtSlide   = it->tgtSlide;
	  tgtDispStart   = it->tgtDispStart;
	  tgtPeak        = it->tgtPeak;
	  break;
	}
      }
    }
  }
  
  
  template<class T>
    void Contact<T>::checkoutTgt(std::vector<ContactTgt>& ContactTgtVec) {
    ContactTgtVec.push_back(ContactTgt(p1->getId(),p2->getId(),
				       tgtForce,
				       tgtDisp, 
				       tgtLoading,
				       tgtDispStart,
				       tgtPeak,
				       tgtSlide) );
  }  
  
  
  template<class T>
    void Contact<T>::contactForce(){
    // isOverlapped() has been called in findContact() in assembly.cpp and information recorded, 
    // now this function is called by internalForce() in assembly.cpp.
    
    if (isInContact) {
      // obtain normal force, using absolute equation instead of stiffness method
      p1->setContactNum(p1->getContactNum() + 1);
      p2->setContactNum(p2->getContactNum() + 1);
      p1->setInContact(true);
      p2->setInContact(true);
      
      R0 = radius1*radius2/(radius1+radius2);
      E0 = 0.5*YOUNG/(1-POISSON*POISSON);
      REAL allowedOverlap = 2.0 * fmin(radius1,radius2) * MAXOVERLAP;
      if (penetr > allowedOverlap) {
	debugInf << "Contact.h: iter=" << iteration 
		 << " ptcl1=" << getP1()->getId()
		 << " ptcl2=" << getP2()->getId()
		 << " penetr=" << penetr 
		 << " allow=" << allowedOverlap << std::endl;
	penetr = allowedOverlap;
      }
#ifdef MEASURE_EPS
      penetr = nearbyint (penetr/MEPS) * MEPS;
#endif
      contactRadius = sqrt(penetr*R0);
      normalDirc = normalize(point1-point2);         // normalDirc points out of particle 1
      normalForce = -sqrt(penetr*penetr*penetr)*sqrt(R0)*4*E0/3* normalDirc; // normalForce pointing to particle 1
      // pow(penetr, 1.5)
      
      // apply cohesion force
      cohesionForce = PI*(penetr*R0)*COHESION*normalDirc;
      p1->addForce(cohesionForce);
      p2->addForce(-cohesionForce);
      
      // apply normal force
      p1->addForce(normalForce);
      p2->addForce(-normalForce);
      p1->addMoment( ( (point1+point2)/2-p1->getCurrPos() ) *   normalForce );
      p2->addMoment( ( (point1+point2)/2-p2->getCurrPos() ) * (-normalForce) );	
      
      /*
      debugInf << "Contact.h: iter=" << iteration
	       << " penetr=" << penetr
	       << " cohesionForce=" << vfabs(cohesionForce)
	       << " normalForce=" << vfabs(normalForce)
	       << " accumulated time=" << iteration*TIMESTEP
	       << std::endl;
      */
      
      // obtain normal damping force
      Vec cp = (point1+point2)/2;        
      Vec veloc1 = p1->getCurrVeloc() + p1->getCurrOmga()*(cp-p1->getCurrPos());
      Vec veloc2 = p2->getCurrVeloc() + p2->getCurrOmga()*(cp-p2->getCurrPos());
      REAL m1 = getP1()->getMass();
      REAL m2 = getP2()->getMass();
      REAL kn = pow(6*vfabs(normalForce)*R0*pow(E0,2),1.0/3.0);
      REAL DMP_CRTC = 2*sqrt(m1*m2/(m1+m2)*kn); // critical damping
      Vec CntDampingForce = DMP_CNT * DMP_CRTC * ((veloc1-veloc2)%normalDirc)*normalDirc;
      vibraTimeStep = 2.0*sqrt( m1*m2 / (m1+m2) /kn );
      impactTimeStep = allowedOverlap / fabs((veloc1-veloc2) % normalDirc);
      
      // apply normal damping force
      p1->addForce(-CntDampingForce);
      p2->addForce(CntDampingForce);
      p1->addMoment( ( (point1+point2)/2-p1->getCurrPos() ) * (-CntDampingForce) );
      p2->addMoment( ( (point1+point2)/2-p2->getCurrPos() ) * CntDampingForce );
      
      if (FRICTION != 0) {
	// obtain tangential force
	G0 = YOUNG/2/(1+POISSON);              // RelaDispInc points along point1's displacement relative to point2
	Vec RelaDispInc = (veloc1-veloc2)*TIMESTEP;
	Vec tgtDispInc = RelaDispInc-(RelaDispInc%normalDirc)*normalDirc;
	tgtDisp = prevTgtDisp + tgtDispInc; // prevTgtDisp read by checkinPreTgt()
	if (vfabs(tgtDisp) == 0)
	  tgtDirc = 0;
	else
	  tgtDirc = normalize(-tgtDisp); // tgtDirc points along Tgtential forces exerted on particle 1
	
	REAL fP = 0;
	REAL ks = 0;
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// linear friction model
	fP = FRICTION*vfabs(normalForce);
	ks = 4*G0*contactRadius/(2-POISSON);
	tgtForce = prevTgtForce + ks*(-tgtDispInc); // prevTgtForce read by CheckinPreTgt()
	if (vfabs(tgtForce) > fP)
	  tgtForce = fP*tgtDirc;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Mindlin's model (loading/unloading condition assumed)
	// This model is not recommended as it is impossible to strictly determine loading/unloading condition
	// unless load is known (the case of pure moment rotation).
#ifdef MINDLIN_ASSUMED
	REAL val = 0;
	fP = FRICTION*vfabs(normalForce);
	tgtLoading = (prevTgtDisp%tgtDispInc >= 0); 
	
	if (tgtLoading) {              // loading
	  if (!prevTgtLoading) {      // pre-step is unloading
	    val = 8*G0*contactRadius*vfabs(tgtDispInc)/(3*(2-POISSON)*fP);
	    tgtDispStart = prevTgtDisp;
	  }
	  else                       // pre-step is loading
	    val = 8*G0*contactRadius*vfabs(tgtDisp-tgtDispStart)/(3*(2-POISSON)*fP);
	  
	  if (val > 1.0)              
	    tgtForce = fP*tgtDirc;
	  else {
	    ks = 4*G0*contactRadius/(2-POISSON)*sqrt(1-val);
	    //incremental method
	    tgtForce = prevTgtForce + ks*(-tgtDispInc); // tgtDispInc determines signs
	    //total value method: tgtForce = fP*(1-pow(1-val, 1.5))*tgtDirc;
	  }
	}
	else {                         // unloading
	  if (prevTgtLoading) {       // pre-step is loading
	    val = 8*G0*contactRadius*vfabs(tgtDisp-tgtDispStart)/(3*(2-POISSON)*fP);
	    tgtPeak = vfabs(prevTgtForce);
	  }
	  else                       // pre-step is unloading
	    val = 8*G0*contactRadius*vfabs(tgtDisp-tgtDispStart)/(3*(2-POISSON)*fP);
	  
	  if (val > 1.0 || tgtPeak > fP)  
	    tgtForce = fP*tgtDirc;
	  else {
	    ks = 2*sqrt(2)*G0*contactRadius/(2-POISSON) * sqrt(1+pow(1-tgtPeak/fP,2.0/3.0)+val);
	    //incremental method
	    tgtForce = prevTgtForce + ks*(-tgtDispInc); // tgtDispInc determines signs
	    //total value method: tgtForce = (tgtPeak-2*fP*(1-sqrt(2)/4*pow(1+ pow(1-tgtPeak/fP,2.0/3.0) + val,1.5)))*tgtDirc;
	  }
	}
	
	if (vfabs(tgtForce) > fP)
	  tgtForce = fP*tgtDirc;
#endif
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Mindlin's model (loading/unloading condition known for pure moment rotation case)
	// As loading/unloading condition is known, both incremental and total value method work well.
	// Herein sliding history is incorporated.
#ifdef MINDLIN_KNOWN
	REAL val = 0;
	fP = FRICTION*vfabs(normalForce);
	if (prevTgtSlide)
	  val = 8*G0*contactRadius*vfabs(tgtDispInc)/(3*(2-POISSON)*fP);
	else
	  val = 8*G0*contactRadius*vfabs(tgtDisp-tgtDispStart)/(3*(2-POISSON)*fP);
	
	if (iteration > 10000 && iteration < 11000) { // loading (and possible sliding)
	  if (val > 1.0) {
	    tgtForce = fP*tgtDirc;
	    tgtSlide = true;
	  }
	  else {
	    if (!prevTgtSlide) {
	      ks = 4*G0*contactRadius/(2-POISSON)*sqrt(1-val);
	      tgtForce = prevTgtForce + ks*(-tgtDispInc); // tgtDispInc determines signs
	      tgtSlide = false;
	    }
	    else {
	      if (vfabs(tgtForce)>vfabs(prevTgtForce))
		tgtSlide = true;
	      else
		tgtSlide = false;
	    }
	  }
	  tgtPeak = vfabs(tgtForce);
	}
	else { // (possible sliding and) unloading
	  if (val > 1.0 || tgtPeak > fP) {  
	    tgtForce = fP*tgtDirc;
	    tgtSlide = true;
	  }
	  else {
	    if (!prevTgtSlide) {
	      ks = 2*sqrt(2)*G0*contactRadius/(2-POISSON) * sqrt(1+pow(1-tgtPeak/fP,2.0/3.0)+val);
	      tgtForce = prevTgtForce + ks*(-tgtDispInc); // tgtDispInc determines signs
	      tgtSlide = false;
	    }
	    else {
	      if (vfabs(tgtForce)>vfabs(prevTgtForce))
		tgtSlide = true;
	      else {
		tgtSlide = false;
		tgtDispStart = tgtDisp;
	      }
	    }
	  }
	}
	
	/*
	debugInf << "Contact.h: iter="<iteration
		 << " prevTgtSlide=" << prevTgtSlide
		 << " tgtSlide=" << tgtSlide
		 << " val=" << val
		 << " ks=" << ks
		 << " tgtDispInc.x=" << tgtDispInc.getx()
		 << " prevTgtForce=" << vfabs(prevTgtForce)
		 << " tgtForce" << vfabs(tgtForce)
		 << std::endl;
	*/
	
	if (vfabs(tgtForce) > fP)
	  tgtForce = fP*tgtDirc;
#endif	    
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// apply tangential force
	p1->addForce(tgtForce);
	p2->addForce(-tgtForce);
	p1->addMoment( ((point1+point2)/2-p1->getCurrPos())*tgtForce);
	p2->addMoment(-((point1+point2)/2-p2->getCurrPos())*tgtForce);
      }
      
    }
    else {
      isInContact = false;
      tgtLoading = false;
      tgtPeak = 0;
      normalForce = 0;
      tgtForce = 0;
      tgtDisp = 0;    //total value
      normalDirc = 0;
      tgtDirc = 0;
      
      penetr = 0;
      contactRadius = 0;
      radius1 = radius2 = 0;
      spin_res = 0;
    }   
    
  }
  
} // namespace dem ends

#endif
