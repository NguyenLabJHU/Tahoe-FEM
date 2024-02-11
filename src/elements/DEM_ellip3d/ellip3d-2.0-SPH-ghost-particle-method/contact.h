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

class cnttgt{
public:
    int  ptcl1;
    int  ptcl2;
    vec  TgtForce;
    vec  TgtDisp;
    bool TgtLoading;
    vec  TgtDispStart;
    REAL TgtPeak;
    bool TgtSlide;

    cnttgt();
    cnttgt(int _ptcl1, int _ptcl2, vec _tf, vec _td, bool _tl, vec _tds, REAL _tp, bool _ts)
      :ptcl1(_ptcl1), ptcl2(_ptcl2), TgtForce(_tf), TgtDisp(_td), 
       TgtLoading(_tl),
       TgtDispStart(_tds),
       TgtPeak(_tp),
       TgtSlide(_ts)
       {};
};


template <class T> class contact{
public:  
    contact();
    contact(T* t1, T* t2);
    
    T*   getP1() const;
    T*   getP2() const;
    vec  getPoint1() const {return point1;}
    vec  getPoint2() const {return point2;}
    REAL getRadius1() const {return radius1;}
    REAL getRadius2() const {return radius2;}
    REAL getR0() const {return R0;}
    REAL getE0() const {return E0;}
    REAL getVibraTimeStep() const {return vibraTimeStep;}
    REAL getImpactTimeStep() const {return impactTimeStep;}
    
    bool isOverlapped();
    void contactForce(bool &exceed);         // calculate normal and tangential force of contact
    REAL getNormalForce() const {return vfabs(NormalForce);}
    REAL getTgtForce()  const {return vfabs(TgtForce);}
    REAL getPenetration() const {return penetration;}
    REAL getContactRadius() const {return contact_radius;}
    REAL getTgtDisp() const {return vfabs(TgtDisp);} // total value during a process of contact
    void checkoutTgt(std::vector<cnttgt>& CntTgtVec);
    void checkinPreTgt(std::vector<cnttgt>& CntTgtVec);
    void print() const;
    vec NormalForceVec() const {return NormalForce;}
    vec TgtForceVec() const {return TgtForce;}
    bool isRedundant(contact<T> other) const;
    
 private:
    T*   p1;              // particle 1
    T*   p2;              // particle 2
    REAL penetration;     // penetration
    REAL contact_radius;  // radius of contact surface
    vec  point1;          // point1 on particle 1, innermost to particle 2
    vec  point2;          // point2 on particle 2, innermost to particle 1
    REAL radius1;         // radius of osculating circles at point1
    REAL radius2;         // radius of osculating circles at point2
    REAL E0;              
    REAL G0;
    REAL R0;
    REAL vibraTimeStep;
    REAL impactTimeStep;

    bool isInContact;

    bool TgtLoading;           // tangential loading or unloading
    vec  NormalForce;          // positive when pointing to paticle 1
    vec  TgtForce;             // TgtrDirc points along tangential forces exerted on particle 1
    vec  TgtDisp;              // tangential relative displacment total vector
    vec  TgtDispStart;         // displacement start value for each loading-unloading loop
    bool TgtSlide;
    vec  NormDirc;    
    vec  TgtDirc;
    vec  CohesionForce;        // cohesion force between particles

    bool PreTgtLoading;        // previous loading-unloading status
    vec  PreNormalForce;
    vec  PreTgtForce;
    vec  PreTgtDisp;           // previous tangential relative displacment total vector
    bool PreTgtSlide;

    REAL TgtPeak;       

    vec  spin_res;
};


template <class T>
contact<T>::contact(){
    p1=NULL;
    p2=NULL;
    isInContact=false;
    TgtLoading=PreTgtLoading=true;
    TgtPeak=0;
    penetration=0;
    contact_radius=0;
    radius1=radius2=0;
    NormalForce=PreNormalForce=0;
    TgtForce=PreTgtForce=0;
    TgtDisp=PreTgtDisp=0;
    TgtDispStart=0;
    NormDirc=0;
    TgtDirc=0;
    spin_res=0;
}


template <class T>
contact<T>::contact(T* t1, T* t2){
    p1=t1;
    p2=t2;
    isInContact=false;
    TgtLoading=PreTgtLoading=true;
    TgtPeak=0;
    penetration=0;
    contact_radius=0;
    radius1=radius2=0;
    NormalForce=PreNormalForce=0;
    TgtForce=PreTgtForce=0;
    TgtDisp=PreTgtDisp=0;
    TgtDispStart=0;
    NormDirc=0;
    TgtDirc=0;
    spin_res=0;
}

template<class T>
bool contact<T>::isRedundant(contact<T> other) const {
  int id1 = getP1() -> getID();
  int id2 = getP2() -> getID();
  int oid1 = ( other.getP1() ) -> getID();
  int oid2 = ( other.getP2() ) -> getID();

  if ( (id2 == oid1 && id1 == oid2) || (id1 == oid1 && id2 == oid2 ) ) {
    //std::cout << id1 << " " << id2 << " " << oid1 << " " << oid2 << " " << std::endl; 
    return true;}
  else 
    return false;
}

template<class T>
T* contact<T>::getP1() const {
    return p1;
}


template<class T>
T* contact<T>::getP2() const {
    return p2;
}


template<class T>
bool contact<T>::isOverlapped(){
    REAL coef1[10],coef2[10];
    p1->getGlobCoef(coef1); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef2);    
    vec v[2];
    bool b1 = root6(coef1,coef2,v[0]);
    bool b2 = root6(coef2,coef1,v[1]);
    point1 = v[1];
    point2 = v[0];
    radius1=p1->getRadius(point1);
    radius2=p2->getRadius(point2);
    penetration=vfabs(point1-point2);
    REAL minRelOverlap = penetration/(2.0*fmax(radius1,radius2));

    if (b1 && b2 
	&& minRelOverlap > MINOVERLAP
	&& nearbyint(penetration/MEPS >= 1) ) {
        isInContact = true;
        return true;
    }
    else {
        isInContact = false;
	return false;
    }
}


template<class T>
void contact<T>::checkinPreTgt(std::vector<cnttgt>& CntTgtVec) {
    if (CntTgtVec.size()>0) {
	for(std::vector<cnttgt>::iterator it=CntTgtVec.begin();it!=CntTgtVec.end();++it) {
	    if (it->ptcl1==p1->getID() && it->ptcl2==p2->getID()) {
		PreTgtForce   = it->TgtForce;
		PreTgtDisp    = it->TgtDisp;
		PreTgtLoading = it->TgtLoading;
		PreTgtSlide   = it->TgtSlide;
		TgtDispStart  = it->TgtDispStart;
		TgtPeak       = it->TgtPeak;
		break;
	    }
	}
    }
}


template<class T>
void contact<T>::checkoutTgt(std::vector<cnttgt>& CntTgtVec) {
    CntTgtVec.push_back(cnttgt(p1->getID(),p2->getID(),
			       TgtForce,
			       TgtDisp, 
			       TgtLoading,
			       TgtDispStart,
			       TgtPeak,
			       TgtSlide) );
}  


template<class T>
void contact<T>::contactForce(bool &exceed){
    // isOverlapped() has been called in findContact() in assembly.cpp and information recorded, 
    // now this function is called by internalForce() in assembly.cpp.

    if (isInContact) {
	// obtain normal force, using absolute equation instead of stiffness method
	p1->cntnum++;
	p2->cntnum++;
	p1->inContact = true;
	p2->inContact = true;

	R0=radius1*radius2/(radius1+radius2);
	E0=0.5*YOUNG/(1-POISSON*POISSON);
	REAL allowedOverlap = 2.0 * fmin(radius1,radius2) * MAXOVERLAP;
	if (penetration > allowedOverlap) {
	  g_debuginf << "contact.h: g_iter=" << g_iteration 
		     << " ptcl1=" << getP1()->getID()
		     << " ptcl2=" << getP2()->getID()
		     << " penetr=" << penetration 
		     << " allow=" << allowedOverlap << std::endl;
	  //if (penetration > 1.0e-3) exceed = true;
	  penetration = allowedOverlap;
	}
#ifdef MEASURE_EPS
	penetration = nearbyint (penetration/MEPS) * MEPS;
#endif
	contact_radius=sqrt(penetration*R0);
	NormDirc=normalize(point1-point2);         // NormDirc points out of particle 1
	NormalForce= -sqrt(penetration*penetration*penetration)*sqrt(R0)*4*E0/3* NormDirc; // NormalForce pointing to particle 1
	// pow(penetration, 1.5)

        // apply cohesion force
	CohesionForce=PI*(penetration*R0)*COHESION*NormDirc;
	p1->addForce(CohesionForce);
	p2->addForce(-CohesionForce);

	// apply normal force
	p1->addForce(NormalForce);
	p2->addForce(-NormalForce);
	p1->addMoment( ( (point1+point2)/2-p1->getCurrPosition() ) *   NormalForce );
	p2->addMoment( ( (point1+point2)/2-p2->getCurrPosition() ) * (-NormalForce) );	
	
	/*
	g_debuginf<<"contact.h: g_iter="<<g_iteration
		  <<" penetr="<<penetration
		  <<" CohesionForce="<<vfabs(CohesionForce)
		  <<" NormalForce="<<vfabs(NormalForce)
		  <<" accumulated time="<<g_iteration*TIMESTEP
		  <<std::endl;
	*/

	// obtain normal damping force
	vec cp = (point1+point2)/2;        
	vec veloc1 = p1->getCurrVelocity() + p1->getCurrOmga()*(cp-p1->getCurrPosition());
	vec veloc2 = p2->getCurrVelocity() + p2->getCurrOmga()*(cp-p2->getCurrPosition());
	REAL m1 = getP1()->getMass();
	REAL m2 = getP2()->getMass();
	REAL kn = pow(6*vfabs(NormalForce)*R0*pow(E0,2),1.0/3.0);
	REAL DMP_CRTC = 2*sqrt(m1*m2/(m1+m2)*kn); // critical damping
	vec CntDampingForce  = DMP_CNT * DMP_CRTC * ((veloc1-veloc2)%NormDirc)*NormDirc;
	vibraTimeStep = 2.0*sqrt( m1*m2 / (m1+m2) /kn );
	impactTimeStep = allowedOverlap / fabs((veloc1-veloc2) % NormDirc);

	// apply normal damping force
	p1->addForce(-CntDampingForce);
	p2->addForce(CntDampingForce);
	p1->addMoment( ( (point1+point2)/2-p1->getCurrPosition() ) * (-CntDampingForce) );
	p2->addMoment( ( (point1+point2)/2-p2->getCurrPosition() ) * CntDampingForce );

	if (FRICTION != 0) {
	    // obtain tangential force
	    G0  = YOUNG/2/(1+POISSON);              // RelaDispInc points along point1's displacement relative to point2
	    vec RelaDispInc  = (veloc1-veloc2)*TIMESTEP;
	    vec TgtDispInc = RelaDispInc-(RelaDispInc%NormDirc)*NormDirc;
	    TgtDisp        = PreTgtDisp + TgtDispInc; // PreTgtDisp read by checkinPreTgt()
	    if (vfabs(TgtDisp) == 0)
		TgtDirc = 0;
	    else
		TgtDirc = normalize(-TgtDisp); // TgtDirc points along Tgtential forces exerted on particle 1

	    REAL fP = 0;
	    REAL ks = 0;

	    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	    // linear friction model
	    fP = FRICTION*vfabs(NormalForce);
	    ks = 4*G0*contact_radius/(2-POISSON);
	    TgtForce = PreTgtForce + ks*(-TgtDispInc); // PreTgtForce read by CheckinPreTgt()
	    if (vfabs(TgtForce) > fP)
		TgtForce = fP*TgtDirc;
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	    
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	    // Mindlin's model (loading/unloading condition assumed)
	    // This model is not recommended as it is impossible to strictly determine loading/unloading condition
            // unless load is known (the case of pure moment rotation).
#ifdef MINDLIN_ASSUMED
	    REAL val = 0;
	    fP = FRICTION*vfabs(NormalForce);
	    TgtLoading = (PreTgtDisp%TgtDispInc >= 0); 
	    
	    if (TgtLoading) {              // loading
		if (!PreTgtLoading) {      // pre-step is unloading
		    val = 8*G0*contact_radius*vfabs(TgtDispInc)/(3*(2-POISSON)*fP);
		    TgtDispStart = PreTgtDisp;
		}
		else                       // pre-step is loading
		    val = 8*G0*contact_radius*vfabs(TgtDisp-TgtDispStart)/(3*(2-POISSON)*fP);
		
		if (val > 1.0)              
		    TgtForce = fP*TgtDirc;
		else {
		    ks = 4*G0*contact_radius/(2-POISSON)*sqrt(1-val);
		    //incremental method
		    TgtForce = PreTgtForce + ks*(-TgtDispInc); // TgtDispInc determines signs
		    //total value method: TgtForce = fP*(1-pow(1-val, 1.5))*TgtDirc;
		}
	    }
	    else {                         // unloading
		if (PreTgtLoading) {       // pre-step is loading
		    val = 8*G0*contact_radius*vfabs(TgtDisp-TgtDispStart)/(3*(2-POISSON)*fP);
		    TgtPeak = vfabs(PreTgtForce);
		}
		else                       // pre-step is unloading
		    val = 8*G0*contact_radius*vfabs(TgtDisp-TgtDispStart)/(3*(2-POISSON)*fP);
		
		if (val > 1.0 || TgtPeak > fP)  
		    TgtForce = fP*TgtDirc;
		else {
		    ks = 2*sqrt(2)*G0*contact_radius/(2-POISSON) * sqrt(1+pow(1-TgtPeak/fP,2.0/3.0)+val);
		    //incremental method
		    TgtForce = PreTgtForce + ks*(-TgtDispInc); // TgtDispInc determines signs
		    //total value method: TgtForce = (TgtPeak-2*fP*(1-sqrt(2)/4*pow(1+ pow(1-TgtPeak/fP,2.0/3.0) + val,1.5)))*TgtDirc;
		}
	    }
	    
	    if (vfabs(TgtForce) > fP)
		TgtForce = fP*TgtDirc;
#endif
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////	

	    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	    // Mindlin's model (loading/unloading condition known for pure moment rotation case)
	    // As loading/unloading condition is known, both incremental and total value method work well.
            // Herein sliding history is incorporated.
#ifdef MINDLIN_KNOWN
	    REAL val = 0;
	    fP = FRICTION*vfabs(NormalForce);
	    if (PreTgtSlide)
		val = 8*G0*contact_radius*vfabs(TgtDispInc)/(3*(2-POISSON)*fP);
	    else
		val = 8*G0*contact_radius*vfabs(TgtDisp-TgtDispStart)/(3*(2-POISSON)*fP);

	    if (g_iteration > 10000 && g_iteration < 11000) { // loading (and possible sliding)
		if (val > 1.0) {
		    TgtForce = fP*TgtDirc;
		    TgtSlide = true;
		}
		else {
		    if (!PreTgtSlide) {
			ks = 4*G0*contact_radius/(2-POISSON)*sqrt(1-val);
			TgtForce = PreTgtForce + ks*(-TgtDispInc); // TgtDispInc determines signs
			TgtSlide = false;
		    }
		    else {
			if (vfabs(TgtForce)>vfabs(PreTgtForce))
			    TgtSlide = true;
			else
			    TgtSlide = false;
		    }
		}
		TgtPeak = vfabs(TgtForce);
	    }
	    else { // (possible sliding and) unloading
		if (val > 1.0 || TgtPeak > fP) {  
		    TgtForce = fP*TgtDirc;
		    TgtSlide = true;
		}
		else {
		    if (!PreTgtSlide) {
			ks = 2*sqrt(2)*G0*contact_radius/(2-POISSON) * sqrt(1+pow(1-TgtPeak/fP,2.0/3.0)+val);
			TgtForce = PreTgtForce + ks*(-TgtDispInc); // TgtDispInc determines signs
			TgtSlide = false;
		    }
		    else {
			if (vfabs(TgtForce)>vfabs(PreTgtForce))
			    TgtSlide = true;
			else {
			    TgtSlide = false;
			    TgtDispStart=TgtDisp;
			}
		    }
		}
	    }
	    
	    /*
	    g_debuginf<<"contact.h: g_iteration="<g_iteration
		      <<" PreTgtSlide="<<PreTgtSlide
		      <<" TgtSlide="<<TgtSlide
		      <<" val="<<val
		      <<" ks="<<ks
		      <<" TgtDispInc.x="<<TgtDispInc.getx()
		      <<" PreTgtForce="<<vfabs(PreTgtForce)
		      <<" TgtForce"<<vfabs(TgtForce)
		      <<std::endl;
	    */

	    if (vfabs(TgtForce) > fP)
		TgtForce = fP*TgtDirc;
#endif	    
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////

	    // apply tangential force
	    p1->addForce(TgtForce);
	    p2->addForce(-TgtForce);
	    p1->addMoment( ((point1+point2)/2-p1->getCurrPosition())*TgtForce);
	    p2->addMoment(-((point1+point2)/2-p2->getCurrPosition())*TgtForce);
	}
	
    }
    else {
	isInContact=false;
	TgtLoading=false;
	TgtPeak=0;
	NormalForce=0;
	TgtForce=0;
	TgtDisp=0;    //total value
	NormDirc=0;
	TgtDirc=0;

	penetration=0;
	contact_radius=0;
	radius1=radius2=0;
	spin_res=0;
	p1->mres=0;p2->mres=0;
    }   

}


template<class T>
void contact<T>::print() const {
    std::cout<<p1->getID()<<' '<<p2->getID()<<std::endl;
    std::cout<<"radius1="<<radius1<<' '<<"radius2="<<radius2<<std::endl;
    std::cout<<"normal force=";
    NormalForce.print();
    PreNormalForce.print();
    std::cout<<"tangential force=";
    TgtForce.print();
    PreTgtForce.print();
}

} // namespace dem ends

#endif

