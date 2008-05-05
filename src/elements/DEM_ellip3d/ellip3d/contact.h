//  This is a template class, for which we have to include the implementation in the header file.
//  As we can't put using statement in a header file, we have to use std::something wherever we
//  need to refer anything from standard namespace.
#ifndef CONTACT_H
#define CONTACT_H

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "const.h"
#include "root6.h"

//#define MINDLIN_ASSUMED
//#define MINDLIN_KNOWN

namespace dem {

extern std::ofstream g_exceptioninf;
extern int g_iteration;
extern long double DMP_CNT;
extern long double TIMESTEP;
extern long double FRICTION;
extern long double COHESION;

class cnttgt{
public:
    int  ptcl1;
    int  ptcl2;
    vec  TgtForce;
    vec  TgtDisp;
    bool TgtLoading;
    vec  TgtDispStart;
    long double TgtPeak;
    bool TgtSlide;

    cnttgt();
    cnttgt(int _ptcl1, int _ptcl2, vec _tf, vec _td, bool _tl, vec _tds, long double _tp, bool _ts)
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
    long double getRadius1() const {return radius1;}
    long double getRadius2() const {return radius2;}
    long double getR0() const {return R0;}
    long double getE0() const {return E0;}
    
    bool isOverlapped();
    void ContactForce();         // calculate normal and tangential force of contact
    long double getNormalForce() const {return vfabsl(NormalForce);}
    long double getTgtForce()  const {return vfabsl(TgtForce);}
    long double getPenetration() const {return penetration;}
    long double getContactRadius() const {return contact_radius;}
    long double getTgtDisp() const {return vfabsl(TgtDisp);} // total value during a process of contact
    void checkoutTgt(std::vector<cnttgt>& CntTgtVec);
    void checkinPreTgt(std::vector<cnttgt>& CntTgtVec);
    void print() const;
    
 private:
    T*   p1;                     // particle 1
    T*   p2;                     // particle 2
    long double penetration;     // penetration
    long double contact_radius;  // radius of contact surface
    vec  point1;                 // point1 on particle 1, innermost to particle 2
    vec  point2;                 // point2 on particle 2, innermost to particle 1
    long double radius1;         // radius of osculating circles at point1
    long double radius2;         // radius of osculating circles at point2
    long double E0;              
    long double G0;
    long double R0;

    bool IsInContact;

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

    long double TgtPeak;       

    vec  spin_res;
};


template <class T>
contact<T>::contact(){
    p1=NULL;
    p2=NULL;
    IsInContact=false;
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
    IsInContact=false;
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
T* contact<T>::getP1() const {
    return p1;
}


template<class T>
T* contact<T>::getP2() const {
    return p2;
}


template<class T>
bool contact<T>::isOverlapped(){
    long double coef1[10],coef2[10];
    p1->getGlobCoef(coef1); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef2);    
    vec v[2];
    if (root6(coef1,coef2,v[0]) && root6(coef2,coef1,v[1])) // a strict detection method
	return true;
    else
	return false;
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
void contact<T>::ContactForce(){
    long double coef1[10];
    long double coef2[10];
    p1->getGlobCoef(coef1); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef2); 
    vec v[2];    
    if (root6(coef1,coef2,v[0]) && root6(coef2,coef1,v[1])) // a strict detection method
	IsInContact=true;
    else
	IsInContact=false;
    point1=v[1];
    point2=v[0]; 

    // As UPDATE_CNT may not be 1, the two particles may not overlap each other when this
    // member function is called or when ContactList is referenced, so we need to judge based
    // on IsInContact. If not overlapped, then the forces are zero, which is reasonable for
    // an "existing" contact. UPDATE_CNT is essential to finding new contacts.
    if (IsInContact) {
	// obtain normal force, using absolute equation instead of stiffness method
	p1->cntnum++;
	p2->cntnum++;
	radius1=p1->getRadius(point1);
	radius2=p2->getRadius(point2);
	R0=radius1*radius2/(radius1+radius2);
	vec dist=point1-point2;
	penetration=vfabsl(dist);
	contact_radius=sqrtl(penetration*R0);
	E0=0.5*YOUNG/(1-POISSON*POISSON);
	NormDirc=normalize(dist);         // NormDirc points out of particle 1
	NormalForce= -sqrtl(penetration*penetration*penetration)*sqrtl(R0)*4*E0/3* NormDirc; // NormalForce pointing to particle 1
	// powl(penetration, 1.5)

        // apply cohesion force
	CohesionForce=PI*(penetration*R0)*COHESION*NormDirc;
	p1->addForce(CohesionForce);
	p2->addForce(-CohesionForce);

	// apply normal force
	p1->addForce(NormalForce);
	p2->addForce(-NormalForce);
	p1->addMoment( ((point1+point2)/2-p1->getCurrPosition())*NormalForce);
	p2->addMoment(-((point1+point2)/2-p2->getCurrPosition())*NormalForce);	
/*
	g_exceptioninf<<std::setw(10)<<g_iteration
		      <<std::setw(16)<<penetration
		      <<std::setw(16)<<vfabsl(CohesionForce)
		      <<std::setw(16)<<vfabsl(NormalForce)
		      <<std::setw(16)<<g_iteration*TIMESTEP
		      <<std::endl;
*/
	// obtain normal damping force
	vec cp = (point1+point2)/2;        
	vec veloc1 = p1->getCurrVelocity() + p1->getCurrOmga()*(cp-p1->getCurrPosition());
	vec veloc2 = p2->getCurrVelocity() + p2->getCurrOmga()*(cp-p2->getCurrPosition());
	long double m1 = getP1()->getMass();
	long double m2 = getP2()->getMass();
	long double kn = powl(6*vfabsl(NormalForce)*R0*powl(E0,2),1.0/3.0);
	long double DMP_CRTC = 2*sqrtl(m1*m2/(m1+m2)*kn); // critical damping
	vec CntDampingForce  = DMP_CNT * DMP_CRTC * ((veloc1-veloc2)%NormDirc)*NormDirc;

	// apply normal damping force
	p1->addForce(-CntDampingForce);
	p2->addForce(CntDampingForce);

	if (FRICTION != 0) {
	    // obtain tangential force
	    G0  = YOUNG/2/(1+POISSON);              // RelaDispInc points along point1's displacement relative to point2
	    vec RelaDispInc  = (veloc1-veloc2)*TIMESTEP;
	    vec TgtDispInc = RelaDispInc-(RelaDispInc%NormDirc)*NormDirc;
	    TgtDisp        = PreTgtDisp + TgtDispInc; // PreTgtDisp read by checkinPreTgt()
	    if (vfabsl(TgtDisp) == 0)
		TgtDirc = 0;
	    else
		TgtDirc = normalize(-TgtDisp); // TgtDirc points along Tgtential forces exerted on particle 1

	    long double fP = 0;
	    long double ks = 0;

	    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	    // linear friction model
	    fP = FRICTION*vfabsl(NormalForce);
	    ks = 4*G0*contact_radius/(2-POISSON);
	    TgtForce = PreTgtForce + ks*(-TgtDispInc); // PreTgtForce read by CheckinPreTgt()
	    if (vfabsl(TgtForce) > fP)
		TgtForce = fP*TgtDirc;
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	    
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	    // Mindlin's model (loading/unloading condition assumed)
	    // This model is not recommended as it is impossible to strictly determine loading/unloading condition
            // unless load is known (the case of pure moment rotation).
#ifdef MINDLIN_ASSUMED
	    long double val = 0;
	    fP = FRICTION*vfabsl(NormalForce);
	    TgtLoading = (PreTgtDisp%TgtDispInc >= 0); 
	    
	    if (TgtLoading) {              // loading
		if (!PreTgtLoading) {      // pre-step is unloading
		    val = 8*G0*contact_radius*vfabsl(TgtDispInc)/(3*(2-POISSON)*fP);
		    TgtDispStart = PreTgtDisp;
		}
		else                       // pre-step is loading
		    val = 8*G0*contact_radius*vfabsl(TgtDisp-TgtDispStart)/(3*(2-POISSON)*fP);
		
		if (val > 1.0)              
		    TgtForce = fP*TgtDirc;
		else {
		    ks = 4*G0*contact_radius/(2-POISSON)*sqrtl(1-val);
		    //incremental method
		    TgtForce = PreTgtForce + ks*(-TgtDispInc); // TgtDispInc determines signs
		    //total value method: TgtForce = fP*(1-powl(1-val, 1.5))*TgtDirc;
		}
	    }
	    else {                         // unloading
		if (PreTgtLoading) {       // pre-step is loading
		    val = 8*G0*contact_radius*vfabsl(TgtDisp-TgtDispStart)/(3*(2-POISSON)*fP);
		    TgtPeak = vfabsl(PreTgtForce);
		}
		else                       // pre-step is unloading
		    val = 8*G0*contact_radius*vfabsl(TgtDisp-TgtDispStart)/(3*(2-POISSON)*fP);
		
		if (val > 1.0 || TgtPeak > fP)  
		    TgtForce = fP*TgtDirc;
		else {
		    ks = 2*sqrtl(2)*G0*contact_radius/(2-POISSON) * sqrtl(1+powl(1-TgtPeak/fP,2.0/3.0)+val);
		    //incremental method
		    TgtForce = PreTgtForce + ks*(-TgtDispInc); // TgtDispInc determines signs
		    //total value method: TgtForce = (TgtPeak-2*fP*(1-sqrtl(2)/4*powl(1+ powl(1-TgtPeak/fP,2.0/3.0) + val,1.5)))*TgtDirc;
		}
	    }
	    
	    if (vfabsl(TgtForce) > fP)
		TgtForce = fP*TgtDirc;
#endif
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////	

	    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	    // Mindlin's model (loading/unloading condition known for pure moment rotation case)
	    // As loading/unloading condition is known, both incremental and total value method work well.
            // Herein sliding history is incorporated.
#ifdef MINDLIN_KNOWN
	    long double val = 0;
	    fP = FRICTION*vfabsl(NormalForce);
	    if (PreTgtSlide)
		val = 8*G0*contact_radius*vfabsl(TgtDispInc)/(3*(2-POISSON)*fP);
	    else
		val = 8*G0*contact_radius*vfabsl(TgtDisp-TgtDispStart)/(3*(2-POISSON)*fP);

	    if (g_iteration > 10000 && g_iteration < 11000) { // loading (and possible sliding)
		if (val > 1.0) {
		    TgtForce = fP*TgtDirc;
		    TgtSlide = true;
		}
		else {
		    if (!PreTgtSlide) {
			ks = 4*G0*contact_radius/(2-POISSON)*sqrtl(1-val);
			TgtForce = PreTgtForce + ks*(-TgtDispInc); // TgtDispInc determines signs
			TgtSlide = false;
		    }
		    else {
			if (vfabsl(TgtForce)>vfabsl(PreTgtForce))
			    TgtSlide = true;
			else
			    TgtSlide = false;
		    }
		}
		TgtPeak = vfabsl(TgtForce);
	    }
	    else { // (possible sliding and) unloading
		if (val > 1.0 || TgtPeak > fP) {  
		    TgtForce = fP*TgtDirc;
		    TgtSlide = true;
		}
		else {
		    if (!PreTgtSlide) {
			ks = 2*sqrtl(2)*G0*contact_radius/(2-POISSON) * sqrtl(1+powl(1-TgtPeak/fP,2.0/3.0)+val);
			TgtForce = PreTgtForce + ks*(-TgtDispInc); // TgtDispInc determines signs
			TgtSlide = false;
		    }
		    else {
			if (vfabsl(TgtForce)>vfabsl(PreTgtForce))
			    TgtSlide = true;
			else {
			    TgtSlide = false;
			    TgtDispStart=TgtDisp;
			}
		    }
		}
	    }
	    /*
	    g_exceptioninf<<std::setw(10)<<g_iteration
			  <<std::setw(05)<<PreTgtSlide
			  <<std::setw(05)<<TgtSlide
			  <<std::setw(16)<<val
			  <<std::setw(16)<<ks
			  <<std::setw(16)<<TgtDispInc.getx()
			  <<std::setw(16)<<vfabsl(PreTgtForce)
			  <<std::setw(16)<<vfabsl(TgtForce)
			  <<std::endl;
	    */
	    if (vfabsl(TgtForce) > fP)
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
	IsInContact=false;
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
