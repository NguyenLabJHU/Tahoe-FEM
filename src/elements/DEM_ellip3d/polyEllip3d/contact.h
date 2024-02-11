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
#include <algorithm>

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
    void Hopkin_contact(vec& oct1, vec& oct2);	// pre-detect which octants in Hopkin's algorithm. August 27, 2013
    void contactForce(bool &exceed);         // calculate normal and tangential force of contact. August 19, 2013
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
bool contact<T>::isOverlapped(){	// August 21, 2013



    /////////////////////////////////////////////////////////////////
    /////////////////// pre-detection step1  ////////////////////////
    /////////////////////////////////////////////////////////////////

    // if all corner points of one particle are outside of the other particle,
    // then this possible contact pair cannot be in contact, we delete this contact pair,
    // March 27, 2014. See case (1) in the notes in page 101

    // local coordinates of corner points of the two particles
    vec p1_ptcl1_local1, p2_ptcl1_local1, p3_ptcl1_local1, p4_ptcl1_local1, p5_ptcl1_local1, p6_ptcl1_local1, p7_ptcl1_local1, p8_ptcl1_local1;
    vec p1_ptcl2_local2, p2_ptcl2_local2, p3_ptcl2_local2, p4_ptcl2_local2, p5_ptcl2_local2, p6_ptcl2_local2, p7_ptcl2_local2, p8_ptcl2_local2;

    // get local coordinates of particle 1
    REAL aplus_ptcl1 = p1->getAplus(); REAL aminus_ptcl1 = p1->getAminus(); 
    REAL bplus_ptcl1 = p1->getBplus(); REAL bminus_ptcl1 = p1->getBminus();  
    REAL cplus_ptcl1 = p1->getCplus(); REAL cminus_ptcl1 = p1->getCminus();

    p1_ptcl1_local1 = vec(aplus_ptcl1, bplus_ptcl1, cplus_ptcl1);
    p2_ptcl1_local1 = vec(-aminus_ptcl1, bplus_ptcl1, cplus_ptcl1);
    p3_ptcl1_local1 = vec(-aminus_ptcl1, -bminus_ptcl1, cplus_ptcl1);
    p4_ptcl1_local1 = vec(aplus_ptcl1, -bminus_ptcl1, cplus_ptcl1);
    p5_ptcl1_local1 = vec(aplus_ptcl1, bplus_ptcl1, -cminus_ptcl1);
    p6_ptcl1_local1 = vec(-aminus_ptcl1, bplus_ptcl1, -cminus_ptcl1);
    p7_ptcl1_local1 = vec(-aminus_ptcl1, -bminus_ptcl1, -cminus_ptcl1);
    p8_ptcl1_local1 = vec(aplus_ptcl1, -bminus_ptcl1, -cminus_ptcl1);

    // get local coordinates of particle 2
    REAL aplus_ptcl2 = p2->getAplus(); REAL aminus_ptcl2 = p2->getAminus(); 
    REAL bplus_ptcl2 = p2->getBplus(); REAL bminus_ptcl2 = p2->getBminus();  
    REAL cplus_ptcl2 = p2->getCplus(); REAL cminus_ptcl2 = p2->getCminus();

    p1_ptcl2_local2 = vec(aplus_ptcl2, bplus_ptcl2, cplus_ptcl2);
    p2_ptcl2_local2 = vec(-aminus_ptcl2, bplus_ptcl2, cplus_ptcl2);
    p3_ptcl2_local2 = vec(-aminus_ptcl2, -bminus_ptcl2, cplus_ptcl2);
    p4_ptcl2_local2 = vec(aplus_ptcl2, -bminus_ptcl2, cplus_ptcl2);
    p5_ptcl2_local2 = vec(aplus_ptcl2, bplus_ptcl2, -cminus_ptcl2);
    p6_ptcl2_local2 = vec(-aminus_ptcl2, bplus_ptcl2, -cminus_ptcl2);
    p7_ptcl2_local2 = vec(-aminus_ptcl2, -bminus_ptcl2, -cminus_ptcl2);
    p8_ptcl2_local2 = vec(aplus_ptcl2, -bminus_ptcl2, -cminus_ptcl2);

    // use particle 1 as the boundary particle and check if particle 2 is outside of particle 1 

    // change the 8 corner points of particle 2 to the local coordinates of particel 1
    vec p1_ptcl2_global = p2->globalVec(p1_ptcl2_local2) + p2->getCurrPosition();
    vec p2_ptcl2_global = p2->globalVec(p2_ptcl2_local2) + p2->getCurrPosition();
    vec p3_ptcl2_global = p2->globalVec(p3_ptcl2_local2) + p2->getCurrPosition();
    vec p4_ptcl2_global = p2->globalVec(p4_ptcl2_local2) + p2->getCurrPosition();
    vec p5_ptcl2_global = p2->globalVec(p5_ptcl2_local2) + p2->getCurrPosition();
    vec p6_ptcl2_global = p2->globalVec(p6_ptcl2_local2) + p2->getCurrPosition();
    vec p7_ptcl2_global = p2->globalVec(p7_ptcl2_local2) + p2->getCurrPosition();
    vec p8_ptcl2_global = p2->globalVec(p8_ptcl2_local2) + p2->getCurrPosition();

    vec p1_ptcl2_local1 = p1->localVec(p1_ptcl2_global - p1->getCurrPosition());
    vec p2_ptcl2_local1 = p1->localVec(p2_ptcl2_global - p1->getCurrPosition());
    vec p3_ptcl2_local1 = p1->localVec(p3_ptcl2_global - p1->getCurrPosition());
    vec p4_ptcl2_local1 = p1->localVec(p4_ptcl2_global - p1->getCurrPosition());
    vec p5_ptcl2_local1 = p1->localVec(p5_ptcl2_global - p1->getCurrPosition());
    vec p6_ptcl2_local1 = p1->localVec(p6_ptcl2_global - p1->getCurrPosition());
    vec p7_ptcl2_local1 = p1->localVec(p7_ptcl2_global - p1->getCurrPosition());
    vec p8_ptcl2_local1 = p1->localVec(p8_ptcl2_global - p1->getCurrPosition());



    // check if all 8 corner points of particle 2 are outside of particle 1
    // x+ boundary, i.e. a+ boundary
    if(p1_ptcl2_local1.getx()>aplus_ptcl1 && p2_ptcl2_local1.getx()>aplus_ptcl1
    && p3_ptcl2_local1.getx()>aplus_ptcl1 && p4_ptcl2_local1.getx()>aplus_ptcl1
    && p5_ptcl2_local1.getx()>aplus_ptcl1 && p6_ptcl2_local1.getx()>aplus_ptcl1
    && p7_ptcl2_local1.getx()>aplus_ptcl1 && p8_ptcl2_local1.getx()>aplus_ptcl1 )	// outside of x+ boundary
	return false;

    // x- boundary, i.e. a- boundary
    if(p1_ptcl2_local1.getx()<-aminus_ptcl1 && p2_ptcl2_local1.getx()<-aminus_ptcl1
    && p3_ptcl2_local1.getx()<-aminus_ptcl1 && p4_ptcl2_local1.getx()<-aminus_ptcl1
    && p5_ptcl2_local1.getx()<-aminus_ptcl1 && p6_ptcl2_local1.getx()<-aminus_ptcl1
    && p7_ptcl2_local1.getx()<-aminus_ptcl1 && p8_ptcl2_local1.getx()<-aminus_ptcl1 )	// outside of x- boundary
	return false;

    // y+ boundary, i.e. b+ boundary
    if(p1_ptcl2_local1.gety()>bplus_ptcl1 && p2_ptcl2_local1.gety()>bplus_ptcl1
    && p3_ptcl2_local1.gety()>bplus_ptcl1 && p4_ptcl2_local1.gety()>bplus_ptcl1
    && p5_ptcl2_local1.gety()>bplus_ptcl1 && p6_ptcl2_local1.gety()>bplus_ptcl1
    && p7_ptcl2_local1.gety()>bplus_ptcl1 && p8_ptcl2_local1.gety()>bplus_ptcl1 )	// outside of y+ boundary
	return false;

    // y- boundary, i.e. b- boundary
    if(p1_ptcl2_local1.gety()<-bminus_ptcl1 && p2_ptcl2_local1.gety()<-bminus_ptcl1
    && p3_ptcl2_local1.gety()<-bminus_ptcl1 && p4_ptcl2_local1.gety()<-bminus_ptcl1
    && p5_ptcl2_local1.gety()<-bminus_ptcl1 && p6_ptcl2_local1.gety()<-bminus_ptcl1
    && p7_ptcl2_local1.gety()<-bminus_ptcl1 && p8_ptcl2_local1.gety()<-bminus_ptcl1 )	// outside of y- boundary
	return false;

    // z+ boundary, i.e. c+ boundary
    if(p1_ptcl2_local1.getz()>cplus_ptcl1 && p2_ptcl2_local1.getz()>cplus_ptcl1
    && p3_ptcl2_local1.getz()>cplus_ptcl1 && p4_ptcl2_local1.getz()>cplus_ptcl1
    && p5_ptcl2_local1.getz()>cplus_ptcl1 && p6_ptcl2_local1.getz()>cplus_ptcl1
    && p7_ptcl2_local1.getz()>cplus_ptcl1 && p8_ptcl2_local1.getz()>cplus_ptcl1 )	// outside of z+ boundary
	return false;

    // z- boundary, i.e. c- boundary
    if(p1_ptcl2_local1.getz()<-cminus_ptcl1 && p2_ptcl2_local1.getz()<-cminus_ptcl1
    && p3_ptcl2_local1.getz()<-cminus_ptcl1 && p4_ptcl2_local1.getz()<-cminus_ptcl1
    && p5_ptcl2_local1.getz()<-cminus_ptcl1 && p6_ptcl2_local1.getz()<-cminus_ptcl1
    && p7_ptcl2_local1.getz()<-cminus_ptcl1 && p8_ptcl2_local1.getz()<-cminus_ptcl1 )	// outside of z- boundary
	return false;



    /////////////////////////////////////////////////////////////////////////////////////////////
    // use particle 2 as the boundary particle and check if particle 1 is outside of particle 2 

    // change the 8 corner points of particle 1 to the local coordinates of particel 2
    vec p1_ptcl1_global = p1->globalVec(p1_ptcl1_local1) + p1->getCurrPosition();
    vec p2_ptcl1_global = p1->globalVec(p2_ptcl1_local1) + p1->getCurrPosition();
    vec p3_ptcl1_global = p1->globalVec(p3_ptcl1_local1) + p1->getCurrPosition();
    vec p4_ptcl1_global = p1->globalVec(p4_ptcl1_local1) + p1->getCurrPosition();
    vec p5_ptcl1_global = p1->globalVec(p5_ptcl1_local1) + p1->getCurrPosition();
    vec p6_ptcl1_global = p1->globalVec(p6_ptcl1_local1) + p1->getCurrPosition();
    vec p7_ptcl1_global = p1->globalVec(p7_ptcl1_local1) + p1->getCurrPosition();
    vec p8_ptcl1_global = p1->globalVec(p8_ptcl1_local1) + p1->getCurrPosition();

    vec p1_ptcl1_local2 = p2->localVec(p1_ptcl1_global - p2->getCurrPosition());
    vec p2_ptcl1_local2 = p2->localVec(p2_ptcl1_global - p2->getCurrPosition());
    vec p3_ptcl1_local2 = p2->localVec(p3_ptcl1_global - p2->getCurrPosition());
    vec p4_ptcl1_local2 = p2->localVec(p4_ptcl1_global - p2->getCurrPosition());
    vec p5_ptcl1_local2 = p2->localVec(p5_ptcl1_global - p2->getCurrPosition());
    vec p6_ptcl1_local2 = p2->localVec(p6_ptcl1_global - p2->getCurrPosition());
    vec p7_ptcl1_local2 = p2->localVec(p7_ptcl1_global - p2->getCurrPosition());
    vec p8_ptcl1_local2 = p2->localVec(p8_ptcl1_global - p2->getCurrPosition());


    // check if all 8 corner points of particle 1 are outside of particle 2
    // x+ boundary, i.e. a+ boundary
    if(p1_ptcl1_local2.getx()>aplus_ptcl2 && p2_ptcl1_local2.getx()>aplus_ptcl2
    && p3_ptcl1_local2.getx()>aplus_ptcl2 && p4_ptcl1_local2.getx()>aplus_ptcl2
    && p5_ptcl1_local2.getx()>aplus_ptcl2 && p6_ptcl1_local2.getx()>aplus_ptcl2
    && p7_ptcl1_local2.getx()>aplus_ptcl2 && p8_ptcl1_local2.getx()>aplus_ptcl2 )	// outside of x+ boundary
	return false;

    // x- boundary, i.e. a- boundary
    if(p1_ptcl1_local2.getx()<-aminus_ptcl2 && p2_ptcl1_local2.getx()<-aminus_ptcl2
    && p3_ptcl1_local2.getx()<-aminus_ptcl2 && p4_ptcl1_local2.getx()<-aminus_ptcl2
    && p5_ptcl1_local2.getx()<-aminus_ptcl2 && p6_ptcl1_local2.getx()<-aminus_ptcl2
    && p7_ptcl1_local2.getx()<-aminus_ptcl2 && p8_ptcl1_local2.getx()<-aminus_ptcl2 )	// outside of x- boundary
	return false;

    // y+ boundary, i.e. b+ boundary
    if(p1_ptcl1_local2.gety()>bplus_ptcl2 && p2_ptcl1_local2.gety()>bplus_ptcl2
    && p3_ptcl1_local2.gety()>bplus_ptcl2 && p4_ptcl1_local2.gety()>bplus_ptcl2
    && p5_ptcl1_local2.gety()>bplus_ptcl2 && p6_ptcl1_local2.gety()>bplus_ptcl2
    && p7_ptcl1_local2.gety()>bplus_ptcl2 && p8_ptcl1_local2.gety()>bplus_ptcl2 )	// outside of y+ boundary
	return false;

    // y- boundary, i.e. b- boundary
    if(p1_ptcl1_local2.gety()<-bminus_ptcl2 && p2_ptcl1_local2.gety()<-bminus_ptcl2
    && p3_ptcl1_local2.gety()<-bminus_ptcl2 && p4_ptcl1_local2.gety()<-bminus_ptcl2
    && p5_ptcl1_local2.gety()<-bminus_ptcl2 && p6_ptcl1_local2.gety()<-bminus_ptcl2
    && p7_ptcl1_local2.gety()<-bminus_ptcl2 && p8_ptcl1_local2.gety()<-bminus_ptcl2 )	// outside of y- boundary
	return false;

    // z+ boundary, i.e. c+ boundary
    if(p1_ptcl1_local2.getz()>cplus_ptcl2 && p2_ptcl1_local2.getz()>cplus_ptcl2
    && p3_ptcl1_local2.getz()>cplus_ptcl2 && p4_ptcl1_local2.getz()>cplus_ptcl2
    && p5_ptcl1_local2.getz()>cplus_ptcl2 && p6_ptcl1_local2.getz()>cplus_ptcl2
    && p7_ptcl1_local2.getz()>cplus_ptcl2 && p8_ptcl1_local2.getz()>cplus_ptcl2 )	// outside of z+ boundary
	return false;

    // z- boundary, i.e. c- boundary
    if(p1_ptcl1_local2.getz()<-cminus_ptcl2 && p2_ptcl1_local2.getz()<-cminus_ptcl2
    && p3_ptcl1_local2.getz()<-cminus_ptcl2 && p4_ptcl1_local2.getz()<-cminus_ptcl2
    && p5_ptcl1_local2.getz()<-cminus_ptcl2 && p6_ptcl1_local2.getz()<-cminus_ptcl2
    && p7_ptcl1_local2.getz()<-cminus_ptcl2 && p8_ptcl1_local2.getz()<-cminus_ptcl2 )	// outside of z- boundary
	return false;



    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     end of pre-detection step1      ///////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////



/*
    ////////////////////////////////////////////////////////////////////////////////////
    /////////   pre-detection step2 in the notes page 101, March 27, 2014     //////////
    ////////////////////////////////////////////////////////////////////////////////////

    std::vector<int> octant_ptcl1;	// store the possible octants of particle 1
    std::vector<int> octant_ptcl2;	// store the possible octants of particle 2

    // view particle 1 as boundary particle, and check which octants of particle 2 are contact with the six boundaries of particle 1

    // get the local coordinates of six points of particle 2 in local particle 1 system, they are centroid of the six faces
    vec aplus_ptcl2_local1  = 0.25*(p1_ptcl2_local1 + p4_ptcl2_local1 + p8_ptcl2_local1 + p5_ptcl2_local1);
    vec aminus_ptcl2_local1 = 0.25*(p2_ptcl2_local1 + p3_ptcl2_local1 + p7_ptcl2_local1 + p6_ptcl2_local1);
    vec bplus_ptcl2_local1  = 0.25*(p1_ptcl2_local1 + p2_ptcl2_local1 + p6_ptcl2_local1 + p5_ptcl2_local1);
    vec bminus_ptcl2_local1 = 0.25*(p4_ptcl2_local1 + p3_ptcl2_local1 + p7_ptcl2_local1 + p8_ptcl2_local1);
    vec cplus_ptcl2_local1  = 0.25*(p1_ptcl2_local1 + p2_ptcl2_local1 + p3_ptcl2_local1 + p4_ptcl2_local1);
    vec cminus_ptcl2_local1 = 0.25*(p5_ptcl2_local1 + p6_ptcl2_local1 + p7_ptcl2_local1 + p8_ptcl2_local1);

    // check with x+ boundary of particle 1
//    if(p1_ptcl2_local1.getx()<aplus_ptcl1 && p2_ptcl2_local1.getx()<aplus_ptcl1
//    && p3_ptcl2_local1.getx()<aplus_ptcl1 && p4_ptcl2_local1.getx()<aplus_ptcl1
//    && p5_ptcl2_local1.getx()<aplus_ptcl1 && p6_ptcl2_local1.getx()<aplus_ptcl1
//    && p7_ptcl2_local1.getx()<aplus_ptcl1 && p8_ptcl2_local1.getx()<aplus_ptcl1 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl2_local1.getx()<aminus_ptcl2_local1.getx())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl2_local1.getx()<bminus_ptcl2_local1.getx())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl2_local1.getx()<cminus_ptcl2_local1.getx())
		zsign = 1;
	else 
		zsign = -1;

	// store this octant number into octant_ptcl2;
	if(xsign==1 && ysign==1 && zsign==1)
	    octant_ptcl2.push_back(1);
	if(xsign==-1 && ysign==1 && zsign==1)
	    octant_ptcl2.push_back(2);
	if(xsign==-1 && ysign==-1 && zsign==1)
	    octant_ptcl2.push_back(3);
	if(xsign==1 && ysign==-1 && zsign==1)
	    octant_ptcl2.push_back(4);
	if(xsign==1 && ysign==1 && zsign==-1)
	    octant_ptcl2.push_back(5);
	if(xsign==-1 && ysign==1 && zsign==-1)
	    octant_ptcl2.push_back(6);
	if(xsign==-1 && ysign==-1 && zsign==-1)
	    octant_ptcl2.push_back(7);
	if(xsign==1 && ysign==-1 && zsign==-1)
	    octant_ptcl2.push_back(8);
    }


    // check with x- boundary of particle 1
//    if(p1_ptcl2_local1.getx()>-aminus_ptcl1 && p2_ptcl2_local1.getx()>-aminus_ptcl1
//    && p3_ptcl2_local1.getx()>-aminus_ptcl1 && p4_ptcl2_local1.getx()>-aminus_ptcl1
//    && p5_ptcl2_local1.getx()>-aminus_ptcl1 && p6_ptcl2_local1.getx()>-aminus_ptcl1
//    && p7_ptcl2_local1.getx()>-aminus_ptcl1 && p8_ptcl2_local1.getx()>-aminus_ptcl1 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl2_local1.getx()>aminus_ptcl2_local1.getx())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl2_local1.getx()>bminus_ptcl2_local1.getx())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl2_local1.getx()>cminus_ptcl2_local1.getx())
		zsign = 1;
	else 
		zsign = -1;

	// store this octant number into octant_ptcl2;
	if(xsign==1 && ysign==1 && zsign==1)
	    octant_ptcl2.push_back(1);
	if(xsign==-1 && ysign==1 && zsign==1)
	    octant_ptcl2.push_back(2);
	if(xsign==-1 && ysign==-1 && zsign==1)
	    octant_ptcl2.push_back(3);
	if(xsign==1 && ysign==-1 && zsign==1)
	    octant_ptcl2.push_back(4);
	if(xsign==1 && ysign==1 && zsign==-1)
	    octant_ptcl2.push_back(5);
	if(xsign==-1 && ysign==1 && zsign==-1)
	    octant_ptcl2.push_back(6);
	if(xsign==-1 && ysign==-1 && zsign==-1)
	    octant_ptcl2.push_back(7);
	if(xsign==1 && ysign==-1 && zsign==-1)
	    octant_ptcl2.push_back(8);
    }


    // check with y+ boundary of particle 1
//    if(p1_ptcl2_local1.gety()<bplus_ptcl1 && p2_ptcl2_local1.gety()<bplus_ptcl1
//    && p3_ptcl2_local1.gety()<bplus_ptcl1 && p4_ptcl2_local1.gety()<bplus_ptcl1
//    && p5_ptcl2_local1.gety()<bplus_ptcl1 && p6_ptcl2_local1.gety()<bplus_ptcl1
//    && p7_ptcl2_local1.gety()<bplus_ptcl1 && p8_ptcl2_local1.gety()<bplus_ptcl1 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl2_local1.gety()<aminus_ptcl2_local1.gety())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl2_local1.gety()<bminus_ptcl2_local1.gety())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl2_local1.gety()<cminus_ptcl2_local1.gety())
		zsign = 1;
	else 
		zsign = -1;

	// store this octant number into octant_ptcl2;
	if(xsign==1 && ysign==1 && zsign==1)
	    octant_ptcl2.push_back(1);
	if(xsign==-1 && ysign==1 && zsign==1)
	    octant_ptcl2.push_back(2);
	if(xsign==-1 && ysign==-1 && zsign==1)
	    octant_ptcl2.push_back(3);
	if(xsign==1 && ysign==-1 && zsign==1)
	    octant_ptcl2.push_back(4);
	if(xsign==1 && ysign==1 && zsign==-1)
	    octant_ptcl2.push_back(5);
	if(xsign==-1 && ysign==1 && zsign==-1)
	    octant_ptcl2.push_back(6);
	if(xsign==-1 && ysign==-1 && zsign==-1)
	    octant_ptcl2.push_back(7);
	if(xsign==1 && ysign==-1 && zsign==-1)
	    octant_ptcl2.push_back(8);
    }


    // check with y- boundary of particle 1
//    if(p1_ptcl2_local1.gety()>-bminus_ptcl1 && p2_ptcl2_local1.gety()>-bminus_ptcl1
//    && p3_ptcl2_local1.gety()>-bminus_ptcl1 && p4_ptcl2_local1.gety()>-bminus_ptcl1
//    && p5_ptcl2_local1.gety()>-bminus_ptcl1 && p6_ptcl2_local1.gety()>-bminus_ptcl1
//    && p7_ptcl2_local1.gety()>-bminus_ptcl1 && p8_ptcl2_local1.gety()>-bminus_ptcl1 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl2_local1.gety()>aminus_ptcl2_local1.gety())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl2_local1.gety()>bminus_ptcl2_local1.gety())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl2_local1.gety()>cminus_ptcl2_local1.gety())
		zsign = 1;
	else 
		zsign = -1;

	// store this octant number into octant_ptcl2;
	if(xsign==1 && ysign==1 && zsign==1)
	    octant_ptcl2.push_back(1);
	if(xsign==-1 && ysign==1 && zsign==1)
	    octant_ptcl2.push_back(2);
	if(xsign==-1 && ysign==-1 && zsign==1)
	    octant_ptcl2.push_back(3);
	if(xsign==1 && ysign==-1 && zsign==1)
	    octant_ptcl2.push_back(4);
	if(xsign==1 && ysign==1 && zsign==-1)
	    octant_ptcl2.push_back(5);
	if(xsign==-1 && ysign==1 && zsign==-1)
	    octant_ptcl2.push_back(6);
	if(xsign==-1 && ysign==-1 && zsign==-1)
	    octant_ptcl2.push_back(7);
	if(xsign==1 && ysign==-1 && zsign==-1)
	    octant_ptcl2.push_back(8);
    }


    // check with z+ boundary of particle 1
//    if(p1_ptcl2_local1.getz()<cplus_ptcl1 && p2_ptcl2_local1.getz()<cplus_ptcl1
//    && p3_ptcl2_local1.getz()<cplus_ptcl1 && p4_ptcl2_local1.getz()<cplus_ptcl1
//    && p5_ptcl2_local1.getz()<cplus_ptcl1 && p6_ptcl2_local1.getz()<cplus_ptcl1
//    && p7_ptcl2_local1.getz()<cplus_ptcl1 && p8_ptcl2_local1.getz()<cplus_ptcl1 )	// do nothing
//	;
//    else{
    {

	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl2_local1.getz()<aminus_ptcl2_local1.getz())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl2_local1.getz()<bminus_ptcl2_local1.getz())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl2_local1.getz()<cminus_ptcl2_local1.getz())
		zsign = 1;
	else 
		zsign = -1;

	// store this octant number into octant_ptcl2;
	if(xsign==1 && ysign==1 && zsign==1)
	    octant_ptcl2.push_back(1);
	if(xsign==-1 && ysign==1 && zsign==1)
	    octant_ptcl2.push_back(2);
	if(xsign==-1 && ysign==-1 && zsign==1)
	    octant_ptcl2.push_back(3);
	if(xsign==1 && ysign==-1 && zsign==1)
	    octant_ptcl2.push_back(4);
	if(xsign==1 && ysign==1 && zsign==-1)
	    octant_ptcl2.push_back(5);
	if(xsign==-1 && ysign==1 && zsign==-1)
	    octant_ptcl2.push_back(6);
	if(xsign==-1 && ysign==-1 && zsign==-1)
	    octant_ptcl2.push_back(7);
	if(xsign==1 && ysign==-1 && zsign==-1)
	    octant_ptcl2.push_back(8);
    }


    // check with z- boundary of particle 1
//    if(p1_ptcl2_local1.getz()>-cminus_ptcl1 && p2_ptcl2_local1.getz()>-cminus_ptcl1
//    && p3_ptcl2_local1.getz()>-cminus_ptcl1 && p4_ptcl2_local1.getz()>-cminus_ptcl1
//    && p5_ptcl2_local1.getz()>-cminus_ptcl1 && p6_ptcl2_local1.getz()>-cminus_ptcl1
//    && p7_ptcl2_local1.getz()>-cminus_ptcl1 && p8_ptcl2_local1.getz()>-cminus_ptcl1 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl2_local1.getz()>aminus_ptcl2_local1.getz())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl2_local1.getz()>bminus_ptcl2_local1.getz())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl2_local1.getz()>cminus_ptcl2_local1.getz())
		zsign = 1;
	else 
		zsign = -1;

	// store this octant number into octant_ptcl2;
	if(xsign==1 && ysign==1 && zsign==1)
	    octant_ptcl2.push_back(1);
	if(xsign==-1 && ysign==1 && zsign==1)
	    octant_ptcl2.push_back(2);
	if(xsign==-1 && ysign==-1 && zsign==1)
	    octant_ptcl2.push_back(3);
	if(xsign==1 && ysign==-1 && zsign==1)
	    octant_ptcl2.push_back(4);
	if(xsign==1 && ysign==1 && zsign==-1)
	    octant_ptcl2.push_back(5);
	if(xsign==-1 && ysign==1 && zsign==-1)
	    octant_ptcl2.push_back(6);
	if(xsign==-1 && ysign==-1 && zsign==-1)
	    octant_ptcl2.push_back(7);
	if(xsign==1 && ysign==-1 && zsign==-1)
	    octant_ptcl2.push_back(8);
    }




    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // view particle 2 as boundary particle, and check which octants of particle 1 are contact with the six boundaries of particle 2

    // get the local coordinates of six points of particle 2 in local particle 1 system, they are centroid of the six faces
    vec aplus_ptcl1_local2  = 0.25*(p1_ptcl1_local2 + p4_ptcl1_local2 + p8_ptcl1_local2 + p5_ptcl1_local2);
    vec aminus_ptcl1_local2 = 0.25*(p2_ptcl1_local2 + p3_ptcl1_local2 + p7_ptcl1_local2 + p6_ptcl1_local2);
    vec bplus_ptcl1_local2  = 0.25*(p1_ptcl1_local2 + p2_ptcl1_local2 + p6_ptcl1_local2 + p5_ptcl1_local2);
    vec bminus_ptcl1_local2 = 0.25*(p4_ptcl1_local2 + p3_ptcl1_local2 + p7_ptcl1_local2 + p8_ptcl1_local2);
    vec cplus_ptcl1_local2  = 0.25*(p1_ptcl1_local2 + p2_ptcl1_local2 + p3_ptcl1_local2 + p4_ptcl1_local2);
    vec cminus_ptcl1_local2 = 0.25*(p5_ptcl1_local2 + p6_ptcl1_local2 + p7_ptcl1_local2 + p8_ptcl1_local2);

    // check with x+ boundary of particle 2
//    if(p1_ptcl1_local2.getx()<aplus_ptcl2 && p2_ptcl1_local2.getx()<aplus_ptcl2
//    && p3_ptcl1_local2.getx()<aplus_ptcl2 && p4_ptcl1_local2.getx()<aplus_ptcl2
//    && p5_ptcl1_local2.getx()<aplus_ptcl2 && p6_ptcl1_local2.getx()<aplus_ptcl2
//    && p7_ptcl1_local2.getx()<aplus_ptcl2 && p8_ptcl1_local2.getx()<aplus_ptcl2 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl1_local2.getx()<aminus_ptcl1_local2.getx())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl1_local2.getx()<bminus_ptcl1_local2.getx())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl1_local2.getx()<cminus_ptcl1_local2.getx())
		zsign = 1;
	else 
		zsign = -1;

	// store this octant number into octant_ptcl1;
	if(xsign==1 && ysign==1 && zsign==1)
	    octant_ptcl1.push_back(1);
	if(xsign==-1 && ysign==1 && zsign==1)
	    octant_ptcl1.push_back(2);
	if(xsign==-1 && ysign==-1 && zsign==1)
	    octant_ptcl1.push_back(3);
	if(xsign==1 && ysign==-1 && zsign==1)
	    octant_ptcl1.push_back(4);
	if(xsign==1 && ysign==1 && zsign==-1)
	    octant_ptcl1.push_back(5);
	if(xsign==-1 && ysign==1 && zsign==-1)
	    octant_ptcl1.push_back(6);
	if(xsign==-1 && ysign==-1 && zsign==-1)
	    octant_ptcl1.push_back(7);
	if(xsign==1 && ysign==-1 && zsign==-1)
	    octant_ptcl1.push_back(8);
    }


    // check with x- boundary of particle 2
//    if(p1_ptcl1_local2.getx()>-aminus_ptcl2 && p2_ptcl1_local2.getx()>-aminus_ptcl2
//    && p3_ptcl1_local2.getx()>-aminus_ptcl2 && p4_ptcl1_local2.getx()>-aminus_ptcl2
//    && p5_ptcl1_local2.getx()>-aminus_ptcl2 && p6_ptcl1_local2.getx()>-aminus_ptcl2
//    && p7_ptcl1_local2.getx()>-aminus_ptcl2 && p8_ptcl1_local2.getx()>-aminus_ptcl2 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl1_local2.getx()>aminus_ptcl1_local2.getx())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl1_local2.getx()>bminus_ptcl1_local2.getx())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl1_local2.getx()>cminus_ptcl1_local2.getx())
		zsign = 1;
	else 
		zsign = -1;

	// store this octant number into octant_ptcl1;
	if(xsign==1 && ysign==1 && zsign==1)
	    octant_ptcl1.push_back(1);
	if(xsign==-1 && ysign==1 && zsign==1)
	    octant_ptcl1.push_back(2);
	if(xsign==-1 && ysign==-1 && zsign==1)
	    octant_ptcl1.push_back(3);
	if(xsign==1 && ysign==-1 && zsign==1)
	    octant_ptcl1.push_back(4);
	if(xsign==1 && ysign==1 && zsign==-1)
	    octant_ptcl1.push_back(5);
	if(xsign==-1 && ysign==1 && zsign==-1)
	    octant_ptcl1.push_back(6);
	if(xsign==-1 && ysign==-1 && zsign==-1)
	    octant_ptcl1.push_back(7);
	if(xsign==1 && ysign==-1 && zsign==-1)
	    octant_ptcl1.push_back(8);
    }


    // check with y+ boundary of particle 2
//    if(p1_ptcl1_local2.gety()<bplus_ptcl2 && p2_ptcl1_local2.gety()<bplus_ptcl2
//    && p3_ptcl1_local2.gety()<bplus_ptcl2 && p4_ptcl1_local2.gety()<bplus_ptcl2
//    && p5_ptcl1_local2.gety()<bplus_ptcl2 && p6_ptcl1_local2.gety()<bplus_ptcl2
//    && p7_ptcl1_local2.gety()<bplus_ptcl2 && p8_ptcl1_local2.gety()<bplus_ptcl2 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl1_local2.gety()<aminus_ptcl1_local2.gety())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl1_local2.gety()<bminus_ptcl1_local2.gety())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl1_local2.gety()<cminus_ptcl1_local2.gety())
		zsign = 1;
	else 
		zsign = -1;

	// store this octant number into octant_ptcl1;
	if(xsign==1 && ysign==1 && zsign==1)
	    octant_ptcl1.push_back(1);
	if(xsign==-1 && ysign==1 && zsign==1)
	    octant_ptcl1.push_back(2);
	if(xsign==-1 && ysign==-1 && zsign==1)
	    octant_ptcl1.push_back(3);
	if(xsign==1 && ysign==-1 && zsign==1)
	    octant_ptcl1.push_back(4);
	if(xsign==1 && ysign==1 && zsign==-1)
	    octant_ptcl1.push_back(5);
	if(xsign==-1 && ysign==1 && zsign==-1)
	    octant_ptcl1.push_back(6);
	if(xsign==-1 && ysign==-1 && zsign==-1)
	    octant_ptcl1.push_back(7);
	if(xsign==1 && ysign==-1 && zsign==-1)
	    octant_ptcl1.push_back(8);
    }


    // check with y- boundary of particle 2
//    if(p1_ptcl1_local2.gety()>-bminus_ptcl2 && p2_ptcl1_local2.gety()>-bminus_ptcl2
//    && p3_ptcl1_local2.gety()>-bminus_ptcl2 && p4_ptcl1_local2.gety()>-bminus_ptcl2
//    && p5_ptcl1_local2.gety()>-bminus_ptcl2 && p6_ptcl1_local2.gety()>-bminus_ptcl2
//    && p7_ptcl1_local2.gety()>-bminus_ptcl2 && p8_ptcl1_local2.gety()>-bminus_ptcl2 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl1_local2.gety()>aminus_ptcl1_local2.gety())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl1_local2.gety()>bminus_ptcl1_local2.gety())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl1_local2.gety()>cminus_ptcl1_local2.gety())
		zsign = 1;
	else 
		zsign = -1;

	// store this octant number into octant_ptcl1;
	if(xsign==1 && ysign==1 && zsign==1)
	    octant_ptcl1.push_back(1);
	if(xsign==-1 && ysign==1 && zsign==1)
	    octant_ptcl1.push_back(2);
	if(xsign==-1 && ysign==-1 && zsign==1)
	    octant_ptcl1.push_back(3);
	if(xsign==1 && ysign==-1 && zsign==1)
	    octant_ptcl1.push_back(4);
	if(xsign==1 && ysign==1 && zsign==-1)
	    octant_ptcl1.push_back(5);
	if(xsign==-1 && ysign==1 && zsign==-1)
	    octant_ptcl1.push_back(6);
	if(xsign==-1 && ysign==-1 && zsign==-1)
	    octant_ptcl1.push_back(7);
	if(xsign==1 && ysign==-1 && zsign==-1)
	    octant_ptcl1.push_back(8);
    }


    // check with z+ boundary of particle 2
//    if(p1_ptcl1_local2.getz()<cplus_ptcl2 && p2_ptcl1_local2.getz()<cplus_ptcl2
//    && p3_ptcl1_local2.getz()<cplus_ptcl2 && p4_ptcl1_local2.getz()<cplus_ptcl2
//    && p5_ptcl1_local2.getz()<cplus_ptcl2 && p6_ptcl1_local2.getz()<cplus_ptcl2
//    && p7_ptcl1_local2.getz()<cplus_ptcl2 && p8_ptcl1_local2.getz()<cplus_ptcl2 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl1_local2.getz()<aminus_ptcl1_local2.getz())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl1_local2.getz()<bminus_ptcl1_local2.getz())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl1_local2.getz()<cminus_ptcl1_local2.getz())
		zsign = 1;
	else 
		zsign = -1;

	// store this octant number into octant_ptcl1;
	if(xsign==1 && ysign==1 && zsign==1)
	    octant_ptcl1.push_back(1);
	if(xsign==-1 && ysign==1 && zsign==1)
	    octant_ptcl1.push_back(2);
	if(xsign==-1 && ysign==-1 && zsign==1)
	    octant_ptcl1.push_back(3);
	if(xsign==1 && ysign==-1 && zsign==1)
	    octant_ptcl1.push_back(4);
	if(xsign==1 && ysign==1 && zsign==-1)
	    octant_ptcl1.push_back(5);
	if(xsign==-1 && ysign==1 && zsign==-1)
	    octant_ptcl1.push_back(6);
	if(xsign==-1 && ysign==-1 && zsign==-1)
	    octant_ptcl1.push_back(7);
	if(xsign==1 && ysign==-1 && zsign==-1)
	    octant_ptcl1.push_back(8);
    }


    // check with z- boundary of particle 2
//    if(p1_ptcl1_local2.getz()>-cminus_ptcl2 && p2_ptcl1_local2.getz()>-cminus_ptcl2
//    && p3_ptcl1_local2.getz()>-cminus_ptcl2 && p4_ptcl1_local2.getz()>-cminus_ptcl2
//    && p5_ptcl1_local2.getz()>-cminus_ptcl2 && p6_ptcl1_local2.getz()>-cminus_ptcl2
//    && p7_ptcl1_local2.getz()>-cminus_ptcl2 && p8_ptcl1_local2.getz()>-cminus_ptcl2 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl1_local2.getz()>aminus_ptcl1_local2.getz())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl1_local2.getz()>bminus_ptcl1_local2.getz())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl1_local2.getz()>cminus_ptcl1_local2.getz())
		zsign = 1;
	else 
		zsign = -1;

	// store this octant number into octant_ptcl1;
	if(xsign==1 && ysign==1 && zsign==1)
	    octant_ptcl1.push_back(1);
	if(xsign==-1 && ysign==1 && zsign==1)
	    octant_ptcl1.push_back(2);
	if(xsign==-1 && ysign==-1 && zsign==1)
	    octant_ptcl1.push_back(3);
	if(xsign==1 && ysign==-1 && zsign==1)
	    octant_ptcl1.push_back(4);
	if(xsign==1 && ysign==1 && zsign==-1)
	    octant_ptcl1.push_back(5);
	if(xsign==-1 && ysign==1 && zsign==-1)
	    octant_ptcl1.push_back(6);
	if(xsign==-1 && ysign==-1 && zsign==-1)
	    octant_ptcl1.push_back(7);
	if(xsign==1 && ysign==-1 && zsign==-1)
	    octant_ptcl1.push_back(8);
    }


    //////////////////////////////////////////////////////////////////////
    /////////////////    end of pre-detection step2   ////////////////////
    //////////////////////////////////////////////////////////////////////

*/



//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////


/*

    ////////////////////////////////////////////////////////////////////////////////////
    /////////   pre-detection method 2 in the notes page 104, April 1, 2014     //////////
    ////////////////////////////////////////////////////////////////////////////////////

    // this is method is to check if any corner points of the particle box are inside of the other particle box
    // if so, we can detect possible octants, if not, the two particles are not in contact 


    // local coordinates of corner points of the two particles
    vec p1_ptcl1_local1, p2_ptcl1_local1, p3_ptcl1_local1, p4_ptcl1_local1, p5_ptcl1_local1, p6_ptcl1_local1, p7_ptcl1_local1, p8_ptcl1_local1;
    vec p1_ptcl2_local2, p2_ptcl2_local2, p3_ptcl2_local2, p4_ptcl2_local2, p5_ptcl2_local2, p6_ptcl2_local2, p7_ptcl2_local2, p8_ptcl2_local2;

    // get local coordinates of particle 1
    REAL aplus_ptcl1 = p1->getAplus(); REAL aminus_ptcl1 = p1->getAminus(); 
    REAL bplus_ptcl1 = p1->getBplus(); REAL bminus_ptcl1 = p1->getBminus();  
    REAL cplus_ptcl1 = p1->getCplus(); REAL cminus_ptcl1 = p1->getCminus();

    p1_ptcl1_local1 = vec(aplus_ptcl1, bplus_ptcl1, cplus_ptcl1);
    p2_ptcl1_local1 = vec(-aminus_ptcl1, bplus_ptcl1, cplus_ptcl1);
    p3_ptcl1_local1 = vec(-aminus_ptcl1, -bminus_ptcl1, cplus_ptcl1);
    p4_ptcl1_local1 = vec(aplus_ptcl1, -bminus_ptcl1, cplus_ptcl1);
    p5_ptcl1_local1 = vec(aplus_ptcl1, bplus_ptcl1, -cminus_ptcl1);
    p6_ptcl1_local1 = vec(-aminus_ptcl1, bplus_ptcl1, -cminus_ptcl1);
    p7_ptcl1_local1 = vec(-aminus_ptcl1, -bminus_ptcl1, -cminus_ptcl1);
    p8_ptcl1_local1 = vec(aplus_ptcl1, -bminus_ptcl1, -cminus_ptcl1);

    // get local coordinates of particle 2
    REAL aplus_ptcl2 = p2->getAplus(); REAL aminus_ptcl2 = p2->getAminus(); 
    REAL bplus_ptcl2 = p2->getBplus(); REAL bminus_ptcl2 = p2->getBminus();  
    REAL cplus_ptcl2 = p2->getCplus(); REAL cminus_ptcl2 = p2->getCminus();

    p1_ptcl2_local2 = vec(aplus_ptcl2, bplus_ptcl2, cplus_ptcl2);
    p2_ptcl2_local2 = vec(-aminus_ptcl2, bplus_ptcl2, cplus_ptcl2);
    p3_ptcl2_local2 = vec(-aminus_ptcl2, -bminus_ptcl2, cplus_ptcl2);
    p4_ptcl2_local2 = vec(aplus_ptcl2, -bminus_ptcl2, cplus_ptcl2);
    p5_ptcl2_local2 = vec(aplus_ptcl2, bplus_ptcl2, -cminus_ptcl2);
    p6_ptcl2_local2 = vec(-aminus_ptcl2, bplus_ptcl2, -cminus_ptcl2);
    p7_ptcl2_local2 = vec(-aminus_ptcl2, -bminus_ptcl2, -cminus_ptcl2);
    p8_ptcl2_local2 = vec(aplus_ptcl2, -bminus_ptcl2, -cminus_ptcl2);


    // change the 8 corner points of particle 2 to the local coordinates of particel 1
    vec p1_ptcl2_global = p2->globalVec(p1_ptcl2_local2) + p2->getCurrPosition();
    vec p2_ptcl2_global = p2->globalVec(p2_ptcl2_local2) + p2->getCurrPosition();
    vec p3_ptcl2_global = p2->globalVec(p3_ptcl2_local2) + p2->getCurrPosition();
    vec p4_ptcl2_global = p2->globalVec(p4_ptcl2_local2) + p2->getCurrPosition();
    vec p5_ptcl2_global = p2->globalVec(p5_ptcl2_local2) + p2->getCurrPosition();
    vec p6_ptcl2_global = p2->globalVec(p6_ptcl2_local2) + p2->getCurrPosition();
    vec p7_ptcl2_global = p2->globalVec(p7_ptcl2_local2) + p2->getCurrPosition();
    vec p8_ptcl2_global = p2->globalVec(p8_ptcl2_local2) + p2->getCurrPosition();

    vec p1_ptcl2_local1 = p1->localVec(p1_ptcl2_global - p1->getCurrPosition());
    vec p2_ptcl2_local1 = p1->localVec(p2_ptcl2_global - p1->getCurrPosition());
    vec p3_ptcl2_local1 = p1->localVec(p3_ptcl2_global - p1->getCurrPosition());
    vec p4_ptcl2_local1 = p1->localVec(p4_ptcl2_global - p1->getCurrPosition());
    vec p5_ptcl2_local1 = p1->localVec(p5_ptcl2_global - p1->getCurrPosition());
    vec p6_ptcl2_local1 = p1->localVec(p6_ptcl2_global - p1->getCurrPosition());
    vec p7_ptcl2_local1 = p1->localVec(p7_ptcl2_global - p1->getCurrPosition());
    vec p8_ptcl2_local1 = p1->localVec(p8_ptcl2_global - p1->getCurrPosition());
    
    vec center_ptcl2_local1 = 0.125*( p1_ptcl2_local1 + p2_ptcl2_local1 + p3_ptcl2_local1 + p4_ptcl2_local1
				    + p5_ptcl2_local1 + p6_ptcl2_local1 + p7_ptcl2_local1 + p8_ptcl2_local1 );

    // change the 8 corner points of particle 1 to the local coordinates of particel 2
    vec p1_ptcl1_global = p1->globalVec(p1_ptcl1_local1) + p1->getCurrPosition();
    vec p2_ptcl1_global = p1->globalVec(p2_ptcl1_local1) + p1->getCurrPosition();
    vec p3_ptcl1_global = p1->globalVec(p3_ptcl1_local1) + p1->getCurrPosition();
    vec p4_ptcl1_global = p1->globalVec(p4_ptcl1_local1) + p1->getCurrPosition();
    vec p5_ptcl1_global = p1->globalVec(p5_ptcl1_local1) + p1->getCurrPosition();
    vec p6_ptcl1_global = p1->globalVec(p6_ptcl1_local1) + p1->getCurrPosition();
    vec p7_ptcl1_global = p1->globalVec(p7_ptcl1_local1) + p1->getCurrPosition();
    vec p8_ptcl1_global = p1->globalVec(p8_ptcl1_local1) + p1->getCurrPosition();

    vec p1_ptcl1_local2 = p2->localVec(p1_ptcl1_global - p2->getCurrPosition());
    vec p2_ptcl1_local2 = p2->localVec(p2_ptcl1_global - p2->getCurrPosition());
    vec p3_ptcl1_local2 = p2->localVec(p3_ptcl1_global - p2->getCurrPosition());
    vec p4_ptcl1_local2 = p2->localVec(p4_ptcl1_global - p2->getCurrPosition());
    vec p5_ptcl1_local2 = p2->localVec(p5_ptcl1_global - p2->getCurrPosition());
    vec p6_ptcl1_local2 = p2->localVec(p6_ptcl1_global - p2->getCurrPosition());
    vec p7_ptcl1_local2 = p2->localVec(p7_ptcl1_global - p2->getCurrPosition());
    vec p8_ptcl1_local2 = p2->localVec(p8_ptcl1_global - p2->getCurrPosition());

    vec center_ptcl1_local2 = 0.125*( p1_ptcl1_local2 + p2_ptcl1_local2 + p3_ptcl1_local2 + p4_ptcl1_local2
				   + p5_ptcl1_local2 + p6_ptcl1_local2 + p7_ptcl1_local2 + p8_ptcl1_local2 );




    /////////////////////////////////////////////////////////////////////////////////////////
    ////////////////  method 2 step 1, as in notes page 104, April 1, 2014   ////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    // this step is to check if there are any corner points of the two particles are inside
    // of the other particle box, if so, then conduct 8x8 possible pairs, if not, return false
    int number_tmp = 0;	// the number of corner points that are inside of the boxes

    for(int num_oct = 1; num_oct != 10; num_oct++){
    
	REAL px, py, pz;

	// use particle 1 as the boundary particle box and check if 
	// corner points of particle 2 are inside of particle 1 box 
	if (num_oct == 1){
	    // check p1_ptcl2_local1
	    px = p1_ptcl2_local1.getx();
	    py = p1_ptcl2_local1.gety();
	    pz = p1_ptcl2_local1.getz();
 
	}
	else if (num_oct == 2){
	    // check p2_ptcl2_local1
	    px = p2_ptcl2_local1.getx();
	    py = p2_ptcl2_local1.gety();
	    pz = p2_ptcl2_local1.getz();
	}
	else if (num_oct == 3){
	    // check p3_ptcl2_local1
	    px = p3_ptcl2_local1.getx();
	    py = p3_ptcl2_local1.gety();
	    pz = p3_ptcl2_local1.getz();
	}
	else if (num_oct == 4){
	    // check p4_ptcl2_local1
	    px = p4_ptcl2_local1.getx();
	    py = p4_ptcl2_local1.gety();
	    pz = p4_ptcl2_local1.getz();
	}
	else if (num_oct == 5){
	    // check p5_ptcl2_local1
	    px = p5_ptcl2_local1.getx();
	    py = p5_ptcl2_local1.gety();
	    pz = p5_ptcl2_local1.getz();
	}
	else if (num_oct == 6){
	    // check p6_ptcl2_local1
	    px = p6_ptcl2_local1.getx();
	    py = p6_ptcl2_local1.gety();
	    pz = p6_ptcl2_local1.getz();
	}
	else if (num_oct == 7){
	    // check p7_ptcl2_local1
	    px = p7_ptcl2_local1.getx();
	    py = p7_ptcl2_local1.gety();
	    pz = p7_ptcl2_local1.getz();
	}
	else if (num_oct == 8){
	    // check p8_ptcl2_local1
	    px = p8_ptcl2_local1.getx();
	    py = p8_ptcl2_local1.gety();
	    pz = p8_ptcl2_local1.getz();
	}
	else{	
	    // chenck center_ptcl2_local1
	    px = center_ptcl2_local1.getx();
	    py = center_ptcl2_local1.gety();
	    pz = center_ptcl2_local1.getz();	
	}



	if( px <= aplus_ptcl1 && px >= -aminus_ptcl1 
	 && py <= bplus_ptcl1 && py >= -bminus_ptcl1 
	 && pz <= cplus_ptcl1 && pz >= -cminus_ptcl1 ){	// inside of particle 1 box

	    number_tmp++;
    	}


	// use particle 2 as the boundary particle box and check if 
	// corner points of particle 1 are inside of particle 2 box 
	if (num_oct == 1){
	    // check p1_ptcl1_local2
	    px = p1_ptcl1_local2.getx();
	    py = p1_ptcl1_local2.gety();
	    pz = p1_ptcl1_local2.getz();
 
	}
	else if (num_oct == 2){
	    // check p2_ptcl1_local2
	    px = p2_ptcl1_local2.getx();
	    py = p2_ptcl1_local2.gety();
	    pz = p2_ptcl1_local2.getz();
	}
	else if (num_oct == 3){
	    // check p3_ptcl1_local2
	    px = p3_ptcl1_local2.getx();
	    py = p3_ptcl1_local2.gety();
	    pz = p3_ptcl1_local2.getz();
	}
	else if (num_oct == 4){
	    // check p4_ptcl1_local2
	    px = p4_ptcl1_local2.getx();
	    py = p4_ptcl1_local2.gety();
	    pz = p4_ptcl1_local2.getz();
	}
	else if (num_oct == 5){
	    // check p5_ptcl1_local2
	    px = p5_ptcl1_local2.getx();
	    py = p5_ptcl1_local2.gety();
	    pz = p5_ptcl1_local2.getz();
	}
	else if (num_oct == 6){
	    // check p6_ptcl1_local2
	    px = p6_ptcl1_local2.getx();
	    py = p6_ptcl1_local2.gety();
	    pz = p6_ptcl1_local2.getz();
	}
	else if (num_oct == 7){
	    // check p7_ptcl1_local2
	    px = p7_ptcl1_local2.getx();
	    py = p7_ptcl1_local2.gety();
	    pz = p7_ptcl1_local2.getz();
	}
	else if (num_oct == 8){
	    // check p8_ptcl1_local2
	    px = p8_ptcl1_local2.getx();
	    py = p8_ptcl1_local2.gety();
	    pz = p8_ptcl1_local2.getz();
	}
	else{
	    // check center_ptcl1_local2
	    px = center_ptcl1_local2.getx();
	    py = center_ptcl1_local2.gety();
	    pz = center_ptcl1_local2.getz();
	}



	if( px <= aplus_ptcl2 && px >= -aminus_ptcl2 
	 && py <= bplus_ptcl2 && py >= -bminus_ptcl2 
	 && pz <= cplus_ptcl2 && pz >= -cminus_ptcl2 ){	// inside of particle 2 box

	    number_tmp++;
    	}

    } // end of num_oct loop

    if(number_tmp == 0)	// means the two particles cannot be contact
	return false;


    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////            end of method 2 step 1           ///////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////







    ////////////////////////////////////////////////////////////////////////////////////////
    //////////////   method 2 step 2, as in notes page 104, April 1, 2014    ///////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    std::vector<int> octant_ptcl1;	// store the possible octants of particle 1
    std::vector<int> octant_ptcl2;	// store the possible octants of particle 2


    for(int num_oct = 1; num_oct != 9; num_oct++){
    
	REAL px, py, pz;

	// use particle 1 as the boundary particle box and check if 
	// corner points of particle 2 are inside of particle 1 box 
	if (num_oct == 1){
	    // check p1_ptcl2_local1
	    px = p1_ptcl2_local1.getx();
	    py = p1_ptcl2_local1.gety();
	    pz = p1_ptcl2_local1.getz();
 
	}
	else if (num_oct == 2){
	    // check p2_ptcl2_local1
	    px = p2_ptcl2_local1.getx();
	    py = p2_ptcl2_local1.gety();
	    pz = p2_ptcl2_local1.getz();
	}
	else if (num_oct == 3){
	    // check p3_ptcl2_local1
	    px = p3_ptcl2_local1.getx();
	    py = p3_ptcl2_local1.gety();
	    pz = p3_ptcl2_local1.getz();
	}
	else if (num_oct == 4){
	    // check p4_ptcl2_local1
	    px = p4_ptcl2_local1.getx();
	    py = p4_ptcl2_local1.gety();
	    pz = p4_ptcl2_local1.getz();
	}
	else if (num_oct == 5){
	    // check p5_ptcl2_local1
	    px = p5_ptcl2_local1.getx();
	    py = p5_ptcl2_local1.gety();
	    pz = p5_ptcl2_local1.getz();
	}
	else if (num_oct == 6){
	    // check p6_ptcl2_local1
	    px = p6_ptcl2_local1.getx();
	    py = p6_ptcl2_local1.gety();
	    pz = p6_ptcl2_local1.getz();
	}
	else if (num_oct == 7){
	    // check p7_ptcl2_local1
	    px = p7_ptcl2_local1.getx();
	    py = p7_ptcl2_local1.gety();
	    pz = p7_ptcl2_local1.getz();
	}
	else{
	    // check p8_ptcl2_local1
	    px = p8_ptcl2_local1.getx();
	    py = p8_ptcl2_local1.gety();
	    pz = p8_ptcl2_local1.getz();
	}


	if( px <= aplus_ptcl1 && px >= -aminus_ptcl1 
	 && py <= bplus_ptcl1 && py >= -bminus_ptcl1 
	 && pz <= cplus_ptcl1 && pz >= -cminus_ptcl1 ){	// inside of particle 1 box

	    octant_ptcl2.push_back(num_oct);
	    // check in which octant of particle 1 box this point is
	    if( px <= aplus_ptcl1 && px >= 0
	     && py <= bplus_ptcl1 && py >= 0 
	     && pz <= cplus_ptcl1 && pz >= 0 ){	// in octant 1 of particle 1 box
		octant_ptcl1.push_back(1);
	    }
	    else if( px <= 0 && px >= -aminus_ptcl1 
		  && py <= bplus_ptcl1 && py >= 0
		  && pz <= cplus_ptcl1 && pz >= 0 ){	// in octant 2 of particle 1 box
		octant_ptcl1.push_back(2);
	    }
	    else if( px <= 0 && px >= -aminus_ptcl1 
		  && py <= 0 && py >= -bminus_ptcl1 
		  && pz <= cplus_ptcl1 && pz >= 0 ){	// in octant 3 of particle 1 box
		octant_ptcl1.push_back(3);
	    }
	    else if( px <= aplus_ptcl1 && px >= 0
		  && py <= 0 && py >= -bminus_ptcl1
		  && pz <= cplus_ptcl1 && px >= 0 ){	// in octant 4 of particle 1 box
		octant_ptcl1.push_back(4);
	    }
	    else if ( px <= aplus_ptcl1 && px >= 0
	     	   && py <= bplus_ptcl1 && py >= 0 
	     	   && pz <= 0 && pz >= -cminus_ptcl1 ){	// in octant 5 of particle 1 box
		octant_ptcl1.push_back(5);
	    }
	    else if( px <= 0 && px >= -aminus_ptcl1 
		  && py <= bplus_ptcl1 && py >= 0
		  && pz <= 0 && pz >= -cminus_ptcl1 ){	// in octant 6 of particle 1 box
		octant_ptcl1.push_back(6);
	    }
	    else if( px <= 0 && px >= -aminus_ptcl1 
		  && py <= 0 && py >= -bminus_ptcl1 
		  && pz <= 0 && pz >= -cminus_ptcl1 ){	// in octant 7 of particle 1 box
		octant_ptcl1.push_back(7);
	    }
	    else{	// in octant 8 of particle 1 box
		octant_ptcl1.push_back(8);
	    }

    	} // end if inside


	// use particle 2 as the boundary particle box and check if 
	// corner points of particle 1 are inside of particle 2 box 
	if (num_oct == 1){
	    // check p1_ptcl1_local2
	    px = p1_ptcl1_local2.getx();
	    py = p1_ptcl1_local2.gety();
	    pz = p1_ptcl1_local2.getz();
 
	}
	else if (num_oct == 2){
	    // check p2_ptcl1_local2
	    px = p2_ptcl1_local2.getx();
	    py = p2_ptcl1_local2.gety();
	    pz = p2_ptcl1_local2.getz();
	}
	else if (num_oct == 3){
	    // check p3_ptcl1_local2
	    px = p3_ptcl1_local2.getx();
	    py = p3_ptcl1_local2.gety();
	    pz = p3_ptcl1_local2.getz();
	}
	else if (num_oct == 4){
	    // check p4_ptcl1_local2
	    px = p4_ptcl1_local2.getx();
	    py = p4_ptcl1_local2.gety();
	    pz = p4_ptcl1_local2.getz();
	}
	else if (num_oct == 5){
	    // check p5_ptcl1_local2
	    px = p5_ptcl1_local2.getx();
	    py = p5_ptcl1_local2.gety();
	    pz = p5_ptcl1_local2.getz();
	}
	else if (num_oct == 6){
	    // check p6_ptcl1_local2
	    px = p6_ptcl1_local2.getx();
	    py = p6_ptcl1_local2.gety();
	    pz = p6_ptcl1_local2.getz();
	}
	else if (num_oct == 7){
	    // check p7_ptcl1_local2
	    px = p7_ptcl1_local2.getx();
	    py = p7_ptcl1_local2.gety();
	    pz = p7_ptcl1_local2.getz();
	}
	else{
	    // check p8_ptcl1_local2
	    px = p8_ptcl1_local2.getx();
	    py = p8_ptcl1_local2.gety();
	    pz = p8_ptcl1_local2.getz();
	}


	if( px <= aplus_ptcl2 && px >= -aminus_ptcl2 
	 && py <= bplus_ptcl2 && py >= -bminus_ptcl2 
	 && pz <= cplus_ptcl2 && pz >= -cminus_ptcl2 ){	// inside of particle 2 box

	    octant_ptcl1.push_back(num_oct);
	    // check in which octant of particle 2 box this point is
	    if( px <= aplus_ptcl2 && px >= 0
	     && py <= bplus_ptcl2 && py >= 0 
	     && pz <= cplus_ptcl2 && pz >= 0 ){	// in octant 1 of particle 2 box
		octant_ptcl2.push_back(1);
	    }
	    else if( px <= 0 && px >= -aminus_ptcl2 
		  && py <= bplus_ptcl2 && py >= 0
		  && pz <= cplus_ptcl2 && pz >= 0 ){	// in octant 2 of particle 2 box
		octant_ptcl2.push_back(2);
	    }
	    else if( px <= 0 && px >= -aminus_ptcl2 
		  && py <= 0 && py >= -bminus_ptcl2 
		  && pz <= cplus_ptcl2 && pz >= 0 ){	// in octant 3 of particle 2 box
		octant_ptcl2.push_back(3);
	    }
	    else if( px <= aplus_ptcl2 && px >= 0
		  && py <= 0 && py >= -bminus_ptcl2
		  && pz <= cplus_ptcl2 && px >= 0 ){	// in octant 4 of particle 2 box
		octant_ptcl2.push_back(4);
	    }
	    else if ( px <= aplus_ptcl2 && px >= 0
	     	   && py <= bplus_ptcl2 && py >= 0 
	     	   && pz <= 0 && pz >= -cminus_ptcl2 ){	// in octant 5 of particle 2 box
		octant_ptcl2.push_back(5);
	    }
	    else if( px <= 0 && px >= -aminus_ptcl2 
		  && py <= bplus_ptcl2 && py >= 0
		  && pz <= 0 && pz >= -cminus_ptcl2 ){	// in octant 6 of particle 2 box
		octant_ptcl2.push_back(6);
	    }
	    else if( px <= 0 && px >= -aminus_ptcl2 
		  && py <= 0 && py >= -bminus_ptcl2 
		  && pz <= 0 && pz >= -cminus_ptcl2 ){	// in octant 7 of particle 2 box
		octant_ptcl2.push_back(7);
	    }
	    else{	// in octant 8 of particle 1 box
		octant_ptcl2.push_back(8);
	    }

    	} // end if inside

    } // end of num_oct loop

    // check center_ptcl2_local1
    if( center_ptcl2_local1.getx() <= aplus_ptcl1 && center_ptcl2_local1.getx() >= -aminus_ptcl1 
     && center_ptcl2_local1.gety() <= bplus_ptcl1 && center_ptcl2_local1.gety() >= -bminus_ptcl1 
     && center_ptcl2_local1.getz() <= cplus_ptcl1 && center_ptcl2_local1.getz() >= -cminus_ptcl1 ){	// inside of particle 1 box

	for(int ii = 1; ii !=9; ii++){
	    octant_ptcl1.push_back(ii);
	    octant_ptcl2.push_back(ii); 
	}
    }

    // check center_ptcl1_local2
    if( center_ptcl1_local2.getx() <= aplus_ptcl2 && center_ptcl1_local2.getx() >= -aminus_ptcl2 
     && center_ptcl1_local2.gety() <= bplus_ptcl2 && center_ptcl1_local2.gety() >= -bminus_ptcl2 
     && center_ptcl1_local2.getz() <= cplus_ptcl2 && center_ptcl1_local2.getz() >= -cminus_ptcl2 ){	// inside of particle 2 box

	for(int ii = 1; ii !=9; ii++){
	    octant_ptcl1.push_back(ii);
	    octant_ptcl2.push_back(ii); 
	}
    }




    // delete dublicated octant number in the octant vectors
    // using default comparison:
    std::vector<int>::iterator it_oct;
    it_oct = std::unique (octant_ptcl1.begin(), octant_ptcl1.end());   
    octant_ptcl1.resize( std::distance(octant_ptcl1.begin(),it_oct) ); 

    it_oct = std::unique (octant_ptcl2.begin(), octant_ptcl2.end());   
    octant_ptcl2.resize( std::distance(octant_ptcl2.begin(),it_oct) ); 



    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////            end of method 2 step 2        /////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////


*/



/*

    int xsign_p1, ysign_p1, zsign_p1, xsign_p2, ysign_p2, zsign_p2;	// octants that p1 and p2 should use

    xsign_p1=0; ysign_p1=0; zsign_p1=0; xsign_p2=0; ysign_p2=0; zsign_p2=0;
    REAL coef_p1[10],coef_p2[10];

    vec v[2];
    bool b1, b2;
    vec local1_point1, local2_point2;



   for(unsigned i=0; i<octant_ptcl1.size(); i++){
	for(unsigned j=0; j<octant_ptcl2.size(); j++){

    // calculate contact points for 8x8 different cases

if(octant_ptcl1[i] == 1 && octant_ptcl2[j] == 1)    // (1,1)
{
    p1->getGlobCoef(coef_p1, 1); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 1);
    xsign_p1=1; ysign_p1=1; zsign_p1=1; xsign_p2=1; ysign_p2=1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() >= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() >= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 1 && octant_ptcl2[j] == 2)    // (1,2)
{
    p1->getGlobCoef(coef_p1, 1); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 2);
    xsign_p1=1; ysign_p1=1; zsign_p1=1; xsign_p2=-1; ysign_p2=1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() >= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() >= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 1 && octant_ptcl2[j] == 3)    // (1,3)
{
    p1->getGlobCoef(coef_p1, 1); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 3);
    xsign_p1=1; ysign_p1=1; zsign_p1=1; xsign_p2=-1; ysign_p2=-1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() >= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() <= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 1 && octant_ptcl2[j] == 4)    // (1,4)
{
    p1->getGlobCoef(coef_p1, 1); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 4);
    xsign_p1=1; ysign_p1=1; zsign_p1=1; xsign_p2=1; ysign_p2=-1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() >= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() <= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 1 && octant_ptcl2[j] == 5)    // (1,5)
{
    p1->getGlobCoef(coef_p1, 1); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 5);
    xsign_p1=1; ysign_p1=1; zsign_p1=1; xsign_p2=1; ysign_p2=1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() >= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() >= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 1 && octant_ptcl2[j] == 6)    // (1,6)
{
    p1->getGlobCoef(coef_p1, 1); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 6);
    xsign_p1=1; ysign_p1=1; zsign_p1=1; xsign_p2=-1; ysign_p2=1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() >= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() >= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 1 && octant_ptcl2[j] == 7)    // (1,7)
{
    p1->getGlobCoef(coef_p1, 1); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 7);
    xsign_p1=1; ysign_p1=1; zsign_p1=1; xsign_p2=-1; ysign_p2=-1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() >= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() <= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 1 && octant_ptcl2[j] == 8)    // (1,8)
{
    p1->getGlobCoef(coef_p1, 1); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 8);
    xsign_p1=1; ysign_p1=1; zsign_p1=1; xsign_p2=1; ysign_p2=-1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() >= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() <= 0 && local2_point2.getz() <= 0 )
	goto end;
}

/////////////////////////////////////////////////////////////////////////////////

else if(octant_ptcl1[i] == 2 && octant_ptcl2[j] == 1)    // (2,1)
{
    p1->getGlobCoef(coef_p1, 2); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 1);
    xsign_p1=-1; ysign_p1=1; zsign_p1=1; xsign_p2=1; ysign_p2=1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() >= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() >= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 2 && octant_ptcl2[j] == 2)    // (2,2)
{
    p1->getGlobCoef(coef_p1, 2); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 2);
    xsign_p1=-1; ysign_p1=1; zsign_p1=1; xsign_p2=-1; ysign_p2=1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() >= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() >= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 2 && octant_ptcl2[j] == 3)    // (2,3)
{
    p1->getGlobCoef(coef_p1, 2); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 3);
    xsign_p1=-1; ysign_p1=1; zsign_p1=1; xsign_p2=-1; ysign_p2=-1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() >= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() <= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 2 && octant_ptcl2[j] == 4)    // (2,4)
{
    p1->getGlobCoef(coef_p1, 2); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 4);
    xsign_p1=-1; ysign_p1=1; zsign_p1=1; xsign_p2=1; ysign_p2=-1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() >= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() <= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 2 && octant_ptcl2[j] == 5)    // (2,5)
{
    p1->getGlobCoef(coef_p1, 2); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 5);
    xsign_p1=-1; ysign_p1=1; zsign_p1=1; xsign_p2=1; ysign_p2=1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() >= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() >= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 2 && octant_ptcl2[j] == 6)    // (2,6)
{
    p1->getGlobCoef(coef_p1, 2); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 6);
    xsign_p1=-1; ysign_p1=1; zsign_p1=1; xsign_p2=-1; ysign_p2=1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() >= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() >= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 2 && octant_ptcl2[j] == 7)    // (2,7)
{
    p1->getGlobCoef(coef_p1, 2); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 7);
    xsign_p1=-1; ysign_p1=1; zsign_p1=1; xsign_p2=-1; ysign_p2=-1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() >= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() <= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 2 && octant_ptcl2[j] == 8)    // (2,8)
{
    p1->getGlobCoef(coef_p1, 2); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 8);
    xsign_p1=-1; ysign_p1=1; zsign_p1=1; xsign_p2=1; ysign_p2=-1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() >= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() <= 0 && local2_point2.getz() <= 0 )
	goto end;
}

////////////////////////////////////////////////////////////////////////////////

else if(octant_ptcl1[i] == 3 && octant_ptcl2[j] == 1)    // (3,1)
{
    p1->getGlobCoef(coef_p1, 3); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 1);
    xsign_p1=-1; ysign_p1=-1; zsign_p1=1; xsign_p2=1; ysign_p2=1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() <= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() >= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 3 && octant_ptcl2[j] == 2)    // (3,2)
{
    p1->getGlobCoef(coef_p1, 3); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 2);
    xsign_p1=-1; ysign_p1=-1; zsign_p1=1; xsign_p2=-1; ysign_p2=1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() <= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() >= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 3 && octant_ptcl2[j] == 3)    // (3,3)
{
    p1->getGlobCoef(coef_p1, 3); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 3);
    xsign_p1=-1; ysign_p1=-1; zsign_p1=1; xsign_p2=-1; ysign_p2=-1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() <= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() <= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 3 && octant_ptcl2[j] == 4)    // (3,4)
{
    p1->getGlobCoef(coef_p1, 3); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 4);
    xsign_p1=-1; ysign_p1=-1; zsign_p1=1; xsign_p2=1; ysign_p2=-1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() <= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() <= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 3 && octant_ptcl2[j] == 5)    // (3,5)
{
    p1->getGlobCoef(coef_p1, 3); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 5);
    xsign_p1=-1; ysign_p1=-1; zsign_p1=1; xsign_p2=1; ysign_p2=1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() <= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() >= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 3 && octant_ptcl2[j] == 6)    // (3,6)
{
    p1->getGlobCoef(coef_p1, 3); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 6);
    xsign_p1=-1; ysign_p1=-1; zsign_p1=1; xsign_p2=-1; ysign_p2=1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() <= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() >= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 3 && octant_ptcl2[j] == 7)    // (3,7)
{
    p1->getGlobCoef(coef_p1, 3); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 7);
    xsign_p1=-1; ysign_p1=-1; zsign_p1=1; xsign_p2=-1; ysign_p2=-1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() <= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() <= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 3 && octant_ptcl2[j] == 8)    // (3,8)
{
    p1->getGlobCoef(coef_p1, 3); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 8);
    xsign_p1=-1; ysign_p1=-1; zsign_p1=1; xsign_p2=1; ysign_p2=-1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() <= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() <= 0 && local2_point2.getz() <= 0 )
	goto end;
}

/////////////////////////////////////////////////////////////////////////////////////////////////

else if(octant_ptcl1[i] == 4 && octant_ptcl2[j] == 1)    // (4,1)
{
    p1->getGlobCoef(coef_p1, 4); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 1);
    xsign_p1=1; ysign_p1=-1; zsign_p1=1; xsign_p2=1; ysign_p2=1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() <= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() >= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 4 && octant_ptcl2[j] == 2)    // (4,2)
{
    p1->getGlobCoef(coef_p1, 4); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 2);
    xsign_p1=1; ysign_p1=-1; zsign_p1=1; xsign_p2=-1; ysign_p2=1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() <= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() >= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 4 && octant_ptcl2[j] == 3)    // (4,3)
{
    p1->getGlobCoef(coef_p1, 4); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 3);
    xsign_p1=1; ysign_p1=-1; zsign_p1=1; xsign_p2=-1; ysign_p2=-1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() <= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() <= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 4 && octant_ptcl2[j] == 4)    // (4,4)
{
    p1->getGlobCoef(coef_p1, 4); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 4);
    xsign_p1=1; ysign_p1=-1; zsign_p1=1; xsign_p2=1; ysign_p2=-1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() <= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() <= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 4 && octant_ptcl2[j] == 5)    // (4,5)
{
    p1->getGlobCoef(coef_p1, 4); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 5);
    xsign_p1=1; ysign_p1=-1; zsign_p1=1; xsign_p2=1; ysign_p2=1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() <= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() >= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 4 && octant_ptcl2[j] == 6)    // (4,6)
{
    p1->getGlobCoef(coef_p1, 4); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 6);
    xsign_p1=1; ysign_p1=-1; zsign_p1=1; xsign_p2=-1; ysign_p2=1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() <= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() >= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 4 && octant_ptcl2[j] == 7)    // (4,7)
{
    p1->getGlobCoef(coef_p1, 4); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 7);
    xsign_p1=1; ysign_p1=-1; zsign_p1=1; xsign_p2=-1; ysign_p2=-1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() <= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() <= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 4 && octant_ptcl2[j] == 8)    // (4,8)
{
    p1->getGlobCoef(coef_p1, 4); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 8);
    xsign_p1=1; ysign_p1=-1; zsign_p1=1; xsign_p2=1; ysign_p2=-1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() <= 0 && local1_point1.getz() >= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() <= 0 && local2_point2.getz() <= 0 )
	goto end;
}

//////////////////////////////////////////////////////////////////////////////////////////

else if(octant_ptcl1[i] == 5 && octant_ptcl2[j] == 1)    // (5,1)
{
    p1->getGlobCoef(coef_p1, 5); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 1);
    xsign_p1=1; ysign_p1=1; zsign_p1=-1; xsign_p2=1; ysign_p2=1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() >= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() >= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 5 && octant_ptcl2[j] == 2)    // (5,2)
{
    p1->getGlobCoef(coef_p1, 5); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 2);
    xsign_p1=1; ysign_p1=1; zsign_p1=-1; xsign_p2=-1; ysign_p2=1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() >= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() >= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 5 && octant_ptcl2[j] == 3)    // (5,3)
{
    p1->getGlobCoef(coef_p1, 5); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 3);
    xsign_p1=1; ysign_p1=1; zsign_p1=-1; xsign_p2=-1; ysign_p2=-1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() >= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() <= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 5 && octant_ptcl2[j] == 4)    // (5,4)
{
    p1->getGlobCoef(coef_p1, 5); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 4);
    xsign_p1=1; ysign_p1=1; zsign_p1=-1; xsign_p2=1; ysign_p2=-1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() >= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() <= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 5 && octant_ptcl2[j] == 5)    // (5,5)
{
    p1->getGlobCoef(coef_p1, 5); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 5);
    xsign_p1=1; ysign_p1=1; zsign_p1=-1; xsign_p2=1; ysign_p2=1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() >= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() >= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 5 && octant_ptcl2[j] == 6)    // (5,6)
{
    p1->getGlobCoef(coef_p1, 5); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 6);
    xsign_p1=1; ysign_p1=1; zsign_p1=-1; xsign_p2=-1; ysign_p2=1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() >= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() >= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 5 && octant_ptcl2[j] == 7)    // (5,7)
{
    p1->getGlobCoef(coef_p1, 5); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 7);
    xsign_p1=1; ysign_p1=1; zsign_p1=-1; xsign_p2=-1; ysign_p2=-1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() >= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() <= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 5 && octant_ptcl2[j] == 8)    // (5,8)
{
    p1->getGlobCoef(coef_p1, 5); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 8);
    xsign_p1=1; ysign_p1=1; zsign_p1=-1; xsign_p2=1; ysign_p2=-1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() >= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() <= 0 && local2_point2.getz() <= 0 )
	goto end;
}

///////////////////////////////////////////////////////////////////////////////////////////////

else if(octant_ptcl1[i] == 6 && octant_ptcl2[j] == 1)    // (6,1)
{
    p1->getGlobCoef(coef_p1, 6); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 1);
    xsign_p1=-1; ysign_p1=1; zsign_p1=-1; xsign_p2=1; ysign_p2=1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() >= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() >= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 6 && octant_ptcl2[j] == 2)    // (6,2)
{
    p1->getGlobCoef(coef_p1, 6); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 2);
    xsign_p1=-1; ysign_p1=1; zsign_p1=-1; xsign_p2=-1; ysign_p2=1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() >= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() >= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 6 && octant_ptcl2[j] == 3)    // (6,3)
{
    p1->getGlobCoef(coef_p1, 6); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 3);
    xsign_p1=-1; ysign_p1=1; zsign_p1=-1; xsign_p2=-1; ysign_p2=-1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() >= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() <= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 6 && octant_ptcl2[j] == 4)    // (6,4)
{
    p1->getGlobCoef(coef_p1, 6); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 4);
    xsign_p1=-1; ysign_p1=1; zsign_p1=-1; xsign_p2=1; ysign_p2=-1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() >= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() <= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 6 && octant_ptcl2[j] == 5)    // (6,5)
{
    p1->getGlobCoef(coef_p1, 6); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 5);
    xsign_p1=-1; ysign_p1=1; zsign_p1=-1; xsign_p2=1; ysign_p2=1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() >= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() >= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 6 && octant_ptcl2[j] == 6)    // (6,6)
{
    p1->getGlobCoef(coef_p1, 6); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 6);
    xsign_p1=-1; ysign_p1=1; zsign_p1=-1; xsign_p2=-1; ysign_p2=1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() >= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() >= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 6 && octant_ptcl2[j] == 7)    // (6,7)
{
    p1->getGlobCoef(coef_p1, 6); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 7);
    xsign_p1=-1; ysign_p1=1; zsign_p1=-1; xsign_p2=-1; ysign_p2=-1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() >= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() <= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 6 && octant_ptcl2[j] == 8)    // (6,8)
{
    p1->getGlobCoef(coef_p1, 6); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 8);
    xsign_p1=-1; ysign_p1=1; zsign_p1=-1; xsign_p2=1; ysign_p2=-1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() >= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() <= 0 && local2_point2.getz() <= 0 )
	goto end;
}

////////////////////////////////////////////////////////////////////////////////////////////////

else if(octant_ptcl1[i] == 7 && octant_ptcl2[j] == 1)    // (7,1)
{
    p1->getGlobCoef(coef_p1, 7); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 1);
    xsign_p1=-1; ysign_p1=-1; zsign_p1=-1; xsign_p2=1; ysign_p2=1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() <= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() >= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 7 && octant_ptcl2[j] == 2)    // (7,2)
{
    p1->getGlobCoef(coef_p1, 7); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 2);
    xsign_p1=-1; ysign_p1=-1; zsign_p1=-1; xsign_p2=-1; ysign_p2=1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() <= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() >= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 7 && octant_ptcl2[j] == 3)    // (7,3)
{
    p1->getGlobCoef(coef_p1, 7); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 3);
    xsign_p1=-1; ysign_p1=-1; zsign_p1=-1; xsign_p2=-1; ysign_p2=-1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() <= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() <= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 7 && octant_ptcl2[j] == 4)    // (7,4)
{
    p1->getGlobCoef(coef_p1, 7); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 4);
    xsign_p1=-1; ysign_p1=-1; zsign_p1=-1; xsign_p2=1; ysign_p2=-1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() <= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() <= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 7 && octant_ptcl2[j] == 5)    // (7,5)
{
    p1->getGlobCoef(coef_p1, 7); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 5);
    xsign_p1=-1; ysign_p1=-1; zsign_p1=-1; xsign_p2=1; ysign_p2=1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() <= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() >= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 7 && octant_ptcl2[j] == 6)    // (7,6)
{
    p1->getGlobCoef(coef_p1, 7); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 6);
    xsign_p1=-1; ysign_p1=-1; zsign_p1=-1; xsign_p2=-1; ysign_p2=1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() <= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() >= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 7 && octant_ptcl2[j] == 7)    // (7,7)
{
    p1->getGlobCoef(coef_p1, 7); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 7);
    xsign_p1=-1; ysign_p1=-1; zsign_p1=-1; xsign_p2=-1; ysign_p2=-1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() <= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() <= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 7 && octant_ptcl2[j] == 8)    // (7,8)
{
    p1->getGlobCoef(coef_p1, 7); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 8);
    xsign_p1=-1; ysign_p1=-1; zsign_p1=-1; xsign_p2=1; ysign_p2=-1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() <= 0 && local1_point1.gety() <= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() <= 0 && local2_point2.getz() <= 0 )
	goto end;
}

//////////////////////////////////////////////////////////////////////////////////////////////

else if(octant_ptcl1[i] == 8 && octant_ptcl2[j] == 1)    // (8,1)
{
    p1->getGlobCoef(coef_p1, 8); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 1);
    xsign_p1=1; ysign_p1=-1; zsign_p1=-1; xsign_p2=1; ysign_p2=1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() <= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() >= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 8 && octant_ptcl2[j] == 2)    // (8,2)
{
    p1->getGlobCoef(coef_p1, 8); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 2);
    xsign_p1=1; ysign_p1=-1; zsign_p1=-1; xsign_p2=-1; ysign_p2=1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() <= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() >= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 8 && octant_ptcl2[j] == 3)    // (8,3)
{
    p1->getGlobCoef(coef_p1, 8); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 3);
    xsign_p1=1; ysign_p1=-1; zsign_p1=-1; xsign_p2=-1; ysign_p2=-1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() <= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() <= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 8 && octant_ptcl2[j] == 4)    // (8,4)
{
    p1->getGlobCoef(coef_p1, 8); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 4);
    xsign_p1=1; ysign_p1=-1; zsign_p1=-1; xsign_p2=1; ysign_p2=-1; zsign_p2=1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() <= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() <= 0 && local2_point2.getz() >= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 8 && octant_ptcl2[j] == 5)    // (8,5)
{
    p1->getGlobCoef(coef_p1, 8); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 5);
    xsign_p1=1; ysign_p1=-1; zsign_p1=-1; xsign_p2=1; ysign_p2=1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() <= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() >= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 8 && octant_ptcl2[j] == 6)    // (8,6)
{
    p1->getGlobCoef(coef_p1, 8); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 6);
    xsign_p1=1; ysign_p1=-1; zsign_p1=-1; xsign_p2=-1; ysign_p2=1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() <= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() >= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 8 && octant_ptcl2[j] == 7)    // (8,7)
{
    p1->getGlobCoef(coef_p1, 8); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 7);
    xsign_p1=1; ysign_p1=-1; zsign_p1=-1; xsign_p2=-1; ysign_p2=-1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() <= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() <= 0 && local2_point2.gety() <= 0 && local2_point2.getz() <= 0 )
	goto end;
}

else if(octant_ptcl1[i] == 8 && octant_ptcl2[j] == 8)    // (8,8)
{
    p1->getGlobCoef(coef_p1, 8); // v[0] is the point on p2, v[1] is the point on p1
    p2->getGlobCoef(coef_p2, 8);
    xsign_p1=1; ysign_p1=-1; zsign_p1=-1; xsign_p2=1; ysign_p2=-1; zsign_p2=-1;   
    b1 = root6(coef_p1,coef_p2,v[0]);
    b2 = root6(coef_p2,coef_p1,v[1]);
    point1 = v[1];
    point2 = v[0];

    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    if(b1 && b2 
    && local1_point1.getx() >= 0 && local1_point1.gety() <= 0 && local1_point1.getz() <= 0
    && local2_point2.getx() >= 0 && local2_point2.gety() <= 0 && local2_point2.getz() <= 0 )
	goto end;
}


	}  // end of j loop
    }	// end of i loop

/////////////////////////////////////////////////////////////////////////////////////////////////




    // if no contact points are suitable
    return false;

end:



*/



//////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////   the step 2 of new method, as in notes page 106, April 2, 2014    //////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

    // get the 27 coordinates of the sub box points of particle 1 in local2
    REAL aplus_ratio_ptcl1 = aplus_ptcl1/(aplus_ptcl1+aminus_ptcl1);
    REAL bplus_ratio_ptcl1 = bplus_ptcl1/(bplus_ptcl1+bminus_ptcl1);
    REAL cplus_ratio_ptcl1 = cplus_ptcl1/(cplus_ptcl1+cminus_ptcl1);

    vec p9_ptcl1_local2  = p1_ptcl1_local2 + aplus_ratio_ptcl1*(p2_ptcl1_local2 - p1_ptcl1_local2);
    vec p10_ptcl1_local2 = p4_ptcl1_local2 + aplus_ratio_ptcl1*(p3_ptcl1_local2 - p4_ptcl1_local2);
    vec p11_ptcl1_local2 = p8_ptcl1_local2 + aplus_ratio_ptcl1*(p7_ptcl1_local2 - p8_ptcl1_local2);
    vec p12_ptcl1_local2 = p5_ptcl1_local2 + aplus_ratio_ptcl1*(p6_ptcl1_local2 - p5_ptcl1_local2);

    vec p13_ptcl1_local2 = p1_ptcl1_local2 + bplus_ratio_ptcl1*(p4_ptcl1_local2 - p1_ptcl1_local2);
    vec p14_ptcl1_local2 = p2_ptcl1_local2 + bplus_ratio_ptcl1*(p3_ptcl1_local2 - p2_ptcl1_local2);
    vec p15_ptcl1_local2 = p6_ptcl1_local2 + bplus_ratio_ptcl1*(p7_ptcl1_local2 - p6_ptcl1_local2);
    vec p16_ptcl1_local2 = p5_ptcl1_local2 + bplus_ratio_ptcl1*(p8_ptcl1_local2 - p5_ptcl1_local2);

    vec p17_ptcl1_local2 = p1_ptcl1_local2 + cplus_ratio_ptcl1*(p5_ptcl1_local2 - p1_ptcl1_local2);
    vec p18_ptcl1_local2 = p2_ptcl1_local2 + cplus_ratio_ptcl1*(p6_ptcl1_local2 - p2_ptcl1_local2);
    vec p19_ptcl1_local2 = p3_ptcl1_local2 + cplus_ratio_ptcl1*(p7_ptcl1_local2 - p3_ptcl1_local2);
    vec p20_ptcl1_local2 = p4_ptcl1_local2 + cplus_ratio_ptcl1*(p8_ptcl1_local2 - p4_ptcl1_local2);

    vec p21_ptcl1_local2 = p13_ptcl1_local2 + aplus_ratio_ptcl1*(p14_ptcl1_local2 - p13_ptcl1_local2);
    vec p22_ptcl1_local2 = p9_ptcl1_local2  + cplus_ratio_ptcl1*(p12_ptcl1_local2 -  p9_ptcl1_local2);
    vec p23_ptcl1_local2 = p14_ptcl1_local2 + cplus_ratio_ptcl1*(p15_ptcl1_local2 - p14_ptcl1_local2);
    vec p24_ptcl1_local2 = p10_ptcl1_local2 + cplus_ratio_ptcl1*(p11_ptcl1_local2 - p10_ptcl1_local2);
    vec p25_ptcl1_local2 = p13_ptcl1_local2 + cplus_ratio_ptcl1*(p16_ptcl1_local2 - p13_ptcl1_local2);
    vec p26_ptcl1_local2 = p16_ptcl1_local2 + aplus_ratio_ptcl1*(p15_ptcl1_local2 - p16_ptcl1_local2);
    vec p27_ptcl1_local2 = p21_ptcl1_local2 + cplus_ratio_ptcl1*(p26_ptcl1_local2 - p21_ptcl1_local2);

    // get the 27 coordinates of the sub box points of particle 2 in local1
    REAL aplus_ratio_ptcl2 = aplus_ptcl2/(aplus_ptcl2+aminus_ptcl2);
    REAL bplus_ratio_ptcl2 = bplus_ptcl2/(bplus_ptcl2+bminus_ptcl2);
    REAL cplus_ratio_ptcl2 = cplus_ptcl2/(cplus_ptcl2+cminus_ptcl2);

    vec p9_ptcl2_local1  = p1_ptcl2_local1 + aplus_ratio_ptcl2*(p2_ptcl2_local1 - p1_ptcl2_local1);
    vec p10_ptcl2_local1 = p4_ptcl2_local1 + aplus_ratio_ptcl2*(p3_ptcl2_local1 - p4_ptcl2_local1);
    vec p11_ptcl2_local1 = p8_ptcl2_local1 + aplus_ratio_ptcl2*(p7_ptcl2_local1 - p8_ptcl2_local1);
    vec p12_ptcl2_local1 = p5_ptcl2_local1 + aplus_ratio_ptcl2*(p6_ptcl2_local1 - p5_ptcl2_local1);

    vec p13_ptcl2_local1 = p1_ptcl2_local1 + bplus_ratio_ptcl2*(p4_ptcl2_local1 - p1_ptcl2_local1);
    vec p14_ptcl2_local1 = p2_ptcl2_local1 + bplus_ratio_ptcl2*(p3_ptcl2_local1 - p2_ptcl2_local1);
    vec p15_ptcl2_local1 = p6_ptcl2_local1 + bplus_ratio_ptcl2*(p7_ptcl2_local1 - p6_ptcl2_local1);
    vec p16_ptcl2_local1 = p5_ptcl2_local1 + bplus_ratio_ptcl2*(p8_ptcl2_local1 - p5_ptcl2_local1);

    vec p17_ptcl2_local1 = p1_ptcl2_local1 + cplus_ratio_ptcl2*(p5_ptcl2_local1 - p1_ptcl2_local1);
    vec p18_ptcl2_local1 = p2_ptcl2_local1 + cplus_ratio_ptcl2*(p6_ptcl2_local1 - p2_ptcl2_local1);
    vec p19_ptcl2_local1 = p3_ptcl2_local1 + cplus_ratio_ptcl2*(p7_ptcl2_local1 - p3_ptcl2_local1);
    vec p20_ptcl2_local1 = p4_ptcl2_local1 + cplus_ratio_ptcl2*(p8_ptcl2_local1 - p4_ptcl2_local1);

    vec p21_ptcl2_local1 = p13_ptcl2_local1 + aplus_ratio_ptcl2*(p14_ptcl2_local1 - p13_ptcl2_local1);
    vec p22_ptcl2_local1 = p9_ptcl2_local1  + cplus_ratio_ptcl2*(p12_ptcl2_local1 - p9_ptcl2_local1);
    vec p23_ptcl2_local1 = p14_ptcl2_local1 + cplus_ratio_ptcl2*(p15_ptcl2_local1 - p14_ptcl2_local1);
    vec p24_ptcl2_local1 = p10_ptcl2_local1 + cplus_ratio_ptcl2*(p11_ptcl2_local1 - p10_ptcl2_local1);
    vec p25_ptcl2_local1 = p13_ptcl2_local1 + cplus_ratio_ptcl2*(p16_ptcl2_local1 - p13_ptcl2_local1);
    vec p26_ptcl2_local1 = p16_ptcl2_local1 + aplus_ratio_ptcl2*(p15_ptcl2_local1 - p16_ptcl2_local1);
    vec p27_ptcl2_local1 = p21_ptcl2_local1 + cplus_ratio_ptcl2*(p26_ptcl2_local1 - p21_ptcl2_local1);




    int xsign_p1, ysign_p1, zsign_p1, xsign_p2, ysign_p2, zsign_p2;	// octants that p1 and p2 should use

    xsign_p1=0; ysign_p1=0; zsign_p1=0; xsign_p2=0; ysign_p2=0; zsign_p2=0;
    REAL coef_p1[10],coef_p2[10];

    vec v[2];
    bool b1, b2;
    vec local1_point1, local2_point2;

    bool flag_contact = false;	// the flag if the two polyellipsoids are contact
    for(int num_ptcl1 = 1; num_ptcl1 !=9; num_ptcl1++){
	for(int num_ptcl2 = 1; num_ptcl2 != 9; num_ptcl2++){

	    // CHECK IF SUB BOX PAIR IS IN CONTACT OR NOT

	    // use sub box of particle 1 as the boundary particle 
	    // and check if sub box of particle 2 is outside of this sub box
	    vec p1_sub_ptcl2_local1, p2_sub_ptcl2_local1, p3_sub_ptcl2_local1, p4_sub_ptcl2_local1;
	    vec p5_sub_ptcl2_local1, p6_sub_ptcl2_local1, p7_sub_ptcl2_local1, p8_sub_ptcl2_local1;
	    REAL aplus_sub_ptcl1, aminus_sub_ptcl1;
	    REAL bplus_sub_ptcl1, bminus_sub_ptcl1;
	    REAL cplus_sub_ptcl1, cminus_sub_ptcl1;

	    // use sub box of particle 2 as the boundary particle 
	    // and check if sub box of particle 1 is outside of this sub box
	    vec p1_sub_ptcl1_local2, p2_sub_ptcl1_local2, p3_sub_ptcl1_local2, p4_sub_ptcl1_local2;
	    vec p5_sub_ptcl1_local2, p6_sub_ptcl1_local2, p7_sub_ptcl1_local2, p8_sub_ptcl1_local2;
	    REAL aplus_sub_ptcl2, aminus_sub_ptcl2;
	    REAL bplus_sub_ptcl2, bminus_sub_ptcl2;
	    REAL cplus_sub_ptcl2, cminus_sub_ptcl2;


	    switch (num_ptcl1){
		case 1:
	            // get local2 coordinates of the eight corner
		    // points of the sub box of particle 1
		    p1_sub_ptcl1_local2 = p1_ptcl1_local2;
		    p2_sub_ptcl1_local2 = p9_ptcl1_local2;
		    p3_sub_ptcl1_local2 = p21_ptcl1_local2;
		    p4_sub_ptcl1_local2 = p13_ptcl1_local2;
		    p5_sub_ptcl1_local2 = p17_ptcl1_local2;
		    p6_sub_ptcl1_local2 = p22_ptcl1_local2;
		    p7_sub_ptcl1_local2 = p27_ptcl1_local2;
		    p8_sub_ptcl1_local2 = p25_ptcl1_local2;

	    	    // get local1 coordinates of the boundaries 
		    // of the sub box of particle 1
		    aplus_sub_ptcl1 = aplus_ptcl1; aminus_sub_ptcl1 = 0;
		    bplus_sub_ptcl1 = bplus_ptcl1; bminus_sub_ptcl1 = 0;
		    cplus_sub_ptcl1 = cplus_ptcl1; cminus_sub_ptcl1 = 0; 

		    xsign_p1 = 1;
		    ysign_p1 = 1;
		    zsign_p1 = 1;

		    break;
		case 2:
	            // get local2 coordinates of the eight corner
		    // points of the sub box of particle 1
		    p1_sub_ptcl1_local2 = p9_ptcl1_local2;
		    p2_sub_ptcl1_local2 = p2_ptcl1_local2;
		    p3_sub_ptcl1_local2 = p14_ptcl1_local2;
		    p4_sub_ptcl1_local2 = p21_ptcl1_local2;
		    p5_sub_ptcl1_local2 = p22_ptcl1_local2;
		    p6_sub_ptcl1_local2 = p18_ptcl1_local2;
		    p7_sub_ptcl1_local2 = p23_ptcl1_local2;
		    p8_sub_ptcl1_local2 = p27_ptcl1_local2;

	    	    // get local1 coordinates of the boundaries 
		    // of the sub box of particle 1
		    aplus_sub_ptcl1 = 0; aminus_sub_ptcl1 = aminus_ptcl1;
		    bplus_sub_ptcl1 = bplus_ptcl1; bminus_sub_ptcl1 = 0;
		    cplus_sub_ptcl1 = cplus_ptcl1; cminus_sub_ptcl1 = 0; 

		    xsign_p1 = -1;
		    ysign_p1 = 1;
		    zsign_p1 = 1;

		    break;
		case 3:
	            // get local2 coordinates of the eight corner
		    // points of the sub box of particle 1
		    p1_sub_ptcl1_local2 = p21_ptcl1_local2;
		    p2_sub_ptcl1_local2 = p14_ptcl1_local2;
		    p3_sub_ptcl1_local2 = p3_ptcl1_local2;
		    p4_sub_ptcl1_local2 = p10_ptcl1_local2;
		    p5_sub_ptcl1_local2 = p27_ptcl1_local2;
		    p6_sub_ptcl1_local2 = p23_ptcl1_local2;
		    p7_sub_ptcl1_local2 = p19_ptcl1_local2;
		    p8_sub_ptcl1_local2 = p24_ptcl1_local2;

	    	    // get local1 coordinates of the boundaries 
		    // of the sub box of particle 1
		    aplus_sub_ptcl1 = 0; aminus_sub_ptcl1 = aminus_ptcl1;
		    bplus_sub_ptcl1 = 0; bminus_sub_ptcl1 = bminus_ptcl1;
		    cplus_sub_ptcl1 = cplus_ptcl1; cminus_sub_ptcl1 = 0; 

		    xsign_p1 = -1;
		    ysign_p1 = -1;
		    zsign_p1 = 1;

		    break;
		case 4:
	            // get local2 coordinates of the eight corner
		    // points of the sub box of particle 1
		    p1_sub_ptcl1_local2 = p13_ptcl1_local2;
		    p2_sub_ptcl1_local2 = p21_ptcl1_local2;
		    p3_sub_ptcl1_local2 = p10_ptcl1_local2;
		    p4_sub_ptcl1_local2 = p4_ptcl1_local2;
		    p5_sub_ptcl1_local2 = p25_ptcl1_local2;
		    p6_sub_ptcl1_local2 = p27_ptcl1_local2;
		    p7_sub_ptcl1_local2 = p24_ptcl1_local2;
		    p8_sub_ptcl1_local2 = p20_ptcl1_local2;

	    	    // get local1 coordinates of the boundaries 
		    // of the sub box of particle 1
		    aplus_sub_ptcl1 = aplus_ptcl1; aminus_sub_ptcl1 = 0;
		    bplus_sub_ptcl1 = 0; bminus_sub_ptcl1 = bminus_ptcl1;
		    cplus_sub_ptcl1 = cplus_ptcl1; cminus_sub_ptcl1 = 0; 

		    xsign_p1 = 1;
		    ysign_p1 = -1;
		    zsign_p1 = 1;

		    break;
		case 5:
	            // get local2 coordinates of the eight corner
		    // points of the sub box of particle 1
		    p1_sub_ptcl1_local2 = p17_ptcl1_local2;
		    p2_sub_ptcl1_local2 = p22_ptcl1_local2;
		    p3_sub_ptcl1_local2 = p27_ptcl1_local2;
		    p4_sub_ptcl1_local2 = p25_ptcl1_local2;
		    p5_sub_ptcl1_local2 = p5_ptcl1_local2;
		    p6_sub_ptcl1_local2 = p12_ptcl1_local2;
		    p7_sub_ptcl1_local2 = p26_ptcl1_local2;
		    p8_sub_ptcl1_local2 = p16_ptcl1_local2;

	    	    // get local1 coordinates of the boundaries 
		    // of the sub box of particle 1
		    aplus_sub_ptcl1 = aplus_ptcl1; aminus_sub_ptcl1 = 0;
		    bplus_sub_ptcl1 = bplus_ptcl1; bminus_sub_ptcl1 = 0;
		    cplus_sub_ptcl1 = 0; cminus_sub_ptcl1 = cminus_ptcl1; 

		    xsign_p1 = 1;
		    ysign_p1 = 1;
		    zsign_p1 = -1;

		    break;
		case 6:
	            // get local2 coordinates of the eight corner
		    // points of the sub box of particle 1
		    p1_sub_ptcl1_local2 = p22_ptcl1_local2;
		    p2_sub_ptcl1_local2 = p18_ptcl1_local2;
		    p3_sub_ptcl1_local2 = p23_ptcl1_local2;
		    p4_sub_ptcl1_local2 = p27_ptcl1_local2;
		    p5_sub_ptcl1_local2 = p12_ptcl1_local2;
		    p6_sub_ptcl1_local2 = p6_ptcl1_local2;
		    p7_sub_ptcl1_local2 = p15_ptcl1_local2;
		    p8_sub_ptcl1_local2 = p26_ptcl1_local2;

	    	    // get local1 coordinates of the boundaries 
		    // of the sub box of particle 1
		    aplus_sub_ptcl1 = 0; aminus_sub_ptcl1 = aminus_ptcl1;
		    bplus_sub_ptcl1 = bplus_ptcl1; bminus_sub_ptcl1 = 0;
		    cplus_sub_ptcl1 = 0; cminus_sub_ptcl1 = cminus_ptcl1; 	

		    xsign_p1 = -1;
		    ysign_p1 = 1;
		    zsign_p1 = -1; 

		    break;
		case 7:
	            // get local2 coordinates of the eight corner
		    // points of the sub box of particle 1
		    p1_sub_ptcl1_local2 = p27_ptcl1_local2;
		    p2_sub_ptcl1_local2 = p23_ptcl1_local2;
		    p3_sub_ptcl1_local2 = p19_ptcl1_local2;
		    p4_sub_ptcl1_local2 = p24_ptcl1_local2;
		    p5_sub_ptcl1_local2 = p26_ptcl1_local2;
		    p6_sub_ptcl1_local2 = p15_ptcl1_local2;
		    p7_sub_ptcl1_local2 = p7_ptcl1_local2;
		    p8_sub_ptcl1_local2 = p11_ptcl1_local2;

	    	    // get local1 coordinates of the boundaries 
		    // of the sub box of particle 1
		    aplus_sub_ptcl1 = 0; aminus_sub_ptcl1 = aminus_ptcl1;
		    bplus_sub_ptcl1 = 0; bminus_sub_ptcl1 = bminus_ptcl1;
		    cplus_sub_ptcl1 = 0; cminus_sub_ptcl1 = cminus_ptcl1; 	

		    xsign_p1 = -1;
		    ysign_p1 = -1;
		    zsign_p1 = -1; 

		    break;
		case 8:
	            // get local2 coordinates of the eight corner
		    // points of the sub box of particle 1
		    p1_sub_ptcl1_local2 = p25_ptcl1_local2;
		    p2_sub_ptcl1_local2 = p27_ptcl1_local2;
		    p3_sub_ptcl1_local2 = p24_ptcl1_local2;
		    p4_sub_ptcl1_local2 = p20_ptcl1_local2;
		    p5_sub_ptcl1_local2 = p16_ptcl1_local2;
		    p6_sub_ptcl1_local2 = p26_ptcl1_local2;
		    p7_sub_ptcl1_local2 = p11_ptcl1_local2;
		    p8_sub_ptcl1_local2 = p8_ptcl1_local2;

	    	    // get local1 coordinates of the boundaries 
		    // of the sub box of particle 1
		    aplus_sub_ptcl1 = aplus_ptcl1; aminus_sub_ptcl1 = 0;
		    bplus_sub_ptcl1 = 0; bminus_sub_ptcl1 = bminus_ptcl1;
		    cplus_sub_ptcl1 = 0; cminus_sub_ptcl1 = cminus_ptcl1;  	 

		    xsign_p1 = 1;
		    ysign_p1 = -1;
		    zsign_p1 = -1;

		    break;

	    } // switch(num_ptcl1)



	    switch (num_ptcl2){
		case 1:
	            // get local1 coordinates of the eight corner
		    // points of the sub box of particle 2
		    p1_sub_ptcl2_local1 = p1_ptcl2_local1;
		    p2_sub_ptcl2_local1 = p9_ptcl2_local1;
		    p3_sub_ptcl2_local1 = p21_ptcl2_local1;
		    p4_sub_ptcl2_local1 = p13_ptcl2_local1;
		    p5_sub_ptcl2_local1 = p17_ptcl2_local1;
		    p6_sub_ptcl2_local1 = p22_ptcl2_local1;
		    p7_sub_ptcl2_local1 = p27_ptcl2_local1;
		    p8_sub_ptcl2_local1 = p25_ptcl2_local1;

	    	    // get local2 coordinates of the boundaries 
		    // of the sub box of particle 2
		    aplus_sub_ptcl2 = aplus_ptcl2; aminus_sub_ptcl2 = 0;
		    bplus_sub_ptcl2 = bplus_ptcl2; bminus_sub_ptcl2 = 0;
		    cplus_sub_ptcl2 = cplus_ptcl2; cminus_sub_ptcl2 = 0; 

		    xsign_p2 = 1;
		    ysign_p2 = 1;
		    zsign_p2 = 1;

		    break;
		case 2:
	            // get local1 coordinates of the eight corner
		    // points of the sub box of particle 2
		    p1_sub_ptcl2_local1 = p9_ptcl2_local1;
		    p2_sub_ptcl2_local1 = p2_ptcl2_local1;
		    p3_sub_ptcl2_local1 = p14_ptcl2_local1;
		    p4_sub_ptcl2_local1 = p21_ptcl2_local1;
		    p5_sub_ptcl2_local1 = p22_ptcl2_local1;
		    p6_sub_ptcl2_local1 = p18_ptcl2_local1;
		    p7_sub_ptcl2_local1 = p23_ptcl2_local1;
		    p8_sub_ptcl2_local1 = p27_ptcl2_local1;

	    	    // get local2 coordinates of the boundaries 
		    // of the sub box of particle 2
		    aplus_sub_ptcl2 = 0; aminus_sub_ptcl2 = aminus_ptcl2;
		    bplus_sub_ptcl2 = bplus_ptcl2; bminus_sub_ptcl2 = 0;
		    cplus_sub_ptcl2 = cplus_ptcl2; cminus_sub_ptcl2 = 0; 

		    xsign_p2 = -1;
		    ysign_p2 = 1;
		    zsign_p2 = 1;

		    break;
		case 3:
	            // get local1 coordinates of the eight corner
		    // points of the sub box of particle 2
		    p1_sub_ptcl2_local1 = p21_ptcl2_local1;
		    p2_sub_ptcl2_local1 = p14_ptcl2_local1;
		    p3_sub_ptcl2_local1 = p3_ptcl2_local1;
		    p4_sub_ptcl2_local1 = p10_ptcl2_local1;
		    p5_sub_ptcl2_local1 = p27_ptcl2_local1;
		    p6_sub_ptcl2_local1 = p23_ptcl2_local1;
		    p7_sub_ptcl2_local1 = p19_ptcl2_local1;
		    p8_sub_ptcl2_local1 = p24_ptcl2_local1;

	    	    // get local2 coordinates of the boundaries 
		    // of the sub box of particle 2
		    aplus_sub_ptcl2 = 0; aminus_sub_ptcl2 = aminus_ptcl2;
		    bplus_sub_ptcl2 = 0; bminus_sub_ptcl2 = bminus_ptcl2;
		    cplus_sub_ptcl2 = cplus_ptcl2; cminus_sub_ptcl2 = 0; 

		    xsign_p2 = -1;
		    ysign_p2 = -1;
		    zsign_p2 = 1;

		    break;
		case 4:
	            // get local1 coordinates of the eight corner
		    // points of the sub box of particle 2
		    p1_sub_ptcl2_local1 = p13_ptcl2_local1;
		    p2_sub_ptcl2_local1 = p21_ptcl2_local1;
		    p3_sub_ptcl2_local1 = p10_ptcl2_local1;
		    p4_sub_ptcl2_local1 = p4_ptcl2_local1;
		    p5_sub_ptcl2_local1 = p25_ptcl2_local1;
		    p6_sub_ptcl2_local1 = p27_ptcl2_local1;
		    p7_sub_ptcl2_local1 = p24_ptcl2_local1;
		    p8_sub_ptcl2_local1 = p20_ptcl2_local1;

	    	    // get local2 coordinates of the boundaries 
		    // of the sub box of particle 2
		    aplus_sub_ptcl2 = aplus_ptcl2; aminus_sub_ptcl2 = 0;
		    bplus_sub_ptcl2 = 0; bminus_sub_ptcl2 = bminus_ptcl2;
		    cplus_sub_ptcl2 = cplus_ptcl2; cminus_sub_ptcl2 = 0; 

		    xsign_p2 = 1;
		    ysign_p2 = -1;
		    zsign_p2 = 1;

		    break;
		case 5:
	            // get local1 coordinates of the eight corner
		    // points of the sub box of particle 2
		    p1_sub_ptcl2_local1 = p17_ptcl2_local1;
		    p2_sub_ptcl2_local1 = p22_ptcl2_local1;
		    p3_sub_ptcl2_local1 = p27_ptcl2_local1;
		    p4_sub_ptcl2_local1 = p25_ptcl2_local1;
		    p5_sub_ptcl2_local1 = p5_ptcl2_local1;
		    p6_sub_ptcl2_local1 = p12_ptcl2_local1;
		    p7_sub_ptcl2_local1 = p26_ptcl2_local1;
		    p8_sub_ptcl2_local1 = p16_ptcl2_local1;

	    	    // get local2 coordinates of the boundaries 
		    // of the sub box of particle 2
		    aplus_sub_ptcl2 = aplus_ptcl2; aminus_sub_ptcl2 = 0;
		    bplus_sub_ptcl2 = bplus_ptcl2; bminus_sub_ptcl2 = 0;
		    cplus_sub_ptcl2 = 0; cminus_sub_ptcl2 = cminus_ptcl2; 

		    xsign_p2 = 1;
		    ysign_p2 = 1;
		    zsign_p2 = -1;

		    break;
		case 6:
	            // get local1 coordinates of the eight corner
		    // points of the sub box of particle 2
		    p1_sub_ptcl2_local1 = p22_ptcl2_local1;
		    p2_sub_ptcl2_local1 = p18_ptcl2_local1;
		    p3_sub_ptcl2_local1 = p23_ptcl2_local1;
		    p4_sub_ptcl2_local1 = p27_ptcl2_local1;
		    p5_sub_ptcl2_local1 = p12_ptcl2_local1;
		    p6_sub_ptcl2_local1 = p6_ptcl2_local1;
		    p7_sub_ptcl2_local1 = p15_ptcl2_local1;
		    p8_sub_ptcl2_local1 = p26_ptcl2_local1;

	    	    // get local2 coordinates of the boundaries 
		    // of the sub box of particle 2
		    aplus_sub_ptcl2 = 0; aminus_sub_ptcl2 = aminus_ptcl2;
		    bplus_sub_ptcl2 = bplus_ptcl2; bminus_sub_ptcl2 = 0;
		    cplus_sub_ptcl2 = 0; cminus_sub_ptcl2 = cminus_ptcl2; 	

		    xsign_p2 = -1;
		    ysign_p2 = 1;
		    zsign_p2 = -1; 

		    break;
		case 7:
	            // get local1 coordinates of the eight corner
		    // points of the sub box of particle 2
		    p1_sub_ptcl2_local1 = p27_ptcl2_local1;
		    p2_sub_ptcl2_local1 = p23_ptcl2_local1;
		    p3_sub_ptcl2_local1 = p19_ptcl2_local1;
		    p4_sub_ptcl2_local1 = p24_ptcl2_local1;
		    p5_sub_ptcl2_local1 = p26_ptcl2_local1;
		    p6_sub_ptcl2_local1 = p15_ptcl2_local1;
		    p7_sub_ptcl2_local1 = p7_ptcl2_local1;
		    p8_sub_ptcl2_local1 = p11_ptcl2_local1;

	    	    // get local2 coordinates of the boundaries 
		    // of the sub box of particle 2
		    aplus_sub_ptcl2 = 0; aminus_sub_ptcl2 = aminus_ptcl2;
		    bplus_sub_ptcl2 = 0; bminus_sub_ptcl2 = bminus_ptcl2;
		    cplus_sub_ptcl2 = 0; cminus_sub_ptcl2 = cminus_ptcl2; 	

		    xsign_p2 = -1;
		    ysign_p2 = -1;
		    zsign_p2 = -1; 

		    break;
		case 8:
	            // get local1 coordinates of the eight corner
		    // points of the sub box of particle 2
		    p1_sub_ptcl2_local1 = p25_ptcl2_local1;
		    p2_sub_ptcl2_local1 = p27_ptcl2_local1;
		    p3_sub_ptcl2_local1 = p24_ptcl2_local1;
		    p4_sub_ptcl2_local1 = p20_ptcl2_local1;
		    p5_sub_ptcl2_local1 = p16_ptcl2_local1;
		    p6_sub_ptcl2_local1 = p26_ptcl2_local1;
		    p7_sub_ptcl2_local1 = p11_ptcl2_local1;
		    p8_sub_ptcl2_local1 = p8_ptcl2_local1;

	    	    // get local2 coordinates of the boundaries 
		    // of the sub box of particle 2
		    aplus_sub_ptcl2 = aplus_ptcl2; aminus_sub_ptcl2 = 0;
		    bplus_sub_ptcl2 = 0; bminus_sub_ptcl2 = bminus_ptcl2;
		    cplus_sub_ptcl2 = 0; cminus_sub_ptcl2 = cminus_ptcl2;  	 

		    xsign_p2 = 1;
		    ysign_p2 = -1;
		    zsign_p2 = -1;

		    break;

	    } // switch(num_ptcl2)



	    //////////////////////////////////////////////////////////////////////////////////////////////////////
    	    // check if all 8 corner points of the sub box of particle 2 are outside of the sub box of particle 1
    	    // x+ boundary, i.e. a+ boundary
    	    if(p1_sub_ptcl2_local1.getx()>aplus_sub_ptcl1 && p2_sub_ptcl2_local1.getx()>aplus_sub_ptcl1
    	    && p3_sub_ptcl2_local1.getx()>aplus_sub_ptcl1 && p4_sub_ptcl2_local1.getx()>aplus_sub_ptcl1
    	    && p5_sub_ptcl2_local1.getx()>aplus_sub_ptcl1 && p6_sub_ptcl2_local1.getx()>aplus_sub_ptcl1
    	    && p7_sub_ptcl2_local1.getx()>aplus_sub_ptcl1 && p8_sub_ptcl2_local1.getx()>aplus_sub_ptcl1 )	// outside of x+ boundary
		continue;	// sub boxes are not contact, go to the next sub box pair

    	    // x- boundary, i.e. a- boundary
    	    if(p1_sub_ptcl2_local1.getx()<-aminus_sub_ptcl1 && p2_sub_ptcl2_local1.getx()<-aminus_sub_ptcl1
    	    && p3_sub_ptcl2_local1.getx()<-aminus_sub_ptcl1 && p4_sub_ptcl2_local1.getx()<-aminus_sub_ptcl1
    	    && p5_sub_ptcl2_local1.getx()<-aminus_sub_ptcl1 && p6_sub_ptcl2_local1.getx()<-aminus_sub_ptcl1
    	    && p7_sub_ptcl2_local1.getx()<-aminus_sub_ptcl1 && p8_sub_ptcl2_local1.getx()<-aminus_sub_ptcl1 )	// outside of x- boundary
		continue;	// sub boxes are not contact, go to the next sub box pair

    	    // y+ boundary, i.e. b+ boundary
    	    if(p1_sub_ptcl2_local1.gety()>bplus_sub_ptcl1 && p2_sub_ptcl2_local1.gety()>bplus_sub_ptcl1
    	    && p3_sub_ptcl2_local1.gety()>bplus_sub_ptcl1 && p4_sub_ptcl2_local1.gety()>bplus_sub_ptcl1
    	    && p5_sub_ptcl2_local1.gety()>bplus_sub_ptcl1 && p6_sub_ptcl2_local1.gety()>bplus_sub_ptcl1
    	    && p7_sub_ptcl2_local1.gety()>bplus_sub_ptcl1 && p8_sub_ptcl2_local1.gety()>bplus_sub_ptcl1 )	// outside of y+ boundary
		continue;	// sub boxes are not contact, go to the next sub box pair

    	    // y- boundary, i.e. b- boundary
    	    if(p1_sub_ptcl2_local1.gety()<-bminus_sub_ptcl1 && p2_sub_ptcl2_local1.gety()<-bminus_sub_ptcl1
    	    && p3_sub_ptcl2_local1.gety()<-bminus_sub_ptcl1 && p4_sub_ptcl2_local1.gety()<-bminus_sub_ptcl1
    	    && p5_sub_ptcl2_local1.gety()<-bminus_sub_ptcl1 && p6_sub_ptcl2_local1.gety()<-bminus_sub_ptcl1
    	    && p7_sub_ptcl2_local1.gety()<-bminus_sub_ptcl1 && p8_sub_ptcl2_local1.gety()<-bminus_sub_ptcl1 )	// outside of y- boundary
		continue;	// sub boxes are not contact, go to the next sub box pair

    	    // z+ boundary, i.e. c+ boundary
    	    if(p1_sub_ptcl2_local1.getz()>cplus_sub_ptcl1 && p2_sub_ptcl2_local1.getz()>cplus_sub_ptcl1
    	    && p3_sub_ptcl2_local1.getz()>cplus_sub_ptcl1 && p4_sub_ptcl2_local1.getz()>cplus_sub_ptcl1
    	    && p5_sub_ptcl2_local1.getz()>cplus_sub_ptcl1 && p6_sub_ptcl2_local1.getz()>cplus_sub_ptcl1
    	    && p7_sub_ptcl2_local1.getz()>cplus_sub_ptcl1 && p8_sub_ptcl2_local1.getz()>cplus_sub_ptcl1 )	// outside of z+ boundary
		continue;	// sub boxes are not contact, go to the next sub box pair

    	    // z- boundary, i.e. c- boundary
    	    if(p1_sub_ptcl2_local1.getz()<-cminus_sub_ptcl1 && p2_sub_ptcl2_local1.getz()<-cminus_sub_ptcl1
    	    && p3_sub_ptcl2_local1.getz()<-cminus_sub_ptcl1 && p4_sub_ptcl2_local1.getz()<-cminus_sub_ptcl1
    	    && p5_sub_ptcl2_local1.getz()<-cminus_sub_ptcl1 && p6_sub_ptcl2_local1.getz()<-cminus_sub_ptcl1
    	    && p7_sub_ptcl2_local1.getz()<-cminus_sub_ptcl1 && p8_sub_ptcl2_local1.getz()<-cminus_sub_ptcl1 )	// outside of z- boundary
		continue;	// sub boxes are not contact, go to the next sub box pair


	    //////////////////////////////////////////////////////////////////////////////////////////////////////
    	    // check if all 8 corner points of the sub box of particle 1 are outside of the sub box of particle 2
    	    // x+ boundary, i.e. a+ boundary
    	    if(p1_sub_ptcl1_local2.getx()>aplus_sub_ptcl2 && p2_sub_ptcl1_local2.getx()>aplus_sub_ptcl2
    	    && p3_sub_ptcl1_local2.getx()>aplus_sub_ptcl2 && p4_sub_ptcl1_local2.getx()>aplus_sub_ptcl2
    	    && p5_sub_ptcl1_local2.getx()>aplus_sub_ptcl2 && p6_sub_ptcl1_local2.getx()>aplus_sub_ptcl2
    	    && p7_sub_ptcl1_local2.getx()>aplus_sub_ptcl2 && p8_sub_ptcl1_local2.getx()>aplus_sub_ptcl2 )	// outside of x+ boundary
		continue;	// sub boxes are not contact, go to the next sub box pair

    	    // x- boundary, i.e. a- boundary
    	    if(p1_sub_ptcl1_local2.getx()<-aminus_sub_ptcl2 && p2_sub_ptcl1_local2.getx()<-aminus_sub_ptcl2
    	    && p3_sub_ptcl1_local2.getx()<-aminus_sub_ptcl2 && p4_sub_ptcl1_local2.getx()<-aminus_sub_ptcl2
    	    && p5_sub_ptcl1_local2.getx()<-aminus_sub_ptcl2 && p6_sub_ptcl1_local2.getx()<-aminus_sub_ptcl2
    	    && p7_sub_ptcl1_local2.getx()<-aminus_sub_ptcl2 && p8_sub_ptcl1_local2.getx()<-aminus_sub_ptcl2 )	// outside of x- boundary
		continue;	// sub boxes are not contact, go to the next sub box pair

    	    // y+ boundary, i.e. b+ boundary
    	    if(p1_sub_ptcl1_local2.gety()>bplus_sub_ptcl2 && p2_sub_ptcl1_local2.gety()>bplus_sub_ptcl2
    	    && p3_sub_ptcl1_local2.gety()>bplus_sub_ptcl2 && p4_sub_ptcl1_local2.gety()>bplus_sub_ptcl2
    	    && p5_sub_ptcl1_local2.gety()>bplus_sub_ptcl2 && p6_sub_ptcl1_local2.gety()>bplus_sub_ptcl2
    	    && p7_sub_ptcl1_local2.gety()>bplus_sub_ptcl2 && p8_sub_ptcl1_local2.gety()>bplus_sub_ptcl2 )	// outside of y+ boundary
		continue;	// sub boxes are not contact, go to the next sub box pair

    	    // y- boundary, i.e. b- boundary
    	    if(p1_sub_ptcl1_local2.gety()<-bminus_sub_ptcl2 && p2_sub_ptcl1_local2.gety()<-bminus_sub_ptcl2
    	    && p3_sub_ptcl1_local2.gety()<-bminus_sub_ptcl2 && p4_sub_ptcl1_local2.gety()<-bminus_sub_ptcl2
    	    && p5_sub_ptcl1_local2.gety()<-bminus_sub_ptcl2 && p6_sub_ptcl1_local2.gety()<-bminus_sub_ptcl2
    	    && p7_sub_ptcl1_local2.gety()<-bminus_sub_ptcl2 && p8_sub_ptcl1_local2.gety()<-bminus_sub_ptcl2 )	// outside of y- boundary
		continue;	// sub boxes are not contact, go to the next sub box pair

    	    // z+ boundary, i.e. c+ boundary
    	    if(p1_sub_ptcl1_local2.getz()>cplus_sub_ptcl2 && p2_sub_ptcl1_local2.getz()>cplus_sub_ptcl2
    	    && p3_sub_ptcl1_local2.getz()>cplus_sub_ptcl2 && p4_sub_ptcl1_local2.getz()>cplus_sub_ptcl2
    	    && p5_sub_ptcl1_local2.getz()>cplus_sub_ptcl2 && p6_sub_ptcl1_local2.getz()>cplus_sub_ptcl2
    	    && p7_sub_ptcl1_local2.getz()>cplus_sub_ptcl2 && p8_sub_ptcl1_local2.getz()>cplus_sub_ptcl2 )	// outside of z+ boundary
		continue;	// sub boxes are not contact, go to the next sub box pair

    	    // z- boundary, i.e. c- boundary
    	    if(p1_sub_ptcl1_local2.getz()<-cminus_sub_ptcl2 && p2_sub_ptcl1_local2.getz()<-cminus_sub_ptcl2
    	    && p3_sub_ptcl1_local2.getz()<-cminus_sub_ptcl2 && p4_sub_ptcl1_local2.getz()<-cminus_sub_ptcl2
    	    && p5_sub_ptcl1_local2.getz()<-cminus_sub_ptcl2 && p6_sub_ptcl1_local2.getz()<-cminus_sub_ptcl2
    	    && p7_sub_ptcl1_local2.getz()<-cminus_sub_ptcl2 && p8_sub_ptcl1_local2.getz()<-cminus_sub_ptcl2 )	// outside of z- boundary
		continue;	// sub boxes are not contact, go to the next sub box pair



	    // USE THE CORRESPONDING OCTANTS TO CALCULATE THE CONTACTS OF ELLIPSOIDS

    	    p1->getGlobCoef(coef_p1, num_ptcl1); // v[0] is the point on p2, v[1] is the point on p1
    	    p2->getGlobCoef(coef_p2, num_ptcl2);

    	    b1 = root6(coef_p1,coef_p2,v[0]);
    	    b2 = root6(coef_p2,coef_p1,v[1]);
    	    point1 = v[1];
    	    point2 = v[0];

    	    local1_point1 = p1->localVec(point1 - p1->getCurrPosition());
    	    local2_point2 = p2->localVec(point2 - p2->getCurrPosition());

    	    if(b1 && b2 
    	    && local1_point1.getx()*xsign_p1 >= 0 && local1_point1.gety()*ysign_p1 >= 0 && local1_point1.getz()*zsign_p1 >= 0
    	    && local2_point2.getx()*xsign_p2 >= 0 && local2_point2.gety()*ysign_p2 >= 0 && local2_point2.getz()*zsign_p2 >= 0 ){
		flag_contact = true;
		break;	// out of loop num_ptcl2
	    }

	} // end num_ptcl2

	if(flag_contact == true) // has already found contacts
	    break;

    } // end num_ptcl1

    if(flag_contact == false) 	// out of loop without finding the contacts
	return false;




//////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////       end of new method, as in notes page 106, April 2, 2014       //////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////









    radius1=p1->getRadius(point1, xsign_p1, ysign_p1, zsign_p1);
    radius2=p2->getRadius(point2, xsign_p2, ysign_p2, zsign_p2);
//std::cout << "radius with octant in isOverlapped():    " << radius1 << "  " << radius2 << std::endl;
//    radius1 = p1->getRadius(point1);
//    radius2 = p2->getRadius(point2);
//std::cout << "radius without octant in isOverlapped(): " << radius1 << "  " << radius2 << std::endl;
    penetration=vfabs(point1-point2);
    REAL minRelOverlap = penetration/(2.0*fmax(radius1,radius2));

    if (b1 && b2 
	&& minRelOverlap > MINOVERLAP
	&& nearbyint(penetration/MEPS >= 1) ) {
        isInContact = true;
/*
std::cout << std::endl;
std::cout << point1.getx() << "   " << point1.gety() << "   " << point1.getz() << " ; " << point2.getx() << "   " << point2.gety() << "   " << point2.getz() << " ;" << std::endl << std::setw(16) << p1->getAplus() << std::setw(16) << p1->getAminus() << std::setw(16) << p1->getBplus() << std::setw(16) << p1->getBminus() << std::setw(16) << p1->getCplus() << std::setw(16) << p1->getCminus() << std::setw(16) << p1->getCurrPosition().getx() << std::setw(16) << p1->getCurrPosition().gety() << std::setw(16) << p1->getCurrPosition().getz() << std::setw(16) << p1->getCurrDirecA().getx() << std::setw(16) << p1->getCurrDirecA().gety() << std::setw(16) << p1->getCurrDirecA().getz() << std::setw(16) << p1->getCurrDirecB().getx() << std::setw(16) << p1->getCurrDirecB().gety() << std::setw(16) << p1->getCurrDirecB().getz() << std::setw(16) << p1->getCurrDirecC().getx() << std::setw(16) << p1->getCurrDirecC().gety() << std::setw(16) << p1->getCurrDirecC().getz() << std::endl << std::setw(16) << p2->getAplus() << std::setw(16) << p2->getAminus() << std::setw(16) << p2->getBplus() << std::setw(16) << p2->getBminus() << std::setw(16) << p2->getCplus() << std::setw(16) << p2->getCminus() << std::setw(16) << p2->getCurrPosition().getx() << std::setw(16) << p2->getCurrPosition().gety() << std::setw(16) << p2->getCurrPosition().getz() << std::setw(16) << p2->getCurrDirecA().getx() << std::setw(16) << p2->getCurrDirecA().gety() << std::setw(16) << p2->getCurrDirecA().getz() << std::setw(16) << p2->getCurrDirecB().getx() << std::setw(16) << p2->getCurrDirecB().gety() << std::setw(16) << p2->getCurrDirecB().getz() << std::setw(16) << p2->getCurrDirecC().getx() << std::setw(16) << p2->getCurrDirecC().gety() << std::setw(16) << p2->getCurrDirecC().getz() << std::endl;
*/
        return true;
    }
    else {
        isInContact = false;
	return false;
    }
}

/*
template<class T>
void contact<T>::Hopkin_contact(vec& oct1, vec& oct2){	// August 27, 2013
	REAL R_ratio = 0.5;	// the ratio of R to the smallest among aplus, aminus, bplus, ...
	REAL d_alpha = 0.05;	// P+alpha*d

	vec local1_n1, local2_n2;	// local coordinates of unit normal
	vec local1_dt, local2_dt;
	vec local1_P1prime, local2_P2prime;
	REAL P1_gamma, P2_gamma;
//	vec prev_P1, prev_P2;	// previous local P1 and P2
	REAL prev_dist_d;

//	REAL dilate_R1, dilate_R2;
//	dilate_R1 = -p1->getMinRadius()*R_ratio;	// a is the longest in ellip3d
//	dilate_R2 = -p2->getMinRadius()*R_ratio;	// to avoid c-dilate_R to be negative

	REAL a1, b1, c1;	// parameters of particle 1, P start from (aplus, bplus, cplus)
	a1 = p1->getAplus();
	b1 = p1->getBplus();
	c1 = p1->getCplus();

	REAL a2, b2, c2;	// parameters of particle 2
	a2 = p2->getAplus();
	b2 = p2->getBplus();
	c2 = p2->getCplus();

	vec center1_geo = p1->getCurrPosition();
	vec center2_geo = p2->getCurrPosition();

	// local coordinates of starting points are the points in 1st octant 
	vec local1_P1 = vec(a1, 0, 0);
	vec local2_P2 = vec(a2, 0, 0);

	vec global_d = center2_geo+p2->globalVec(local2_P2)-center1_geo-p1->globalVec(local1_P1);	// from particle 1 to 2
	REAL dist_d = vfabs(global_d);	// distance of vector d
	prev_dist_d = dist_d+2.0;	

	while(dist_d < prev_dist_d){
//		prev_P1 = local1_P1; 
//		prev_P2 = local2_P2;
		prev_dist_d = dist_d;
		// update p in particle 1

		// local coordinates of unit normal
		local1_n1 = vec(local1_P1.getx()/a1/a1, local1_P1.gety()/b1/b1, local1_P1.getz()/c1/c1);
		local1_n1 = normalize(local1_n1);
		// local tangential vector dt
		local1_dt = p1->localVec(global_d)-(p1->localVec(global_d)*local1_n1)*local1_n1;
		
		local1_P1prime = local1_P1+d_alpha*local1_dt;
		// update a1, b1, c1 
		if(local1_P1prime.getx()>0)
			a1 = p1->getAplus();
		else
			a1 = p1->getAminus();
		if(local1_P1prime.gety()>0)
			b1 = p1->getBplus();
		else
			b1 = p1->getBminus();
		if(local1_P1prime.getz()>0)
			c1 = p1->getCplus();
		else
			c1 = p1->getCminus();

		P1_gamma = 1.0/sqrt(pow(local1_P1prime.getx()/a1,2)+pow(local1_P1prime.gety()/b1,2)+pow(local1_P1prime.getz()/c1,2));
		// update local1_P1 for next iteration
		local1_P1 = P1_gamma*local1_P1prime;

		// update global d, Octorber 6, 2013
		global_d = center2_geo+p2->globalVec(local2_P2)-center1_geo-p1->globalVec(local1_P1);

		// update p in particle 2

		// local coordinates of unit normal
		local2_n2 = vec(local2_P2.getx()/a2/a2, local2_P2.gety()/b2/b2, local2_P2.getz()/c2/c2);
		local2_n2 = normalize(local2_n2);
		// local tangential vector dt
		local2_dt = p2->localVec(-global_d)-(p2->localVec(-global_d)*local2_n2)*local2_n2;
		
		local2_P2prime = local2_P2+d_alpha*local2_dt;
		// update a2, b2, c2 
		if(local2_P2prime.getx()>0)
			a2 = p2->getAplus();
		else
			a2 = p2->getAminus();
		if(local2_P2prime.gety()>0)
			b2 = p2->getBplus();
		else
			b2 = p2->getBminus();
		if(local2_P2prime.getz()>0)
			c2 = p2->getCplus();
		else
			c2 = p2->getCminus();

		P2_gamma = 1.0/sqrt(pow(local2_P2prime.getx()/a2,2)+pow(local2_P2prime.gety()/b2,2)/+pow(local2_P2prime.getz()/c2,2));
		// update local2_P2 for next iteration
		local2_P2 = P2_gamma*local2_P2prime;

		// distance of new vector d
		global_d = center2_geo+p2->globalVec(local2_P2)-center1_geo-p1->globalVec(local1_P1);
		dist_d = vfabs(global_d);	
	
	}

	// find the octants
	// particle 1
	if(prev_P1.getx()>0)
		oct1.setx(1);
	else
		oct1.setx(-1);
	if(prev_P1.gety()>0)
		oct1.sety(1);
	else
		oct1.sety(-1);
	if(prev_P1.getz()>0)
		oct1.setz(1);
	else
		oct1.setz(-1);
	// particle 2
	if(prev_P2.getx()>0)
		oct2.setx(1);
	else
		oct2.setx(-1);
	if(prev_P2.gety()>0)
		oct2.sety(1);
	else
		oct2.sety(-1);
	if(prev_P2.getz()>0)
		oct2.setz(1);
	else
		oct2.setz(-1);


	// find the octants using current vector
	// particle 1
	if(local1_P1.getx()>0)
		oct1.setx(1);
	else
		oct1.setx(-1);
	if(local1_P1.gety()>0)
		oct1.sety(1);
	else
		oct1.sety(-1);
	if(local1_P1.getz()>0)
		oct1.setz(1);
	else
		oct1.setz(-1);
	// particle 2
	if(local2_P2.getx()>0)
		oct2.setx(1);
	else
		oct2.setx(-1);
	if(local2_P2.gety()>0)
		oct2.sety(1);
	else
		oct2.sety(-1);
	if(local2_P2.getz()>0)
		oct2.setz(1);
	else
		oct2.setz(-1);

}
*/



template<class T>
void contact<T>::Hopkin_contact(vec& oct1, vec& oct2){	// August 26, 2013
	REAL R_ratio = 0.7;	// the ratio of R to the smallest among aplus, aminus, bplus, ...
	REAL d_alpha = 0.05;	// P+alpha*d

	vec local1_grad1, local2_grad2;	// local coordinates of gradient
	vec local1_n1, local2_n2;	// local coordinates of unit normal
	vec local1_dt, local2_dt;
	vec local1_P1prime, local2_P2prime;
	REAL P1_gamma, P2_gamma;
	vec prev_P1, prev_P2;	// previous local P1 and P2
//	REAL prev_dist_d;

//	REAL dilate_R1, dilate_R2;
//	dilate_R1 = p1->getMinRadius()*R_ratio;	// a is the longest in ellip3d
//	dilate_R2 = p2->getMinRadius()*R_ratio;	// to avoid c-dilate_R to be negative

	REAL a1, b1, c1;	// parameters of particle 1
	a1 = p1->getAplus();
	b1 = p1->getBplus();
	c1 = p1->getCplus();

	REAL a2, b2, c2;	// parameters of particle 2
	a2 = p2->getAplus();
	b2 = p2->getBplus();
	c2 = p2->getCplus();

	vec center1_geo = p1->getCurrPosition();
	vec center2_geo = p2->getCurrPosition();

	// local coordinates of starting points are the points in 1st octant 
	vec local1_P1 = vec(a1, 0, 0);
	vec local2_P2 = vec(a2, 0, 0);

	vec global_d = center2_geo+p2->globalVec(local2_P2)-center1_geo-p1->globalVec(local1_P1);	// from particle 1 to 2
//	REAL dist_d = vfabs(global_d);	// distance of vector d
//	prev_dist_d = dist_d+2.0;

	REAL dtol_1 = p1->getMaxRadius()*0.00001;
	REAL dtol_2 = p2->getMaxRadius()*0.00001;

	REAL change_P1 = p1->getMaxRadius();
	REAL change_P2 = p2->getMaxRadius();	

//	while(dist_d < prev_dist_d){
	while(change_P1 > dtol_1 || change_P2 > dtol_2){
		prev_P1 = local1_P1; 
		prev_P2 = local2_P2;
//		prev_dist_d = dist_d;
		// update p in particle 1

		// local coordinates of unit normal
		local1_grad1 = vec(local1_P1.getx()/a1/a1, local1_P1.gety()/b1/b1, local1_P1.getz()/c1/c1);
		local1_n1 = normalize(local1_grad1);
		// local tangential vector dt
		local1_dt = p1->localVec(global_d)-(p1->localVec(global_d)*local1_n1)*local1_n1;
		
		local1_P1prime = local1_P1+d_alpha*local1_dt;
		// update a1, b1, c1 
		if(local1_P1prime.getx()>0)
			a1 = p1->getAplus();
		else
			a1 = p1->getAminus();
		if(local1_P1prime.gety()>0)
			b1 = p1->getBplus();
		else
			b1 = p1->getBminus();
		if(local1_P1prime.getz()>0)
			c1 = p1->getCplus();
		else
			c1 = p1->getCminus();

		P1_gamma = 1.0/sqrt(pow(local1_P1prime.getx()/a1,2)+pow(local1_P1prime.gety()/b1,2)+pow(local1_P1prime.getz()/c1,2));
		// update local1_P1 for next iteration
		local1_P1 = P1_gamma*local1_P1prime;

		// update global d, Octorber 6, 2013
		global_d = center2_geo+p2->globalVec(local2_P2)-center1_geo-p1->globalVec(local1_P1);

		// update p in particle 2

		// local coordinates of unit normal
		local2_grad2 = vec(local2_P2.getx()/a2/a2, local2_P2.gety()/b2/b2, local2_P2.getz()/c2/c2);
		local2_n2 = normalize(local2_grad2);
		// local tangential vector dt
		local2_dt = p2->localVec(-global_d)-(p2->localVec(-global_d)*local2_n2)*local2_n2;
		
		local2_P2prime = local2_P2+d_alpha*local2_dt;
		// update a2, b2, c2 
		if(local2_P2prime.getx()>0)
			a2 = p2->getAplus();
		else
			a2 = p2->getAminus();
		if(local2_P2prime.gety()>0)
			b2 = p2->getBplus();
		else
			b2 = p2->getBminus();
		if(local2_P2prime.getz()>0)
			c2 = p2->getCplus();
		else
			c2 = p2->getCminus();

		P2_gamma = 1.0/sqrt(pow(local2_P2prime.getx()/a2,2)+pow(local2_P2prime.gety()/b2,2)+pow(local2_P2prime.getz()/c2,2));
		// update local2_P2 for next iteration
		local2_P2 = P2_gamma*local2_P2prime;

		// distance of new vector d
		global_d = center2_geo+p2->globalVec(local2_P2)-center1_geo-p1->globalVec(local1_P1);
//		dist_d = vfabs(global_d);	

		change_P1 = vfabs(local1_P1-prev_P1);
		change_P2 = vfabs(local2_P2-prev_P2);	
	}

	// find the octants
	// particle 1
	if(local1_P1.getx()>0)
		oct1.setx(1);
	else
		oct1.setx(-1);
	if(local1_P1.gety()>0)
		oct1.sety(1);
	else
		oct1.sety(-1);
	if(local1_P1.getz()>0)
		oct1.setz(1);
	else
		oct1.setz(-1);
	// particle 2
	if(local2_P2.getx()>0)
		oct2.setx(1);
	else
		oct2.setx(-1);
	if(local2_P2.gety()>0)
		oct2.sety(1);
	else
		oct2.sety(-1);
	if(local2_P2.getz()>0)
		oct2.setz(1);
	else
		oct2.setz(-1);


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
void contact<T>::contactForce(bool &exceed){	// exceed August, 19. 2013 apply moment on the mass center
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
	p1->addMoment( ( (point1+point2)/2-p1->getCurrCenterMass() ) *   NormalForce );
	p2->addMoment( ( (point1+point2)/2-p2->getCurrCenterMass() ) * (-NormalForce) );	
	
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
	vec veloc1 = p1->getCurrVelocity() + p1->getCurrOmga()*(cp-p1->getCurrCenterMass());
	vec veloc2 = p2->getCurrVelocity() + p2->getCurrOmga()*(cp-p2->getCurrCenterMass());
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
	p1->addMoment( ( (point1+point2)/2-p1->getCurrCenterMass() ) * (-CntDampingForce) );
	p2->addMoment( ( (point1+point2)/2-p2->getCurrCenterMass() ) * CntDampingForce );

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
	    p1->addMoment( ((point1+point2)/2-p1->getCurrCenterMass())*TgtForce);
	    p2->addMoment(-((point1+point2)/2-p2->getCurrCenterMass())*TgtForce);
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

