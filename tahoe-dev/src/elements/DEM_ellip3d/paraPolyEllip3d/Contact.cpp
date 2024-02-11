#include "Contact.h"
#include "Parameter.h"
#include "const.h"
#include "root6.h"
#include <iostream>
#include <fstream>
#include <iomanip>

//#define MINDLIN_ASSUMED
//#define MINDLIN_KNOWN

namespace dem {
  
  Contact::Contact() {
    p1 = NULL;
    p2 = NULL;
    penetr = 0;
    contactRadius = 0;
    point1 = point2 = 0;
    radius1 = radius2 = 0;
    normalDirc = tgtDirc = 0;
    
    isInContact = false;
    tgtLoading = prevTgtLoading = true;
    normalForce = prevNormalForce = 0;
    tgtForce = prevTgtForce = 0;
    tgtDisp = prevTgtDisp = 0;
    tgtDispStart = 0;
    tgtSlide = prevTgtSlide = false;
    tgtPeak = 0;

    cohesionForce = 0;
    spinResist = 0;

    E0 = G0 = R0 = 0;
  }
  
  
  Contact::Contact(Particle* t1, Particle* t2){
    p1 = t1;
    p2 = t2;
    penetr = 0;
    contactRadius = 0;
    point1 = point2 = 0;
    radius1 = radius2 = 0;
    normalDirc = tgtDirc = 0;

    isInContact = false;
    tgtLoading = prevTgtLoading = true;
    normalForce = prevNormalForce = 0;
    tgtForce = prevTgtForce = 0;
    tgtDisp = prevTgtDisp = 0;
    tgtDispStart = 0;
    tgtSlide = prevTgtSlide = false;
    tgtPeak = 0;

    cohesionForce = 0;
    spinResist = 0;

    E0 = G0 = R0 = 0;
  }
  

  bool Contact::isRedundant(const Contact &other) const {
    std::size_t id1 = getP1() -> getId();
    std::size_t id2 = getP2() -> getId();
    std::size_t oId1 = ( other.getP1() ) -> getId();
    std::size_t oId2 = ( other.getP2() ) -> getId();
    
    return ( (id2 == oId1 && id1 == oId2) || (id1 == oId1 && id2 == oId2) );
  }


  bool Contact::operator==(const Contact &other) const {
    std::size_t id1 = getP1() -> getId();
    std::size_t id2 = getP2() -> getId();
    std::size_t oId1 = ( other.getP1() ) -> getId();
    std::size_t oId2 = ( other.getP2() ) -> getId();
    
    return ( (id2 == oId1 && id1 == oId2) || (id1 == oId1 && id2 == oId2) );
  }

  
  Particle* Contact::getP1() const {
    return p1;
  }
  

  Particle* Contact::getP2() const {
    return p2;
  }
  

  bool Contact::isOverlapped(){

    /////////////////////////////////////////////////////////////////
    /////////////////// pre-detection step1  ////////////////////////
    /////////////////////////////////////////////////////////////////

    // if all corner points of one particle are outside of the other particle,
    // then this possible contact pair cannot be in contact, we delete this contact pair,
    // March 27, 2014. See case (1) in the notes in page 101

    // local coordinates of corner points of the two particles
    Vec p1_ptcl1_local1, p2_ptcl1_local1, p3_ptcl1_local1, p4_ptcl1_local1, p5_ptcl1_local1, p6_ptcl1_local1, p7_ptcl1_local1, p8_ptcl1_local1;
    Vec p1_ptcl2_local2, p2_ptcl2_local2, p3_ptcl2_local2, p4_ptcl2_local2, p5_ptcl2_local2, p6_ptcl2_local2, p7_ptcl2_local2, p8_ptcl2_local2;

    // get local coordinates of particle 1
    REAL aplus_ptcl1 = p1->getAplus(); REAL aminus_ptcl1 = p1->getAminus(); 
    REAL bplus_ptcl1 = p1->getBplus(); REAL bminus_ptcl1 = p1->getBminus();  
    REAL cplus_ptcl1 = p1->getCplus(); REAL cminus_ptcl1 = p1->getCminus();

    p1_ptcl1_local1 = Vec(aplus_ptcl1, bplus_ptcl1, cplus_ptcl1);
    p2_ptcl1_local1 = Vec(-aminus_ptcl1, bplus_ptcl1, cplus_ptcl1);
    p3_ptcl1_local1 = Vec(-aminus_ptcl1, -bminus_ptcl1, cplus_ptcl1);
    p4_ptcl1_local1 = Vec(aplus_ptcl1, -bminus_ptcl1, cplus_ptcl1);
    p5_ptcl1_local1 = Vec(aplus_ptcl1, bplus_ptcl1, -cminus_ptcl1);
    p6_ptcl1_local1 = Vec(-aminus_ptcl1, bplus_ptcl1, -cminus_ptcl1);
    p7_ptcl1_local1 = Vec(-aminus_ptcl1, -bminus_ptcl1, -cminus_ptcl1);
    p8_ptcl1_local1 = Vec(aplus_ptcl1, -bminus_ptcl1, -cminus_ptcl1);

    // get local coordinates of particle 2
    REAL aplus_ptcl2 = p2->getAplus(); REAL aminus_ptcl2 = p2->getAminus(); 
    REAL bplus_ptcl2 = p2->getBplus(); REAL bminus_ptcl2 = p2->getBminus();  
    REAL cplus_ptcl2 = p2->getCplus(); REAL cminus_ptcl2 = p2->getCminus();

    p1_ptcl2_local2 = Vec(aplus_ptcl2, bplus_ptcl2, cplus_ptcl2);
    p2_ptcl2_local2 = Vec(-aminus_ptcl2, bplus_ptcl2, cplus_ptcl2);
    p3_ptcl2_local2 = Vec(-aminus_ptcl2, -bminus_ptcl2, cplus_ptcl2);
    p4_ptcl2_local2 = Vec(aplus_ptcl2, -bminus_ptcl2, cplus_ptcl2);
    p5_ptcl2_local2 = Vec(aplus_ptcl2, bplus_ptcl2, -cminus_ptcl2);
    p6_ptcl2_local2 = Vec(-aminus_ptcl2, bplus_ptcl2, -cminus_ptcl2);
    p7_ptcl2_local2 = Vec(-aminus_ptcl2, -bminus_ptcl2, -cminus_ptcl2);
    p8_ptcl2_local2 = Vec(aplus_ptcl2, -bminus_ptcl2, -cminus_ptcl2);

    // use particle 1 as the boundary particle and check if particle 2 is outside of particle 1 

    // change the 8 corner points of particle 2 to the local coordinates of particel 1
    Vec p1_ptcl2_global = p2->localToGlobal(p1_ptcl2_local2) + p2->getCurrPos();
    Vec p2_ptcl2_global = p2->localToGlobal(p2_ptcl2_local2) + p2->getCurrPos();
    Vec p3_ptcl2_global = p2->localToGlobal(p3_ptcl2_local2) + p2->getCurrPos();
    Vec p4_ptcl2_global = p2->localToGlobal(p4_ptcl2_local2) + p2->getCurrPos();
    Vec p5_ptcl2_global = p2->localToGlobal(p5_ptcl2_local2) + p2->getCurrPos();
    Vec p6_ptcl2_global = p2->localToGlobal(p6_ptcl2_local2) + p2->getCurrPos();
    Vec p7_ptcl2_global = p2->localToGlobal(p7_ptcl2_local2) + p2->getCurrPos();
    Vec p8_ptcl2_global = p2->localToGlobal(p8_ptcl2_local2) + p2->getCurrPos();

    Vec p1_ptcl2_local1 = p1->globalToLocal(p1_ptcl2_global - p1->getCurrPos());
    Vec p2_ptcl2_local1 = p1->globalToLocal(p2_ptcl2_global - p1->getCurrPos());
    Vec p3_ptcl2_local1 = p1->globalToLocal(p3_ptcl2_global - p1->getCurrPos());
    Vec p4_ptcl2_local1 = p1->globalToLocal(p4_ptcl2_global - p1->getCurrPos());
    Vec p5_ptcl2_local1 = p1->globalToLocal(p5_ptcl2_global - p1->getCurrPos());
    Vec p6_ptcl2_local1 = p1->globalToLocal(p6_ptcl2_global - p1->getCurrPos());
    Vec p7_ptcl2_local1 = p1->globalToLocal(p7_ptcl2_global - p1->getCurrPos());
    Vec p8_ptcl2_local1 = p1->globalToLocal(p8_ptcl2_global - p1->getCurrPos());



    // check if all 8 corner points of particle 2 are outside of particle 1
    // x+ boundary, i.e. a+ boundary
    if(p1_ptcl2_local1.getX()>aplus_ptcl1 && p2_ptcl2_local1.getX()>aplus_ptcl1
    && p3_ptcl2_local1.getX()>aplus_ptcl1 && p4_ptcl2_local1.getX()>aplus_ptcl1
    && p5_ptcl2_local1.getX()>aplus_ptcl1 && p6_ptcl2_local1.getX()>aplus_ptcl1
    && p7_ptcl2_local1.getX()>aplus_ptcl1 && p8_ptcl2_local1.getX()>aplus_ptcl1 )	// outside of x+ boundary
	return false;

    // x- boundary, i.e. a- boundary
    if(p1_ptcl2_local1.getX()<-aminus_ptcl1 && p2_ptcl2_local1.getX()<-aminus_ptcl1
    && p3_ptcl2_local1.getX()<-aminus_ptcl1 && p4_ptcl2_local1.getX()<-aminus_ptcl1
    && p5_ptcl2_local1.getX()<-aminus_ptcl1 && p6_ptcl2_local1.getX()<-aminus_ptcl1
    && p7_ptcl2_local1.getX()<-aminus_ptcl1 && p8_ptcl2_local1.getX()<-aminus_ptcl1 )	// outside of x- boundary
	return false;

    // y+ boundary, i.e. b+ boundary
    if(p1_ptcl2_local1.getY()>bplus_ptcl1 && p2_ptcl2_local1.getY()>bplus_ptcl1
    && p3_ptcl2_local1.getY()>bplus_ptcl1 && p4_ptcl2_local1.getY()>bplus_ptcl1
    && p5_ptcl2_local1.getY()>bplus_ptcl1 && p6_ptcl2_local1.getY()>bplus_ptcl1
    && p7_ptcl2_local1.getY()>bplus_ptcl1 && p8_ptcl2_local1.getY()>bplus_ptcl1 )	// outside of y+ boundary
	return false;

    // y- boundary, i.e. b- boundary
    if(p1_ptcl2_local1.getY()<-bminus_ptcl1 && p2_ptcl2_local1.getY()<-bminus_ptcl1
    && p3_ptcl2_local1.getY()<-bminus_ptcl1 && p4_ptcl2_local1.getY()<-bminus_ptcl1
    && p5_ptcl2_local1.getY()<-bminus_ptcl1 && p6_ptcl2_local1.getY()<-bminus_ptcl1
    && p7_ptcl2_local1.getY()<-bminus_ptcl1 && p8_ptcl2_local1.getY()<-bminus_ptcl1 )	// outside of y- boundary
	return false;

    // z+ boundary, i.e. c+ boundary
    if(p1_ptcl2_local1.getZ()>cplus_ptcl1 && p2_ptcl2_local1.getZ()>cplus_ptcl1
    && p3_ptcl2_local1.getZ()>cplus_ptcl1 && p4_ptcl2_local1.getZ()>cplus_ptcl1
    && p5_ptcl2_local1.getZ()>cplus_ptcl1 && p6_ptcl2_local1.getZ()>cplus_ptcl1
    && p7_ptcl2_local1.getZ()>cplus_ptcl1 && p8_ptcl2_local1.getZ()>cplus_ptcl1 )	// outside of z+ boundary
	return false;

    // z- boundary, i.e. c- boundary
    if(p1_ptcl2_local1.getZ()<-cminus_ptcl1 && p2_ptcl2_local1.getZ()<-cminus_ptcl1
    && p3_ptcl2_local1.getZ()<-cminus_ptcl1 && p4_ptcl2_local1.getZ()<-cminus_ptcl1
    && p5_ptcl2_local1.getZ()<-cminus_ptcl1 && p6_ptcl2_local1.getZ()<-cminus_ptcl1
    && p7_ptcl2_local1.getZ()<-cminus_ptcl1 && p8_ptcl2_local1.getZ()<-cminus_ptcl1 )	// outside of z- boundary
	return false;



    /////////////////////////////////////////////////////////////////////////////////////////////
    // use particle 2 as the boundary particle and check if particle 1 is outside of particle 2 

    // change the 8 corner points of particle 1 to the local coordinates of particel 2
    Vec p1_ptcl1_global = p1->localToGlobal(p1_ptcl1_local1) + p1->getCurrPos();
    Vec p2_ptcl1_global = p1->localToGlobal(p2_ptcl1_local1) + p1->getCurrPos();
    Vec p3_ptcl1_global = p1->localToGlobal(p3_ptcl1_local1) + p1->getCurrPos();
    Vec p4_ptcl1_global = p1->localToGlobal(p4_ptcl1_local1) + p1->getCurrPos();
    Vec p5_ptcl1_global = p1->localToGlobal(p5_ptcl1_local1) + p1->getCurrPos();
    Vec p6_ptcl1_global = p1->localToGlobal(p6_ptcl1_local1) + p1->getCurrPos();
    Vec p7_ptcl1_global = p1->localToGlobal(p7_ptcl1_local1) + p1->getCurrPos();
    Vec p8_ptcl1_global = p1->localToGlobal(p8_ptcl1_local1) + p1->getCurrPos();

    Vec p1_ptcl1_local2 = p2->globalToLocal(p1_ptcl1_global - p2->getCurrPos());
    Vec p2_ptcl1_local2 = p2->globalToLocal(p2_ptcl1_global - p2->getCurrPos());
    Vec p3_ptcl1_local2 = p2->globalToLocal(p3_ptcl1_global - p2->getCurrPos());
    Vec p4_ptcl1_local2 = p2->globalToLocal(p4_ptcl1_global - p2->getCurrPos());
    Vec p5_ptcl1_local2 = p2->globalToLocal(p5_ptcl1_global - p2->getCurrPos());
    Vec p6_ptcl1_local2 = p2->globalToLocal(p6_ptcl1_global - p2->getCurrPos());
    Vec p7_ptcl1_local2 = p2->globalToLocal(p7_ptcl1_global - p2->getCurrPos());
    Vec p8_ptcl1_local2 = p2->globalToLocal(p8_ptcl1_global - p2->getCurrPos());


    // check if all 8 corner points of particle 1 are outside of particle 2
    // x+ boundary, i.e. a+ boundary
    if(p1_ptcl1_local2.getX()>aplus_ptcl2 && p2_ptcl1_local2.getX()>aplus_ptcl2
    && p3_ptcl1_local2.getX()>aplus_ptcl2 && p4_ptcl1_local2.getX()>aplus_ptcl2
    && p5_ptcl1_local2.getX()>aplus_ptcl2 && p6_ptcl1_local2.getX()>aplus_ptcl2
    && p7_ptcl1_local2.getX()>aplus_ptcl2 && p8_ptcl1_local2.getX()>aplus_ptcl2 )	// outside of x+ boundary
	return false;

    // x- boundary, i.e. a- boundary
    if(p1_ptcl1_local2.getX()<-aminus_ptcl2 && p2_ptcl1_local2.getX()<-aminus_ptcl2
    && p3_ptcl1_local2.getX()<-aminus_ptcl2 && p4_ptcl1_local2.getX()<-aminus_ptcl2
    && p5_ptcl1_local2.getX()<-aminus_ptcl2 && p6_ptcl1_local2.getX()<-aminus_ptcl2
    && p7_ptcl1_local2.getX()<-aminus_ptcl2 && p8_ptcl1_local2.getX()<-aminus_ptcl2 )	// outside of x- boundary
	return false;

    // y+ boundary, i.e. b+ boundary
    if(p1_ptcl1_local2.getY()>bplus_ptcl2 && p2_ptcl1_local2.getY()>bplus_ptcl2
    && p3_ptcl1_local2.getY()>bplus_ptcl2 && p4_ptcl1_local2.getY()>bplus_ptcl2
    && p5_ptcl1_local2.getY()>bplus_ptcl2 && p6_ptcl1_local2.getY()>bplus_ptcl2
    && p7_ptcl1_local2.getY()>bplus_ptcl2 && p8_ptcl1_local2.getY()>bplus_ptcl2 )	// outside of y+ boundary
	return false;

    // y- boundary, i.e. b- boundary
    if(p1_ptcl1_local2.getY()<-bminus_ptcl2 && p2_ptcl1_local2.getY()<-bminus_ptcl2
    && p3_ptcl1_local2.getY()<-bminus_ptcl2 && p4_ptcl1_local2.getY()<-bminus_ptcl2
    && p5_ptcl1_local2.getY()<-bminus_ptcl2 && p6_ptcl1_local2.getY()<-bminus_ptcl2
    && p7_ptcl1_local2.getY()<-bminus_ptcl2 && p8_ptcl1_local2.getY()<-bminus_ptcl2 )	// outside of y- boundary
	return false;

    // z+ boundary, i.e. c+ boundary
    if(p1_ptcl1_local2.getZ()>cplus_ptcl2 && p2_ptcl1_local2.getZ()>cplus_ptcl2
    && p3_ptcl1_local2.getZ()>cplus_ptcl2 && p4_ptcl1_local2.getZ()>cplus_ptcl2
    && p5_ptcl1_local2.getZ()>cplus_ptcl2 && p6_ptcl1_local2.getZ()>cplus_ptcl2
    && p7_ptcl1_local2.getZ()>cplus_ptcl2 && p8_ptcl1_local2.getZ()>cplus_ptcl2 )	// outside of z+ boundary
	return false;

    // z- boundary, i.e. c- boundary
    if(p1_ptcl1_local2.getZ()<-cminus_ptcl2 && p2_ptcl1_local2.getZ()<-cminus_ptcl2
    && p3_ptcl1_local2.getZ()<-cminus_ptcl2 && p4_ptcl1_local2.getZ()<-cminus_ptcl2
    && p5_ptcl1_local2.getZ()<-cminus_ptcl2 && p6_ptcl1_local2.getZ()<-cminus_ptcl2
    && p7_ptcl1_local2.getZ()<-cminus_ptcl2 && p8_ptcl1_local2.getZ()<-cminus_ptcl2 )	// outside of z- boundary
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
    Vec aplus_ptcl2_local1  = 0.25*(p1_ptcl2_local1 + p4_ptcl2_local1 + p8_ptcl2_local1 + p5_ptcl2_local1);
    Vec aminus_ptcl2_local1 = 0.25*(p2_ptcl2_local1 + p3_ptcl2_local1 + p7_ptcl2_local1 + p6_ptcl2_local1);
    Vec bplus_ptcl2_local1  = 0.25*(p1_ptcl2_local1 + p2_ptcl2_local1 + p6_ptcl2_local1 + p5_ptcl2_local1);
    Vec bminus_ptcl2_local1 = 0.25*(p4_ptcl2_local1 + p3_ptcl2_local1 + p7_ptcl2_local1 + p8_ptcl2_local1);
    Vec cplus_ptcl2_local1  = 0.25*(p1_ptcl2_local1 + p2_ptcl2_local1 + p3_ptcl2_local1 + p4_ptcl2_local1);
    Vec cminus_ptcl2_local1 = 0.25*(p5_ptcl2_local1 + p6_ptcl2_local1 + p7_ptcl2_local1 + p8_ptcl2_local1);

    // check with x+ boundary of particle 1
//    if(p1_ptcl2_local1.getX()<aplus_ptcl1 && p2_ptcl2_local1.getX()<aplus_ptcl1
//    && p3_ptcl2_local1.getX()<aplus_ptcl1 && p4_ptcl2_local1.getX()<aplus_ptcl1
//    && p5_ptcl2_local1.getX()<aplus_ptcl1 && p6_ptcl2_local1.getX()<aplus_ptcl1
//    && p7_ptcl2_local1.getX()<aplus_ptcl1 && p8_ptcl2_local1.getX()<aplus_ptcl1 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl2_local1.getX()<aminus_ptcl2_local1.getX())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl2_local1.getX()<bminus_ptcl2_local1.getX())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl2_local1.getX()<cminus_ptcl2_local1.getX())
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
//    if(p1_ptcl2_local1.getX()>-aminus_ptcl1 && p2_ptcl2_local1.getX()>-aminus_ptcl1
//    && p3_ptcl2_local1.getX()>-aminus_ptcl1 && p4_ptcl2_local1.getX()>-aminus_ptcl1
//    && p5_ptcl2_local1.getX()>-aminus_ptcl1 && p6_ptcl2_local1.getX()>-aminus_ptcl1
//    && p7_ptcl2_local1.getX()>-aminus_ptcl1 && p8_ptcl2_local1.getX()>-aminus_ptcl1 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl2_local1.getX()>aminus_ptcl2_local1.getX())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl2_local1.getX()>bminus_ptcl2_local1.getX())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl2_local1.getX()>cminus_ptcl2_local1.getX())
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
//    if(p1_ptcl2_local1.getY()<bplus_ptcl1 && p2_ptcl2_local1.getY()<bplus_ptcl1
//    && p3_ptcl2_local1.getY()<bplus_ptcl1 && p4_ptcl2_local1.getY()<bplus_ptcl1
//    && p5_ptcl2_local1.getY()<bplus_ptcl1 && p6_ptcl2_local1.getY()<bplus_ptcl1
//    && p7_ptcl2_local1.getY()<bplus_ptcl1 && p8_ptcl2_local1.getY()<bplus_ptcl1 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl2_local1.getY()<aminus_ptcl2_local1.getY())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl2_local1.getY()<bminus_ptcl2_local1.getY())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl2_local1.getY()<cminus_ptcl2_local1.getY())
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
//    if(p1_ptcl2_local1.getY()>-bminus_ptcl1 && p2_ptcl2_local1.getY()>-bminus_ptcl1
//    && p3_ptcl2_local1.getY()>-bminus_ptcl1 && p4_ptcl2_local1.getY()>-bminus_ptcl1
//    && p5_ptcl2_local1.getY()>-bminus_ptcl1 && p6_ptcl2_local1.getY()>-bminus_ptcl1
//    && p7_ptcl2_local1.getY()>-bminus_ptcl1 && p8_ptcl2_local1.getY()>-bminus_ptcl1 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl2_local1.getY()>aminus_ptcl2_local1.getY())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl2_local1.getY()>bminus_ptcl2_local1.getY())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl2_local1.getY()>cminus_ptcl2_local1.getY())
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
//    if(p1_ptcl2_local1.getZ()<cplus_ptcl1 && p2_ptcl2_local1.getZ()<cplus_ptcl1
//    && p3_ptcl2_local1.getZ()<cplus_ptcl1 && p4_ptcl2_local1.getZ()<cplus_ptcl1
//    && p5_ptcl2_local1.getZ()<cplus_ptcl1 && p6_ptcl2_local1.getZ()<cplus_ptcl1
//    && p7_ptcl2_local1.getZ()<cplus_ptcl1 && p8_ptcl2_local1.getZ()<cplus_ptcl1 )	// do nothing
//	;
//    else{
    {

	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl2_local1.getZ()<aminus_ptcl2_local1.getZ())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl2_local1.getZ()<bminus_ptcl2_local1.getZ())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl2_local1.getZ()<cminus_ptcl2_local1.getZ())
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
//    if(p1_ptcl2_local1.getZ()>-cminus_ptcl1 && p2_ptcl2_local1.getZ()>-cminus_ptcl1
//    && p3_ptcl2_local1.getZ()>-cminus_ptcl1 && p4_ptcl2_local1.getZ()>-cminus_ptcl1
//    && p5_ptcl2_local1.getZ()>-cminus_ptcl1 && p6_ptcl2_local1.getZ()>-cminus_ptcl1
//    && p7_ptcl2_local1.getZ()>-cminus_ptcl1 && p8_ptcl2_local1.getZ()>-cminus_ptcl1 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl2_local1.getZ()>aminus_ptcl2_local1.getZ())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl2_local1.getZ()>bminus_ptcl2_local1.getZ())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl2_local1.getZ()>cminus_ptcl2_local1.getZ())
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
    Vec aplus_ptcl1_local2  = 0.25*(p1_ptcl1_local2 + p4_ptcl1_local2 + p8_ptcl1_local2 + p5_ptcl1_local2);
    Vec aminus_ptcl1_local2 = 0.25*(p2_ptcl1_local2 + p3_ptcl1_local2 + p7_ptcl1_local2 + p6_ptcl1_local2);
    Vec bplus_ptcl1_local2  = 0.25*(p1_ptcl1_local2 + p2_ptcl1_local2 + p6_ptcl1_local2 + p5_ptcl1_local2);
    Vec bminus_ptcl1_local2 = 0.25*(p4_ptcl1_local2 + p3_ptcl1_local2 + p7_ptcl1_local2 + p8_ptcl1_local2);
    Vec cplus_ptcl1_local2  = 0.25*(p1_ptcl1_local2 + p2_ptcl1_local2 + p3_ptcl1_local2 + p4_ptcl1_local2);
    Vec cminus_ptcl1_local2 = 0.25*(p5_ptcl1_local2 + p6_ptcl1_local2 + p7_ptcl1_local2 + p8_ptcl1_local2);

    // check with x+ boundary of particle 2
//    if(p1_ptcl1_local2.getX()<aplus_ptcl2 && p2_ptcl1_local2.getX()<aplus_ptcl2
//    && p3_ptcl1_local2.getX()<aplus_ptcl2 && p4_ptcl1_local2.getX()<aplus_ptcl2
//    && p5_ptcl1_local2.getX()<aplus_ptcl2 && p6_ptcl1_local2.getX()<aplus_ptcl2
//    && p7_ptcl1_local2.getX()<aplus_ptcl2 && p8_ptcl1_local2.getX()<aplus_ptcl2 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl1_local2.getX()<aminus_ptcl1_local2.getX())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl1_local2.getX()<bminus_ptcl1_local2.getX())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl1_local2.getX()<cminus_ptcl1_local2.getX())
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
//    if(p1_ptcl1_local2.getX()>-aminus_ptcl2 && p2_ptcl1_local2.getX()>-aminus_ptcl2
//    && p3_ptcl1_local2.getX()>-aminus_ptcl2 && p4_ptcl1_local2.getX()>-aminus_ptcl2
//    && p5_ptcl1_local2.getX()>-aminus_ptcl2 && p6_ptcl1_local2.getX()>-aminus_ptcl2
//    && p7_ptcl1_local2.getX()>-aminus_ptcl2 && p8_ptcl1_local2.getX()>-aminus_ptcl2 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl1_local2.getX()>aminus_ptcl1_local2.getX())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl1_local2.getX()>bminus_ptcl1_local2.getX())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl1_local2.getX()>cminus_ptcl1_local2.getX())
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
//    if(p1_ptcl1_local2.getY()<bplus_ptcl2 && p2_ptcl1_local2.getY()<bplus_ptcl2
//    && p3_ptcl1_local2.getY()<bplus_ptcl2 && p4_ptcl1_local2.getY()<bplus_ptcl2
//    && p5_ptcl1_local2.getY()<bplus_ptcl2 && p6_ptcl1_local2.getY()<bplus_ptcl2
//    && p7_ptcl1_local2.getY()<bplus_ptcl2 && p8_ptcl1_local2.getY()<bplus_ptcl2 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl1_local2.getY()<aminus_ptcl1_local2.getY())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl1_local2.getY()<bminus_ptcl1_local2.getY())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl1_local2.getY()<cminus_ptcl1_local2.getY())
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
//    if(p1_ptcl1_local2.getY()>-bminus_ptcl2 && p2_ptcl1_local2.getY()>-bminus_ptcl2
//    && p3_ptcl1_local2.getY()>-bminus_ptcl2 && p4_ptcl1_local2.getY()>-bminus_ptcl2
//    && p5_ptcl1_local2.getY()>-bminus_ptcl2 && p6_ptcl1_local2.getY()>-bminus_ptcl2
//    && p7_ptcl1_local2.getY()>-bminus_ptcl2 && p8_ptcl1_local2.getY()>-bminus_ptcl2 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl1_local2.getY()>aminus_ptcl1_local2.getY())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl1_local2.getY()>bminus_ptcl1_local2.getY())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl1_local2.getY()>cminus_ptcl1_local2.getY())
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
//    if(p1_ptcl1_local2.getZ()<cplus_ptcl2 && p2_ptcl1_local2.getZ()<cplus_ptcl2
//    && p3_ptcl1_local2.getZ()<cplus_ptcl2 && p4_ptcl1_local2.getZ()<cplus_ptcl2
//    && p5_ptcl1_local2.getZ()<cplus_ptcl2 && p6_ptcl1_local2.getZ()<cplus_ptcl2
//    && p7_ptcl1_local2.getZ()<cplus_ptcl2 && p8_ptcl1_local2.getZ()<cplus_ptcl2 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl1_local2.getZ()<aminus_ptcl1_local2.getZ())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl1_local2.getZ()<bminus_ptcl1_local2.getZ())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl1_local2.getZ()<cminus_ptcl1_local2.getZ())
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
//    if(p1_ptcl1_local2.getZ()>-cminus_ptcl2 && p2_ptcl1_local2.getZ()>-cminus_ptcl2
//    && p3_ptcl1_local2.getZ()>-cminus_ptcl2 && p4_ptcl1_local2.getZ()>-cminus_ptcl2
//    && p5_ptcl1_local2.getZ()>-cminus_ptcl2 && p6_ptcl1_local2.getZ()>-cminus_ptcl2
//    && p7_ptcl1_local2.getZ()>-cminus_ptcl2 && p8_ptcl1_local2.getZ()>-cminus_ptcl2 )	// do nothing
//	;
//    else{
    {
	// check which octant is possible contact with this boundary
	int xsign, ysign, zsign;

	if(aplus_ptcl1_local2.getZ()>aminus_ptcl1_local2.getZ())
		xsign = 1;
	else
		xsign = -1;
	if(bplus_ptcl1_local2.getZ()>bminus_ptcl1_local2.getZ())
		ysign = 1;
	else 
		ysign = -1;
	if(cplus_ptcl1_local2.getZ()>cminus_ptcl1_local2.getZ())
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
    Vec p1_ptcl1_local1, p2_ptcl1_local1, p3_ptcl1_local1, p4_ptcl1_local1, p5_ptcl1_local1, p6_ptcl1_local1, p7_ptcl1_local1, p8_ptcl1_local1;
    Vec p1_ptcl2_local2, p2_ptcl2_local2, p3_ptcl2_local2, p4_ptcl2_local2, p5_ptcl2_local2, p6_ptcl2_local2, p7_ptcl2_local2, p8_ptcl2_local2;

    // get local coordinates of particle 1
    REAL aplus_ptcl1 = p1->getAplus(); REAL aminus_ptcl1 = p1->getAminus(); 
    REAL bplus_ptcl1 = p1->getBplus(); REAL bminus_ptcl1 = p1->getBminus();  
    REAL cplus_ptcl1 = p1->getCplus(); REAL cminus_ptcl1 = p1->getCminus();

    p1_ptcl1_local1 = Vec(aplus_ptcl1, bplus_ptcl1, cplus_ptcl1);
    p2_ptcl1_local1 = Vec(-aminus_ptcl1, bplus_ptcl1, cplus_ptcl1);
    p3_ptcl1_local1 = Vec(-aminus_ptcl1, -bminus_ptcl1, cplus_ptcl1);
    p4_ptcl1_local1 = Vec(aplus_ptcl1, -bminus_ptcl1, cplus_ptcl1);
    p5_ptcl1_local1 = Vec(aplus_ptcl1, bplus_ptcl1, -cminus_ptcl1);
    p6_ptcl1_local1 = Vec(-aminus_ptcl1, bplus_ptcl1, -cminus_ptcl1);
    p7_ptcl1_local1 = Vec(-aminus_ptcl1, -bminus_ptcl1, -cminus_ptcl1);
    p8_ptcl1_local1 = Vec(aplus_ptcl1, -bminus_ptcl1, -cminus_ptcl1);

    // get local coordinates of particle 2
    REAL aplus_ptcl2 = p2->getAplus(); REAL aminus_ptcl2 = p2->getAminus(); 
    REAL bplus_ptcl2 = p2->getBplus(); REAL bminus_ptcl2 = p2->getBminus();  
    REAL cplus_ptcl2 = p2->getCplus(); REAL cminus_ptcl2 = p2->getCminus();

    p1_ptcl2_local2 = Vec(aplus_ptcl2, bplus_ptcl2, cplus_ptcl2);
    p2_ptcl2_local2 = Vec(-aminus_ptcl2, bplus_ptcl2, cplus_ptcl2);
    p3_ptcl2_local2 = Vec(-aminus_ptcl2, -bminus_ptcl2, cplus_ptcl2);
    p4_ptcl2_local2 = Vec(aplus_ptcl2, -bminus_ptcl2, cplus_ptcl2);
    p5_ptcl2_local2 = Vec(aplus_ptcl2, bplus_ptcl2, -cminus_ptcl2);
    p6_ptcl2_local2 = Vec(-aminus_ptcl2, bplus_ptcl2, -cminus_ptcl2);
    p7_ptcl2_local2 = Vec(-aminus_ptcl2, -bminus_ptcl2, -cminus_ptcl2);
    p8_ptcl2_local2 = Vec(aplus_ptcl2, -bminus_ptcl2, -cminus_ptcl2);


    // change the 8 corner points of particle 2 to the local coordinates of particel 1
    Vec p1_ptcl2_global = p2->localToGlobal(p1_ptcl2_local2) + p2->getCurrPos();
    Vec p2_ptcl2_global = p2->localToGlobal(p2_ptcl2_local2) + p2->getCurrPos();
    Vec p3_ptcl2_global = p2->localToGlobal(p3_ptcl2_local2) + p2->getCurrPos();
    Vec p4_ptcl2_global = p2->localToGlobal(p4_ptcl2_local2) + p2->getCurrPos();
    Vec p5_ptcl2_global = p2->localToGlobal(p5_ptcl2_local2) + p2->getCurrPos();
    Vec p6_ptcl2_global = p2->localToGlobal(p6_ptcl2_local2) + p2->getCurrPos();
    Vec p7_ptcl2_global = p2->localToGlobal(p7_ptcl2_local2) + p2->getCurrPos();
    Vec p8_ptcl2_global = p2->localToGlobal(p8_ptcl2_local2) + p2->getCurrPos();

    Vec p1_ptcl2_local1 = p1->globalToLocal(p1_ptcl2_global - p1->getCurrPos());
    Vec p2_ptcl2_local1 = p1->globalToLocal(p2_ptcl2_global - p1->getCurrPos());
    Vec p3_ptcl2_local1 = p1->globalToLocal(p3_ptcl2_global - p1->getCurrPos());
    Vec p4_ptcl2_local1 = p1->globalToLocal(p4_ptcl2_global - p1->getCurrPos());
    Vec p5_ptcl2_local1 = p1->globalToLocal(p5_ptcl2_global - p1->getCurrPos());
    Vec p6_ptcl2_local1 = p1->globalToLocal(p6_ptcl2_global - p1->getCurrPos());
    Vec p7_ptcl2_local1 = p1->globalToLocal(p7_ptcl2_global - p1->getCurrPos());
    Vec p8_ptcl2_local1 = p1->globalToLocal(p8_ptcl2_global - p1->getCurrPos());
    
    Vec center_ptcl2_local1 = 0.125*( p1_ptcl2_local1 + p2_ptcl2_local1 + p3_ptcl2_local1 + p4_ptcl2_local1
				    + p5_ptcl2_local1 + p6_ptcl2_local1 + p7_ptcl2_local1 + p8_ptcl2_local1 );

    // change the 8 corner points of particle 1 to the local coordinates of particel 2
    Vec p1_ptcl1_global = p1->localToGlobal(p1_ptcl1_local1) + p1->getCurrPos();
    Vec p2_ptcl1_global = p1->localToGlobal(p2_ptcl1_local1) + p1->getCurrPos();
    Vec p3_ptcl1_global = p1->localToGlobal(p3_ptcl1_local1) + p1->getCurrPos();
    Vec p4_ptcl1_global = p1->localToGlobal(p4_ptcl1_local1) + p1->getCurrPos();
    Vec p5_ptcl1_global = p1->localToGlobal(p5_ptcl1_local1) + p1->getCurrPos();
    Vec p6_ptcl1_global = p1->localToGlobal(p6_ptcl1_local1) + p1->getCurrPos();
    Vec p7_ptcl1_global = p1->localToGlobal(p7_ptcl1_local1) + p1->getCurrPos();
    Vec p8_ptcl1_global = p1->localToGlobal(p8_ptcl1_local1) + p1->getCurrPos();

    Vec p1_ptcl1_local2 = p2->globalToLocal(p1_ptcl1_global - p2->getCurrPos());
    Vec p2_ptcl1_local2 = p2->globalToLocal(p2_ptcl1_global - p2->getCurrPos());
    Vec p3_ptcl1_local2 = p2->globalToLocal(p3_ptcl1_global - p2->getCurrPos());
    Vec p4_ptcl1_local2 = p2->globalToLocal(p4_ptcl1_global - p2->getCurrPos());
    Vec p5_ptcl1_local2 = p2->globalToLocal(p5_ptcl1_global - p2->getCurrPos());
    Vec p6_ptcl1_local2 = p2->globalToLocal(p6_ptcl1_global - p2->getCurrPos());
    Vec p7_ptcl1_local2 = p2->globalToLocal(p7_ptcl1_global - p2->getCurrPos());
    Vec p8_ptcl1_local2 = p2->globalToLocal(p8_ptcl1_global - p2->getCurrPos());

    Vec center_ptcl1_local2 = 0.125*( p1_ptcl1_local2 + p2_ptcl1_local2 + p3_ptcl1_local2 + p4_ptcl1_local2
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
	    px = p1_ptcl2_local1.getX();
	    py = p1_ptcl2_local1.getY();
	    pz = p1_ptcl2_local1.getZ();
 
	}
	else if (num_oct == 2){
	    // check p2_ptcl2_local1
	    px = p2_ptcl2_local1.getX();
	    py = p2_ptcl2_local1.getY();
	    pz = p2_ptcl2_local1.getZ();
	}
	else if (num_oct == 3){
	    // check p3_ptcl2_local1
	    px = p3_ptcl2_local1.getX();
	    py = p3_ptcl2_local1.getY();
	    pz = p3_ptcl2_local1.getZ();
	}
	else if (num_oct == 4){
	    // check p4_ptcl2_local1
	    px = p4_ptcl2_local1.getX();
	    py = p4_ptcl2_local1.getY();
	    pz = p4_ptcl2_local1.getZ();
	}
	else if (num_oct == 5){
	    // check p5_ptcl2_local1
	    px = p5_ptcl2_local1.getX();
	    py = p5_ptcl2_local1.getY();
	    pz = p5_ptcl2_local1.getZ();
	}
	else if (num_oct == 6){
	    // check p6_ptcl2_local1
	    px = p6_ptcl2_local1.getX();
	    py = p6_ptcl2_local1.getY();
	    pz = p6_ptcl2_local1.getZ();
	}
	else if (num_oct == 7){
	    // check p7_ptcl2_local1
	    px = p7_ptcl2_local1.getX();
	    py = p7_ptcl2_local1.getY();
	    pz = p7_ptcl2_local1.getZ();
	}
	else if (num_oct == 8){
	    // check p8_ptcl2_local1
	    px = p8_ptcl2_local1.getX();
	    py = p8_ptcl2_local1.getY();
	    pz = p8_ptcl2_local1.getZ();
	}
	else{	
	    // chenck center_ptcl2_local1
	    px = center_ptcl2_local1.getX();
	    py = center_ptcl2_local1.getY();
	    pz = center_ptcl2_local1.getZ();	
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
	    px = p1_ptcl1_local2.getX();
	    py = p1_ptcl1_local2.getY();
	    pz = p1_ptcl1_local2.getZ();
 
	}
	else if (num_oct == 2){
	    // check p2_ptcl1_local2
	    px = p2_ptcl1_local2.getX();
	    py = p2_ptcl1_local2.getY();
	    pz = p2_ptcl1_local2.getZ();
	}
	else if (num_oct == 3){
	    // check p3_ptcl1_local2
	    px = p3_ptcl1_local2.getX();
	    py = p3_ptcl1_local2.getY();
	    pz = p3_ptcl1_local2.getZ();
	}
	else if (num_oct == 4){
	    // check p4_ptcl1_local2
	    px = p4_ptcl1_local2.getX();
	    py = p4_ptcl1_local2.getY();
	    pz = p4_ptcl1_local2.getZ();
	}
	else if (num_oct == 5){
	    // check p5_ptcl1_local2
	    px = p5_ptcl1_local2.getX();
	    py = p5_ptcl1_local2.getY();
	    pz = p5_ptcl1_local2.getZ();
	}
	else if (num_oct == 6){
	    // check p6_ptcl1_local2
	    px = p6_ptcl1_local2.getX();
	    py = p6_ptcl1_local2.getY();
	    pz = p6_ptcl1_local2.getZ();
	}
	else if (num_oct == 7){
	    // check p7_ptcl1_local2
	    px = p7_ptcl1_local2.getX();
	    py = p7_ptcl1_local2.getY();
	    pz = p7_ptcl1_local2.getZ();
	}
	else if (num_oct == 8){
	    // check p8_ptcl1_local2
	    px = p8_ptcl1_local2.getX();
	    py = p8_ptcl1_local2.getY();
	    pz = p8_ptcl1_local2.getZ();
	}
	else{
	    // check center_ptcl1_local2
	    px = center_ptcl1_local2.getX();
	    py = center_ptcl1_local2.getY();
	    pz = center_ptcl1_local2.getZ();
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
	    px = p1_ptcl2_local1.getX();
	    py = p1_ptcl2_local1.getY();
	    pz = p1_ptcl2_local1.getZ();
 
	}
	else if (num_oct == 2){
	    // check p2_ptcl2_local1
	    px = p2_ptcl2_local1.getX();
	    py = p2_ptcl2_local1.getY();
	    pz = p2_ptcl2_local1.getZ();
	}
	else if (num_oct == 3){
	    // check p3_ptcl2_local1
	    px = p3_ptcl2_local1.getX();
	    py = p3_ptcl2_local1.getY();
	    pz = p3_ptcl2_local1.getZ();
	}
	else if (num_oct == 4){
	    // check p4_ptcl2_local1
	    px = p4_ptcl2_local1.getX();
	    py = p4_ptcl2_local1.getY();
	    pz = p4_ptcl2_local1.getZ();
	}
	else if (num_oct == 5){
	    // check p5_ptcl2_local1
	    px = p5_ptcl2_local1.getX();
	    py = p5_ptcl2_local1.getY();
	    pz = p5_ptcl2_local1.getZ();
	}
	else if (num_oct == 6){
	    // check p6_ptcl2_local1
	    px = p6_ptcl2_local1.getX();
	    py = p6_ptcl2_local1.getY();
	    pz = p6_ptcl2_local1.getZ();
	}
	else if (num_oct == 7){
	    // check p7_ptcl2_local1
	    px = p7_ptcl2_local1.getX();
	    py = p7_ptcl2_local1.getY();
	    pz = p7_ptcl2_local1.getZ();
	}
	else{
	    // check p8_ptcl2_local1
	    px = p8_ptcl2_local1.getX();
	    py = p8_ptcl2_local1.getY();
	    pz = p8_ptcl2_local1.getZ();
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
	    px = p1_ptcl1_local2.getX();
	    py = p1_ptcl1_local2.getY();
	    pz = p1_ptcl1_local2.getZ();
 
	}
	else if (num_oct == 2){
	    // check p2_ptcl1_local2
	    px = p2_ptcl1_local2.getX();
	    py = p2_ptcl1_local2.getY();
	    pz = p2_ptcl1_local2.getZ();
	}
	else if (num_oct == 3){
	    // check p3_ptcl1_local2
	    px = p3_ptcl1_local2.getX();
	    py = p3_ptcl1_local2.getY();
	    pz = p3_ptcl1_local2.getZ();
	}
	else if (num_oct == 4){
	    // check p4_ptcl1_local2
	    px = p4_ptcl1_local2.getX();
	    py = p4_ptcl1_local2.getY();
	    pz = p4_ptcl1_local2.getZ();
	}
	else if (num_oct == 5){
	    // check p5_ptcl1_local2
	    px = p5_ptcl1_local2.getX();
	    py = p5_ptcl1_local2.getY();
	    pz = p5_ptcl1_local2.getZ();
	}
	else if (num_oct == 6){
	    // check p6_ptcl1_local2
	    px = p6_ptcl1_local2.getX();
	    py = p6_ptcl1_local2.getY();
	    pz = p6_ptcl1_local2.getZ();
	}
	else if (num_oct == 7){
	    // check p7_ptcl1_local2
	    px = p7_ptcl1_local2.getX();
	    py = p7_ptcl1_local2.getY();
	    pz = p7_ptcl1_local2.getZ();
	}
	else{
	    // check p8_ptcl1_local2
	    px = p8_ptcl1_local2.getX();
	    py = p8_ptcl1_local2.getY();
	    pz = p8_ptcl1_local2.getZ();
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
    if( center_ptcl2_local1.getX() <= aplus_ptcl1 && center_ptcl2_local1.getX() >= -aminus_ptcl1 
     && center_ptcl2_local1.getY() <= bplus_ptcl1 && center_ptcl2_local1.getY() >= -bminus_ptcl1 
     && center_ptcl2_local1.getZ() <= cplus_ptcl1 && center_ptcl2_local1.getZ() >= -cminus_ptcl1 ){	// inside of particle 1 box

	for(int ii = 1; ii !=9; ii++){
	    octant_ptcl1.push_back(ii);
	    octant_ptcl2.push_back(ii); 
	}
    }

    // check center_ptcl1_local2
    if( center_ptcl1_local2.getX() <= aplus_ptcl2 && center_ptcl1_local2.getX() >= -aminus_ptcl2 
     && center_ptcl1_local2.getY() <= bplus_ptcl2 && center_ptcl1_local2.getY() >= -bminus_ptcl2 
     && center_ptcl1_local2.getZ() <= cplus_ptcl2 && center_ptcl1_local2.getZ() >= -cminus_ptcl2 ){	// inside of particle 2 box

	for(int ii = 1; ii !=9; ii++){
	    octant_ptcl1.push_back(ii);
	    octant_ptcl2.push_back(ii); 
	}
    }




    // delete dublicated octant number in the octant Vectors
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

    Vec v[2];
    bool b1, b2;
    Vec local1_point1, local2_point2;



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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() >= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() >= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() <= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() >= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() >= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() <= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() <= 0 )
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

    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    if(b1 && b2 
    && local1_point1.getX() >= 0 && local1_point1.getY() <= 0 && local1_point1.getZ() <= 0
    && local2_point2.getX() >= 0 && local2_point2.getY() <= 0 && local2_point2.getZ() <= 0 )
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

    Vec p9_ptcl1_local2  = p1_ptcl1_local2 + aplus_ratio_ptcl1*(p2_ptcl1_local2 - p1_ptcl1_local2);
    Vec p10_ptcl1_local2 = p4_ptcl1_local2 + aplus_ratio_ptcl1*(p3_ptcl1_local2 - p4_ptcl1_local2);
    Vec p11_ptcl1_local2 = p8_ptcl1_local2 + aplus_ratio_ptcl1*(p7_ptcl1_local2 - p8_ptcl1_local2);
    Vec p12_ptcl1_local2 = p5_ptcl1_local2 + aplus_ratio_ptcl1*(p6_ptcl1_local2 - p5_ptcl1_local2);

    Vec p13_ptcl1_local2 = p1_ptcl1_local2 + bplus_ratio_ptcl1*(p4_ptcl1_local2 - p1_ptcl1_local2);
    Vec p14_ptcl1_local2 = p2_ptcl1_local2 + bplus_ratio_ptcl1*(p3_ptcl1_local2 - p2_ptcl1_local2);
    Vec p15_ptcl1_local2 = p6_ptcl1_local2 + bplus_ratio_ptcl1*(p7_ptcl1_local2 - p6_ptcl1_local2);
    Vec p16_ptcl1_local2 = p5_ptcl1_local2 + bplus_ratio_ptcl1*(p8_ptcl1_local2 - p5_ptcl1_local2);

    Vec p17_ptcl1_local2 = p1_ptcl1_local2 + cplus_ratio_ptcl1*(p5_ptcl1_local2 - p1_ptcl1_local2);
    Vec p18_ptcl1_local2 = p2_ptcl1_local2 + cplus_ratio_ptcl1*(p6_ptcl1_local2 - p2_ptcl1_local2);
    Vec p19_ptcl1_local2 = p3_ptcl1_local2 + cplus_ratio_ptcl1*(p7_ptcl1_local2 - p3_ptcl1_local2);
    Vec p20_ptcl1_local2 = p4_ptcl1_local2 + cplus_ratio_ptcl1*(p8_ptcl1_local2 - p4_ptcl1_local2);

    Vec p21_ptcl1_local2 = p13_ptcl1_local2 + aplus_ratio_ptcl1*(p14_ptcl1_local2 - p13_ptcl1_local2);
    Vec p22_ptcl1_local2 = p9_ptcl1_local2  + cplus_ratio_ptcl1*(p12_ptcl1_local2 -  p9_ptcl1_local2);
    Vec p23_ptcl1_local2 = p14_ptcl1_local2 + cplus_ratio_ptcl1*(p15_ptcl1_local2 - p14_ptcl1_local2);
    Vec p24_ptcl1_local2 = p10_ptcl1_local2 + cplus_ratio_ptcl1*(p11_ptcl1_local2 - p10_ptcl1_local2);
    Vec p25_ptcl1_local2 = p13_ptcl1_local2 + cplus_ratio_ptcl1*(p16_ptcl1_local2 - p13_ptcl1_local2);
    Vec p26_ptcl1_local2 = p16_ptcl1_local2 + aplus_ratio_ptcl1*(p15_ptcl1_local2 - p16_ptcl1_local2);
    Vec p27_ptcl1_local2 = p21_ptcl1_local2 + cplus_ratio_ptcl1*(p26_ptcl1_local2 - p21_ptcl1_local2);

    // get the 27 coordinates of the sub box points of particle 2 in local1
    REAL aplus_ratio_ptcl2 = aplus_ptcl2/(aplus_ptcl2+aminus_ptcl2);
    REAL bplus_ratio_ptcl2 = bplus_ptcl2/(bplus_ptcl2+bminus_ptcl2);
    REAL cplus_ratio_ptcl2 = cplus_ptcl2/(cplus_ptcl2+cminus_ptcl2);

    Vec p9_ptcl2_local1  = p1_ptcl2_local1 + aplus_ratio_ptcl2*(p2_ptcl2_local1 - p1_ptcl2_local1);
    Vec p10_ptcl2_local1 = p4_ptcl2_local1 + aplus_ratio_ptcl2*(p3_ptcl2_local1 - p4_ptcl2_local1);
    Vec p11_ptcl2_local1 = p8_ptcl2_local1 + aplus_ratio_ptcl2*(p7_ptcl2_local1 - p8_ptcl2_local1);
    Vec p12_ptcl2_local1 = p5_ptcl2_local1 + aplus_ratio_ptcl2*(p6_ptcl2_local1 - p5_ptcl2_local1);

    Vec p13_ptcl2_local1 = p1_ptcl2_local1 + bplus_ratio_ptcl2*(p4_ptcl2_local1 - p1_ptcl2_local1);
    Vec p14_ptcl2_local1 = p2_ptcl2_local1 + bplus_ratio_ptcl2*(p3_ptcl2_local1 - p2_ptcl2_local1);
    Vec p15_ptcl2_local1 = p6_ptcl2_local1 + bplus_ratio_ptcl2*(p7_ptcl2_local1 - p6_ptcl2_local1);
    Vec p16_ptcl2_local1 = p5_ptcl2_local1 + bplus_ratio_ptcl2*(p8_ptcl2_local1 - p5_ptcl2_local1);

    Vec p17_ptcl2_local1 = p1_ptcl2_local1 + cplus_ratio_ptcl2*(p5_ptcl2_local1 - p1_ptcl2_local1);
    Vec p18_ptcl2_local1 = p2_ptcl2_local1 + cplus_ratio_ptcl2*(p6_ptcl2_local1 - p2_ptcl2_local1);
    Vec p19_ptcl2_local1 = p3_ptcl2_local1 + cplus_ratio_ptcl2*(p7_ptcl2_local1 - p3_ptcl2_local1);
    Vec p20_ptcl2_local1 = p4_ptcl2_local1 + cplus_ratio_ptcl2*(p8_ptcl2_local1 - p4_ptcl2_local1);

    Vec p21_ptcl2_local1 = p13_ptcl2_local1 + aplus_ratio_ptcl2*(p14_ptcl2_local1 - p13_ptcl2_local1);
    Vec p22_ptcl2_local1 = p9_ptcl2_local1  + cplus_ratio_ptcl2*(p12_ptcl2_local1 - p9_ptcl2_local1);
    Vec p23_ptcl2_local1 = p14_ptcl2_local1 + cplus_ratio_ptcl2*(p15_ptcl2_local1 - p14_ptcl2_local1);
    Vec p24_ptcl2_local1 = p10_ptcl2_local1 + cplus_ratio_ptcl2*(p11_ptcl2_local1 - p10_ptcl2_local1);
    Vec p25_ptcl2_local1 = p13_ptcl2_local1 + cplus_ratio_ptcl2*(p16_ptcl2_local1 - p13_ptcl2_local1);
    Vec p26_ptcl2_local1 = p16_ptcl2_local1 + aplus_ratio_ptcl2*(p15_ptcl2_local1 - p16_ptcl2_local1);
    Vec p27_ptcl2_local1 = p21_ptcl2_local1 + cplus_ratio_ptcl2*(p26_ptcl2_local1 - p21_ptcl2_local1);




    int xsign_p1, ysign_p1, zsign_p1, xsign_p2, ysign_p2, zsign_p2;	// octants that p1 and p2 should use

    xsign_p1=0; ysign_p1=0; zsign_p1=0; xsign_p2=0; ysign_p2=0; zsign_p2=0;
    REAL coef_p1[10],coef_p2[10];

    Vec v[2];
    bool b1, b2;
    Vec local1_point1, local2_point2;

    bool flag_contact = false;	// the flag if the two polyellipsoids are contact
    for(int num_ptcl1 = 1; num_ptcl1 !=9; num_ptcl1++){
	for(int num_ptcl2 = 1; num_ptcl2 != 9; num_ptcl2++){

	    // CHECK IF SUB BOX PAIR IS IN CONTACT OR NOT

	    // use sub box of particle 1 as the boundary particle 
	    // and check if sub box of particle 2 is outside of this sub box
	    Vec p1_sub_ptcl2_local1, p2_sub_ptcl2_local1, p3_sub_ptcl2_local1, p4_sub_ptcl2_local1;
	    Vec p5_sub_ptcl2_local1, p6_sub_ptcl2_local1, p7_sub_ptcl2_local1, p8_sub_ptcl2_local1;
	    REAL aplus_sub_ptcl1, aminus_sub_ptcl1;
	    REAL bplus_sub_ptcl1, bminus_sub_ptcl1;
	    REAL cplus_sub_ptcl1, cminus_sub_ptcl1;

	    // use sub box of particle 2 as the boundary particle 
	    // and check if sub box of particle 1 is outside of this sub box
	    Vec p1_sub_ptcl1_local2, p2_sub_ptcl1_local2, p3_sub_ptcl1_local2, p4_sub_ptcl1_local2;
	    Vec p5_sub_ptcl1_local2, p6_sub_ptcl1_local2, p7_sub_ptcl1_local2, p8_sub_ptcl1_local2;
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
    	    if(p1_sub_ptcl2_local1.getX()>aplus_sub_ptcl1 && p2_sub_ptcl2_local1.getX()>aplus_sub_ptcl1
    	    && p3_sub_ptcl2_local1.getX()>aplus_sub_ptcl1 && p4_sub_ptcl2_local1.getX()>aplus_sub_ptcl1
    	    && p5_sub_ptcl2_local1.getX()>aplus_sub_ptcl1 && p6_sub_ptcl2_local1.getX()>aplus_sub_ptcl1
    	    && p7_sub_ptcl2_local1.getX()>aplus_sub_ptcl1 && p8_sub_ptcl2_local1.getX()>aplus_sub_ptcl1 )	// outside of x+ boundary
		continue;	// sub boxes are not contact, go to the next sub box pair

    	    // x- boundary, i.e. a- boundary
    	    if(p1_sub_ptcl2_local1.getX()<-aminus_sub_ptcl1 && p2_sub_ptcl2_local1.getX()<-aminus_sub_ptcl1
    	    && p3_sub_ptcl2_local1.getX()<-aminus_sub_ptcl1 && p4_sub_ptcl2_local1.getX()<-aminus_sub_ptcl1
    	    && p5_sub_ptcl2_local1.getX()<-aminus_sub_ptcl1 && p6_sub_ptcl2_local1.getX()<-aminus_sub_ptcl1
    	    && p7_sub_ptcl2_local1.getX()<-aminus_sub_ptcl1 && p8_sub_ptcl2_local1.getX()<-aminus_sub_ptcl1 )	// outside of x- boundary
		continue;	// sub boxes are not contact, go to the next sub box pair

    	    // y+ boundary, i.e. b+ boundary
    	    if(p1_sub_ptcl2_local1.getY()>bplus_sub_ptcl1 && p2_sub_ptcl2_local1.getY()>bplus_sub_ptcl1
    	    && p3_sub_ptcl2_local1.getY()>bplus_sub_ptcl1 && p4_sub_ptcl2_local1.getY()>bplus_sub_ptcl1
    	    && p5_sub_ptcl2_local1.getY()>bplus_sub_ptcl1 && p6_sub_ptcl2_local1.getY()>bplus_sub_ptcl1
    	    && p7_sub_ptcl2_local1.getY()>bplus_sub_ptcl1 && p8_sub_ptcl2_local1.getY()>bplus_sub_ptcl1 )	// outside of y+ boundary
		continue;	// sub boxes are not contact, go to the next sub box pair

    	    // y- boundary, i.e. b- boundary
    	    if(p1_sub_ptcl2_local1.getY()<-bminus_sub_ptcl1 && p2_sub_ptcl2_local1.getY()<-bminus_sub_ptcl1
    	    && p3_sub_ptcl2_local1.getY()<-bminus_sub_ptcl1 && p4_sub_ptcl2_local1.getY()<-bminus_sub_ptcl1
    	    && p5_sub_ptcl2_local1.getY()<-bminus_sub_ptcl1 && p6_sub_ptcl2_local1.getY()<-bminus_sub_ptcl1
    	    && p7_sub_ptcl2_local1.getY()<-bminus_sub_ptcl1 && p8_sub_ptcl2_local1.getY()<-bminus_sub_ptcl1 )	// outside of y- boundary
		continue;	// sub boxes are not contact, go to the next sub box pair

    	    // z+ boundary, i.e. c+ boundary
    	    if(p1_sub_ptcl2_local1.getZ()>cplus_sub_ptcl1 && p2_sub_ptcl2_local1.getZ()>cplus_sub_ptcl1
    	    && p3_sub_ptcl2_local1.getZ()>cplus_sub_ptcl1 && p4_sub_ptcl2_local1.getZ()>cplus_sub_ptcl1
    	    && p5_sub_ptcl2_local1.getZ()>cplus_sub_ptcl1 && p6_sub_ptcl2_local1.getZ()>cplus_sub_ptcl1
    	    && p7_sub_ptcl2_local1.getZ()>cplus_sub_ptcl1 && p8_sub_ptcl2_local1.getZ()>cplus_sub_ptcl1 )	// outside of z+ boundary
		continue;	// sub boxes are not contact, go to the next sub box pair

    	    // z- boundary, i.e. c- boundary
    	    if(p1_sub_ptcl2_local1.getZ()<-cminus_sub_ptcl1 && p2_sub_ptcl2_local1.getZ()<-cminus_sub_ptcl1
    	    && p3_sub_ptcl2_local1.getZ()<-cminus_sub_ptcl1 && p4_sub_ptcl2_local1.getZ()<-cminus_sub_ptcl1
    	    && p5_sub_ptcl2_local1.getZ()<-cminus_sub_ptcl1 && p6_sub_ptcl2_local1.getZ()<-cminus_sub_ptcl1
    	    && p7_sub_ptcl2_local1.getZ()<-cminus_sub_ptcl1 && p8_sub_ptcl2_local1.getZ()<-cminus_sub_ptcl1 )	// outside of z- boundary
		continue;	// sub boxes are not contact, go to the next sub box pair


	    //////////////////////////////////////////////////////////////////////////////////////////////////////
    	    // check if all 8 corner points of the sub box of particle 1 are outside of the sub box of particle 2
    	    // x+ boundary, i.e. a+ boundary
    	    if(p1_sub_ptcl1_local2.getX()>aplus_sub_ptcl2 && p2_sub_ptcl1_local2.getX()>aplus_sub_ptcl2
    	    && p3_sub_ptcl1_local2.getX()>aplus_sub_ptcl2 && p4_sub_ptcl1_local2.getX()>aplus_sub_ptcl2
    	    && p5_sub_ptcl1_local2.getX()>aplus_sub_ptcl2 && p6_sub_ptcl1_local2.getX()>aplus_sub_ptcl2
    	    && p7_sub_ptcl1_local2.getX()>aplus_sub_ptcl2 && p8_sub_ptcl1_local2.getX()>aplus_sub_ptcl2 )	// outside of x+ boundary
		continue;	// sub boxes are not contact, go to the next sub box pair

    	    // x- boundary, i.e. a- boundary
    	    if(p1_sub_ptcl1_local2.getX()<-aminus_sub_ptcl2 && p2_sub_ptcl1_local2.getX()<-aminus_sub_ptcl2
    	    && p3_sub_ptcl1_local2.getX()<-aminus_sub_ptcl2 && p4_sub_ptcl1_local2.getX()<-aminus_sub_ptcl2
    	    && p5_sub_ptcl1_local2.getX()<-aminus_sub_ptcl2 && p6_sub_ptcl1_local2.getX()<-aminus_sub_ptcl2
    	    && p7_sub_ptcl1_local2.getX()<-aminus_sub_ptcl2 && p8_sub_ptcl1_local2.getX()<-aminus_sub_ptcl2 )	// outside of x- boundary
		continue;	// sub boxes are not contact, go to the next sub box pair

    	    // y+ boundary, i.e. b+ boundary
    	    if(p1_sub_ptcl1_local2.getY()>bplus_sub_ptcl2 && p2_sub_ptcl1_local2.getY()>bplus_sub_ptcl2
    	    && p3_sub_ptcl1_local2.getY()>bplus_sub_ptcl2 && p4_sub_ptcl1_local2.getY()>bplus_sub_ptcl2
    	    && p5_sub_ptcl1_local2.getY()>bplus_sub_ptcl2 && p6_sub_ptcl1_local2.getY()>bplus_sub_ptcl2
    	    && p7_sub_ptcl1_local2.getY()>bplus_sub_ptcl2 && p8_sub_ptcl1_local2.getY()>bplus_sub_ptcl2 )	// outside of y+ boundary
		continue;	// sub boxes are not contact, go to the next sub box pair

    	    // y- boundary, i.e. b- boundary
    	    if(p1_sub_ptcl1_local2.getY()<-bminus_sub_ptcl2 && p2_sub_ptcl1_local2.getY()<-bminus_sub_ptcl2
    	    && p3_sub_ptcl1_local2.getY()<-bminus_sub_ptcl2 && p4_sub_ptcl1_local2.getY()<-bminus_sub_ptcl2
    	    && p5_sub_ptcl1_local2.getY()<-bminus_sub_ptcl2 && p6_sub_ptcl1_local2.getY()<-bminus_sub_ptcl2
    	    && p7_sub_ptcl1_local2.getY()<-bminus_sub_ptcl2 && p8_sub_ptcl1_local2.getY()<-bminus_sub_ptcl2 )	// outside of y- boundary
		continue;	// sub boxes are not contact, go to the next sub box pair

    	    // z+ boundary, i.e. c+ boundary
    	    if(p1_sub_ptcl1_local2.getZ()>cplus_sub_ptcl2 && p2_sub_ptcl1_local2.getZ()>cplus_sub_ptcl2
    	    && p3_sub_ptcl1_local2.getZ()>cplus_sub_ptcl2 && p4_sub_ptcl1_local2.getZ()>cplus_sub_ptcl2
    	    && p5_sub_ptcl1_local2.getZ()>cplus_sub_ptcl2 && p6_sub_ptcl1_local2.getZ()>cplus_sub_ptcl2
    	    && p7_sub_ptcl1_local2.getZ()>cplus_sub_ptcl2 && p8_sub_ptcl1_local2.getZ()>cplus_sub_ptcl2 )	// outside of z+ boundary
		continue;	// sub boxes are not contact, go to the next sub box pair

    	    // z- boundary, i.e. c- boundary
    	    if(p1_sub_ptcl1_local2.getZ()<-cminus_sub_ptcl2 && p2_sub_ptcl1_local2.getZ()<-cminus_sub_ptcl2
    	    && p3_sub_ptcl1_local2.getZ()<-cminus_sub_ptcl2 && p4_sub_ptcl1_local2.getZ()<-cminus_sub_ptcl2
    	    && p5_sub_ptcl1_local2.getZ()<-cminus_sub_ptcl2 && p6_sub_ptcl1_local2.getZ()<-cminus_sub_ptcl2
    	    && p7_sub_ptcl1_local2.getZ()<-cminus_sub_ptcl2 && p8_sub_ptcl1_local2.getZ()<-cminus_sub_ptcl2 )	// outside of z- boundary
		continue;	// sub boxes are not contact, go to the next sub box pair



	    // USE THE CORRESPONDING OCTANTS TO CALCULATE THE CONTACTS OF ELLIPSOIDS

    	    p1->getGlobalCoef(coef_p1, num_ptcl1); // v[0] is the point on p2, v[1] is the point on p1
    	    p2->getGlobalCoef(coef_p2, num_ptcl2);

    	    b1 = root6(coef_p1,coef_p2,v[0]);
    	    b2 = root6(coef_p2,coef_p1,v[1]);
    	    point1 = v[1];
    	    point2 = v[0];

    	    local1_point1 = p1->globalToLocal(point1 - p1->getCurrPos());
    	    local2_point2 = p2->globalToLocal(point2 - p2->getCurrPos());

    	    if(b1 && b2 
    	    && local1_point1.getX()*xsign_p1 >= 0 && local1_point1.getY()*ysign_p1 >= 0 && local1_point1.getZ()*zsign_p1 >= 0
    	    && local2_point2.getX()*xsign_p2 >= 0 && local2_point2.getY()*ysign_p2 >= 0 && local2_point2.getZ()*zsign_p2 >= 0 ){
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




//    REAL coef1[10],coef2[10];
//    p1->getGlobalCoef(coef1); // v[0] is the point on p2, v[1] is the point on p1
//    p2->getGlobalCoef(coef2);    
//    Vec v[2];
//    bool b1 = root6(coef1,coef2,v[0]);
//    bool b2 = root6(coef2,coef1,v[1]);
//    point1 = v[1];
//    point2 = v[0];
    radius1=p1->getRadius(point1, xsign_p1, ysign_p1, zsign_p1);
    radius2=p2->getRadius(point2, xsign_p2, ysign_p2, zsign_p2);
    penetr = vfabs(point1-point2);

    if (b1 && b2 
	&& penetr/(2.0*fmax(radius1,radius2)) > dem::Parameter::getSingleton().parameter["minRelaOverlap"]
	&& nearbyint(penetr/dem::Parameter::getSingleton().parameter["measureOverlap"]) >= 1) { // a strict detection method
        isInContact = true;
        return true;
    }
    else {
        isInContact = false;
	return false;
    }
  }
  
  
  void Contact::checkinPrevTgt(std::vector<ContactTgt>& contactTgtVec) {
    for(std::vector<ContactTgt>::iterator it = contactTgtVec.begin();it != contactTgtVec.end(); ++it) {
      if (it->ptcl1 == p1->getId() && it->ptcl2 == p2->getId()) {
	prevTgtForce   = it->tgtForce;
	prevTgtDisp    = it->tgtDisp;
	prevTgtLoading = it->tgtLoading;
	tgtDispStart   = it->tgtDispStart;
	tgtPeak        = it->tgtPeak;
	prevTgtSlide   = it->tgtSlide;
	break;
      }
    }
  }
  
  
  void Contact::checkoutTgt(std::vector<ContactTgt>& contactTgtVec) {
    contactTgtVec.push_back(ContactTgt(p1->getId(),
				       p2->getId(),
				       tgtForce,
				       tgtDisp, 
				       tgtLoading,
				       tgtDispStart,
				       tgtPeak,
				       tgtSlide) );
  }  
  
  
  void Contact::contactForce() {
    // isOverlapped() has been called in findContact() in assembly.cpp and information recorded, 
    // now this function is called by internalForce() in assembly.cpp.
    
    if (isInContact) {

      REAL young = dem::Parameter::getSingleton().parameter["young"];
      REAL poisson = dem::Parameter::getSingleton().parameter["poisson"];
      REAL maxRelaOverlap = dem::Parameter::getSingleton().parameter["maxRelaOverlap"];
      REAL measureOverlap = dem::Parameter::getSingleton().parameter["measureOverlap"];
      REAL contactCohesion = dem::Parameter::getSingleton().parameter["contactCohesion"];
      REAL contactDamp = dem::Parameter::getSingleton().parameter["contactDamp"];
      REAL contactFric = dem::Parameter::getSingleton().parameter["contactFric"];

      // obtain normal force, using absolute equation instead of stiffness method
      p1->setContactNum(p1->getContactNum() + 1);
      p2->setContactNum(p2->getContactNum() + 1);
      p1->setInContact(true);
      p2->setInContact(true);
      
      R0 = radius1*radius2/(radius1+radius2);
      E0 = 0.5*young/(1-poisson*poisson);
      REAL allowedOverlap = 2.0 * fmin(radius1,radius2) * maxRelaOverlap;
      if (penetr > allowedOverlap) {
#ifndef NDEBUG
	std::stringstream inf;
	inf.setf(std::ios::scientific, std::ios::floatfield);
	inf << " Contact.cpp: iter=" << std::setw(8) << iteration 
	    << " ptcl1=" << std::setw(8) << getP1()->getId()
	    << " ptcl2=" << std::setw(8) << getP2()->getId()
	    << " penetr="<< std::setw(OWID) << penetr 
	    << " allow=" << std::setw(OWID) << allowedOverlap
	    << std::endl;
	MPI_Status status;
	int length = OWID*2 + 8*3 + 19 + 7*3 + 8 + 1;
	MPI_File_write_shared(overlapInf, const_cast<char*> (inf.str().c_str()), length, MPI_CHAR, &status);
#endif
	penetr = allowedOverlap;
      }

      penetr = nearbyint (penetr/measureOverlap) * measureOverlap;
      contactRadius = sqrt(penetr*R0);
      normalDirc = normalize(point1-point2); // normalDirc points from particle 1 to particle 2
      normalForce = -sqrt(penetr*penetr*penetr)*sqrt(R0)*4*E0/3.0* normalDirc; // normalForce pointing to particle 1
      // pow(penetr, 1.5)
      
      // apply cohesion force
      cohesionForce = Pi * (penetr*R0) * contactCohesion * normalDirc;
      p1->addForce(cohesionForce);
      p2->addForce(-cohesionForce);
      
      // apply normal force
      p1->addForce(normalForce);
      p2->addForce(-normalForce);
      p1->addMoment( ( (point1+point2)*0.5-p1->getCurrCenterMass() ) % normalForce );
      p2->addMoment( ( (point1+point2)*0.5-p2->getCurrCenterMass() ) % (-normalForce) );	
      
      /*
      debugInf << "Contact.h: iter=" << iteration
	       << " penetr=" << penetr
	       << " cohesionForce=" << vfabs(cohesionForce)
	       << " normalForce=" << vfabs(normalForce)
	       << " accumulated time=" << iteration * timeStep
	       << std::endl;
      */
      
      // obtain normal damping force
      Vec cp = (point1+point2)*0.5;        
      Vec veloc1 = p1->getCurrVeloc() + p1->getCurrOmga() % (cp-p1->getCurrCenterMass());
      Vec veloc2 = p2->getCurrVeloc() + p2->getCurrOmga() % (cp-p2->getCurrCenterMass());
      REAL m1 = getP1()->getMass();
      REAL m2 = getP2()->getMass();
      REAL kn = pow(6*vfabs(normalForce)*R0*pow(E0,2),1.0/3.0);
      REAL dampCritical = 2*sqrt(m1*m2/(m1+m2)*kn); // critical damping
      Vec cntDampingForce = contactDamp * dampCritical * ((veloc1-veloc2) * normalDirc)*normalDirc;
      vibraTimeStep = 2.0*sqrt( m1*m2 / (m1+m2) /kn );
      impactTimeStep = allowedOverlap / fabs((veloc1-veloc2) * normalDirc);

      // apply normal damping force
      p1->addForce(-cntDampingForce);
      p2->addForce(cntDampingForce);
      p1->addMoment( ( (point1+point2)*0.5-p1->getCurrCenterMass() ) % (-cntDampingForce) );
      p2->addMoment( ( (point1+point2)*0.5-p2->getCurrCenterMass() ) % cntDampingForce );
      
      if (contactFric != 0) {
	// obtain tangential force
	G0 = young*0.5/(1+poisson);              // RelaDispInc points along point1's displacement relative to point2
	Vec RelaDispInc = (veloc1-veloc2) * timeStep;
	Vec tgtDispInc = RelaDispInc - (RelaDispInc * normalDirc)*normalDirc;
	tgtDisp = prevTgtDisp + tgtDispInc; // prevTgtDisp read by checkinPrevTgt()
	if (vfabs(tgtDisp) == 0)
	  tgtDirc = 0;
	else
	  tgtDirc = normalize(-tgtDisp); // tgtDirc points along Tgtential forces exerted on particle 1
	
	REAL fP = 0;
	REAL ks = 0;
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// linear friction model
	fP = contactFric*vfabs(normalForce);
	ks = 4*G0*contactRadius/(2-poisson);
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
	fP = contactFric*vfabs(normalForce);
	tgtLoading = (prevTgtDisp * tgtDispInc >= 0); 
	
	if (tgtLoading) {              // loading
	  if (!prevTgtLoading) {      // pre-step is unloading
	    val = 8*G0*contactRadius*vfabs(tgtDispInc)/(3*(2-poisson)*fP);
	    tgtDispStart = prevTgtDisp;
	  }
	  else                       // pre-step is loading
	    val = 8*G0*contactRadius*vfabs(tgtDisp-tgtDispStart)/(3*(2-poisson)*fP);
	  
	  if (val > 1.0)              
	    tgtForce = fP*tgtDirc;
	  else {
	    ks = 4*G0*contactRadius/(2-poisson)*sqrt(1-val);
	    //incremental method
	    tgtForce = prevTgtForce + ks*(-tgtDispInc); // tgtDispInc determines signs
	    //total value method: tgtForce = fP*(1-pow(1-val, 1.5))*tgtDirc;
	  }
	}
	else {                         // unloading
	  if (prevTgtLoading) {       // pre-step is loading
	    val = 8*G0*contactRadius*vfabs(tgtDisp-tgtDispStart)/(3*(2-poisson)*fP);
	    tgtPeak = vfabs(prevTgtForce);
	  }
	  else                       // pre-step is unloading
	    val = 8*G0*contactRadius*vfabs(tgtDisp-tgtDispStart)/(3*(2-poisson)*fP);
	  
	  if (val > 1.0 || tgtPeak > fP)  
	    tgtForce = fP*tgtDirc;
	  else {
	    ks = 2*sqrt(2)*G0*contactRadius/(2-poisson) * sqrt(1+pow(1-tgtPeak/fP,2.0/3.0)+val);
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
	fP = contactFric*vfabs(normalForce);
	if (prevTgtSlide)
	  val = 8*G0*contactRadius*vfabs(tgtDispInc)/(3*(2-poisson)*fP);
	else
	  val = 8*G0*contactRadius*vfabs(tgtDisp-tgtDispStart)/(3*(2-poisson)*fP);
	
	if (iteration > 10000 && iteration < 11000) { // loading (and possible sliding)
	  if (val > 1.0) {
	    tgtForce = fP*tgtDirc;
	    tgtSlide = true;
	  }
	  else {
	    if (!prevTgtSlide) {
	      ks = 4*G0*contactRadius/(2-poisson)*sqrt(1-val);
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
	      ks = 2*sqrt(2)*G0*contactRadius/(2-poisson) * sqrt(1+pow(1-tgtPeak/fP,2.0/3.0)+val);
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
		 << " tgtDispInc.x=" << tgtDispInc.getX()
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
	p1->addMoment( ((point1+point2)*0.5-p1->getCurrCenterMass()) % tgtForce);
	p2->addMoment(-((point1+point2)*0.5-p2->getCurrCenterMass()) % tgtForce);
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
      spinResist = 0;
    }   
    
  }
  
} // namespace dem ends
