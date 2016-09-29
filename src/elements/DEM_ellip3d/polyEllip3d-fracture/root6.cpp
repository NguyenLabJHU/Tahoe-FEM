//    <<Numerical Recipes in C>>
//    void zrhqr(float a[], int m, float rtr[], float rti[])
//    Find all the roots of a polynomial with real coefficients, i=0~m, a(i)x^i, given the degree m
//    and the coefficients a[0..m]. The method is to construct an upper Hessenberg matrix whose
//    eigenvalues are the desired roots, and then use the routines balanc and hqr. The real and
//    imaginary parts of the roots are returned in rtr[1..m] and rti[1..m], respectively

//    The purpose of this function is to find the point on surface of particle 2, which penetrates
//    deepest into surface of particle 1. However, the algorithm in the thesis can't guarantee that
//    the point found is inside particle 1, which means, even if the two particles are completely 
//    separated (not-in-touch), the algorithm could still find a point!
//
//    Input: coef1[0] and coef2[0] are always 1.0, which has been guaranteed by the caller. 
//
//    Return values:
//      false - non-overlapped
//      true  - overlapped and vector point returns the deepest penetrated point.

#include "root6.h"
#include "parameter.h"
#include "nr.h"
#include <iostream>
#include <fstream>
#include <iomanip>

namespace dem {

bool root6(REAL coef1[],REAL coef2[],vec& point){
	if((coef1[0]==1.0&&coef1[1]==1.0&&coef1[2]==1.0&&
		coef1[3]==0.0&&coef1[4]==0.0&&coef1[5]==0.0)&&
		(coef2[0]==1.0&&coef2[1]==1.0&&coef2[2]==1.0&&
		coef2[3]==0.0&&coef2[4]==0.0&&coef2[5]==0.0)){
		REAL x1=-coef1[6]/2;
		REAL y1=-coef1[7]/2;
		REAL z1=-coef1[8]/2;
		REAL R1=sqrt(x1*x1+y1*y1+z1*z1-coef1[9]);
		REAL x2=-coef2[6]/2;
		REAL y2=-coef2[7]/2;
		REAL z2=-coef2[8]/2;
		REAL R2=sqrt(x2*x2+y2*y2+z2*z2-coef2[9]);
		vec  dirc=vec(x1-x2,y1-y2,z1-z2);
		REAL dist=vfabs(dirc);
		dirc=normalize(dirc);
		point=vec(x2,y2,z2)+R2*dirc;
		vec judge=point-vec(x1,y1,z1);
		if(dist < R1 + R2  &&
		   (R1-vfabs(judge)) / (2.0*fmax(R1,R2)) > MINOVERLAP ) // satisfy minimum relative overlap
		    return true;  // overlapped
		else
		    return false; // non-overlapped
	}

	REAL a1=coef1[0];
	REAL b1=coef1[1];
	REAL c1=coef1[2];
	REAL d1=coef1[3];
	REAL e1=coef1[4];
	REAL f1=coef1[5];
	REAL g1=coef1[6];
	REAL h1=coef1[7];
	REAL i1=coef1[8];
	REAL j1=coef1[9];

	REAL a2=coef2[0];
	REAL b2=coef2[1];
	REAL c2=coef2[2];
	REAL d2=coef2[3];
	REAL e2=coef2[4];
	REAL f2=coef2[5];
	REAL g2=coef2[6];
	REAL h2=coef2[7];
	REAL i2=coef2[8];
	REAL j2=coef2[9];
	REAL rtc[7], rti[7], rtr[7];
	rtc[0]=
		-8*b1*c1*d1*e1*f1*g1*g2 + 8*a1*c1*d1*e1*e2*g1*h1 + 
	   8*a2*b1*c1*e1*f1*g1*h1 + 8*a1*b2*c1*e1*f1*g1*h1 + 
	   8*a1*b1*c2*e1*f1*g1*h1 - 4*c1*d1*d2*e1*f1*g1*h1 - 
	   8*a1*b1*c1*e2*f1*g1*h1 - 8*a1*b1*c1*e1*f2*g1*h1 + 
	   8*b1*c1*d1*f1*f2*g1*h1 - 8*a1*b1*c1*e1*f1*g2*h1 - 
	   8*a1*b1*c1*e1*f1*g1*h2 - 8*a1*c1*d1*e1*f1*h1*h2 + 
	   8*a2*b1*c1*d1*e1*g1*i1 + 8*a1*b2*c1*d1*e1*g1*i1 + 
	   8*a1*b1*c2*d1*e1*g1*i1 - 8*a1*b1*c1*d2*e1*g1*i1 - 
	   8*a1*b1*c1*d1*e2*g1*i1 + 8*b1*c1*d1*d2*f1*g1*i1 + 
	   8*a1*b1*e1*e2*f1*g1*i1 - 4*b1*d1*e1*f1*f2*g1*i1 - 
	   8*a1*b1*c1*d1*e1*g2*i1 + 8*a1*c1*d1*d2*e1*h1*i1 + 
	   8*a2*b1*c1*d1*f1*h1*i1 + 8*a1*b2*c1*d1*f1*h1*i1 + 
	   8*a1*b1*c2*d1*f1*h1*i1 - 8*a1*b1*c1*d2*f1*h1*i1 - 
	   4*a1*d1*e1*e2*f1*h1*i1 - 8*a1*b1*c1*d1*f2*h1*i1 + 
	   8*a1*b1*e1*f1*f2*h1*i1 - 8*a1*b1*c1*d1*f1*h2*i1 - 
	   8*a1*b1*c1*d1*e1*g1*i2 - 8*a1*b1*c1*d1*f1*h1*i2 - 
	   8*a1*b1*d1*e1*f1*i1*i2 + 32*a1*b1*c1*d1*e1*f1*j2 - 
	   16*b2*c1*e1*h1*i1*pow(a1,2) - 
	   16*b1*c2*e1*h1*i1*pow(a1,2) + 
	   16*b1*c1*e2*h1*i1*pow(a1,2) + 
	   16*b1*c1*e1*h2*i1*pow(a1,2) + 
	   16*b1*c1*e1*h1*i2*pow(a1,2) - 
	   16*a2*c1*f1*g1*i1*pow(b1,2) - 
	   16*a1*c2*f1*g1*i1*pow(b1,2) + 
	   16*a1*c1*f2*g1*i1*pow(b1,2) + 
	   16*a1*c1*f1*g2*i1*pow(b1,2) + 
	   16*a1*c1*f1*g1*i2*pow(b1,2) - 
	   32*c1*i1*i2*pow(a1,2)*pow(b1,2) - 
	   16*a2*b1*d1*g1*h1*pow(c1,2) - 
	   16*a1*b2*d1*g1*h1*pow(c1,2) + 
	   16*a1*b1*d2*g1*h1*pow(c1,2) + 
	   16*a1*b1*d1*g2*h1*pow(c1,2) + 
	   16*a1*b1*d1*g1*h2*pow(c1,2) - 
	   32*b1*h1*h2*pow(a1,2)*pow(c1,2) - 
	   32*a1*g1*g2*pow(b1,2)*pow(c1,2) + 
	   64*j2*pow(a1,2)*pow(b1,2)*pow(c1,2) + 
	   2*c2*e1*f1*g1*h1*pow(d1,2) - 
	   2*c1*e2*f1*g1*h1*pow(d1,2) - 
	   2*c1*e1*f2*g1*h1*pow(d1,2) + 
	   6*c1*e1*f1*g2*h1*pow(d1,2) + 
	   6*c1*e1*f1*g1*h2*pow(d1,2) - 
	   2*c1*d2*e1*g1*i1*pow(d1,2) - 
	   4*b2*c1*f1*g1*i1*pow(d1,2) + 
	   4*b1*c2*f1*g1*i1*pow(d1,2) - 
	   4*b1*c1*f2*g1*i1*pow(d1,2) - 
	   4*b1*c1*f1*g2*i1*pow(d1,2) - 
	   4*a2*c1*e1*h1*i1*pow(d1,2) + 
	   4*a1*c2*e1*h1*i1*pow(d1,2) - 
	   4*a1*c1*e2*h1*i1*pow(d1,2) - 
	   2*c1*d2*f1*h1*i1*pow(d1,2) - 
	   4*a1*c1*e1*h2*i1*pow(d1,2) - 
	   4*b1*c1*f1*g1*i2*pow(d1,2) - 
	   4*a1*c1*e1*h1*i2*pow(d1,2) + 
	   16*a1*b1*c1*i1*i2*pow(d1,2) + 
	   8*b1*g1*g2*pow(c1,2)*pow(d1,2) + 
	   4*d2*g1*h1*pow(c1,2)*pow(d1,2) + 
	   8*a1*h1*h2*pow(c1,2)*pow(d1,2) - 
	   32*a1*b1*j2*pow(c1,2)*pow(d1,2) - 
	   2*c2*e1*g1*i1*pow(d1,3) + 2*c1*e2*g1*i1*pow(d1,3) + 
	   2*c1*e1*g2*i1*pow(d1,3) - 2*c2*f1*h1*i1*pow(d1,3) + 
	   2*c1*f2*h1*i1*pow(d1,3) + 2*c1*f1*h2*i1*pow(d1,3) + 
	   2*c1*e1*g1*i2*pow(d1,3) + 2*c1*f1*h1*i2*pow(d1,3) + 
	   2*e1*f1*i1*i2*pow(d1,3) - 8*c1*e1*f1*j2*pow(d1,3) - 
	   4*g2*h1*pow(c1,2)*pow(d1,3) - 
	   4*g1*h2*pow(c1,2)*pow(d1,3) - 2*c1*i1*i2*pow(d1,4) + 
	   4*j2*pow(c1,2)*pow(d1,4) + 
	   16*a1*b1*c1*g1*g2*pow(e1,2) + 
	   4*a2*c1*d1*g1*h1*pow(e1,2) - 
	   4*a1*c2*d1*g1*h1*pow(e1,2) - 
	   4*a1*c1*d2*g1*h1*pow(e1,2) - 
	   2*a1*e2*f1*g1*h1*pow(e1,2) - 
	   4*a1*c1*d1*g2*h1*pow(e1,2) - 
	   4*a1*c1*d1*g1*h2*pow(e1,2) - 
	   2*a1*d1*e2*g1*i1*pow(e1,2) + 
	   4*a2*b1*f1*g1*i1*pow(e1,2) - 
	   4*a1*b2*f1*g1*i1*pow(e1,2) - 
	   4*a1*b1*f2*g1*i1*pow(e1,2) - 
	   4*a1*b1*f1*g2*i1*pow(e1,2) + 
	   2*a2*d1*f1*h1*i1*pow(e1,2) - 
	   2*a1*d2*f1*h1*i1*pow(e1,2) - 
	   2*a1*d1*f2*h1*i1*pow(e1,2) + 
	   6*a1*d1*f1*h2*i1*pow(e1,2) - 
	   4*a1*b1*f1*g1*i2*pow(e1,2) + 
	   6*a1*d1*f1*h1*i2*pow(e1,2) + 
	   8*c1*h1*h2*pow(a1,2)*pow(e1,2) + 
	   4*e2*h1*i1*pow(a1,2)*pow(e1,2) + 
	   8*b1*i1*i2*pow(a1,2)*pow(e1,2) - 
	   32*b1*c1*j2*pow(a1,2)*pow(e1,2) - 
	   2*c1*g1*g2*pow(d1,2)*pow(e1,2) + 
	   2*f2*g1*i1*pow(d1,2)*pow(e1,2) - 
	   2*f1*g2*i1*pow(d1,2)*pow(e1,2) - 
	   2*f1*g1*i2*pow(d1,2)*pow(e1,2) - 
	   2*a1*i1*i2*pow(d1,2)*pow(e1,2) + 
	   8*a1*c1*j2*pow(d1,2)*pow(e1,2) + 
	   2*d1*f1*g1*g2*pow(e1,3) - 2*a2*f1*g1*h1*pow(e1,3) + 
	   2*a1*f2*g1*h1*pow(e1,3) + 2*a1*f1*g2*h1*pow(e1,3) + 
	   2*a1*f1*g1*h2*pow(e1,3) - 2*a2*d1*g1*i1*pow(e1,3) + 
	   2*a1*d2*g1*i1*pow(e1,3) + 2*a1*d1*g2*i1*pow(e1,3) + 
	   2*a1*d1*g1*i2*pow(e1,3) - 8*a1*d1*f1*j2*pow(e1,3) - 
	   4*h2*i1*pow(a1,2)*pow(e1,3) - 
	   4*h1*i2*pow(a1,2)*pow(e1,3) - 2*a1*g1*g2*pow(e1,4) + 
	   4*j2*pow(a1,2)*pow(e1,4) + 4*b2*c1*d1*g1*h1*pow(f1,2) - 
	   4*b1*c2*d1*g1*h1*pow(f1,2) - 
	   4*b1*c1*d2*g1*h1*pow(f1,2) - 
	   2*b1*e1*f2*g1*h1*pow(f1,2) - 
	   4*b1*c1*d1*g2*h1*pow(f1,2) - 
	   4*b1*c1*d1*g1*h2*pow(f1,2) + 
	   16*a1*b1*c1*h1*h2*pow(f1,2) + 
	   2*b2*d1*e1*g1*i1*pow(f1,2) - 
	   2*b1*d2*e1*g1*i1*pow(f1,2) - 
	   2*b1*d1*e2*g1*i1*pow(f1,2) + 
	   6*b1*d1*e1*g2*i1*pow(f1,2) - 
	   4*a2*b1*e1*h1*i1*pow(f1,2) + 
	   4*a1*b2*e1*h1*i1*pow(f1,2) - 
	   4*a1*b1*e2*h1*i1*pow(f1,2) - 
	   2*b1*d1*f2*h1*i1*pow(f1,2) - 
	   4*a1*b1*e1*h2*i1*pow(f1,2) + 
	   6*b1*d1*e1*g1*i2*pow(f1,2) - 
	   4*a1*b1*e1*h1*i2*pow(f1,2) + 
	   8*c1*g1*g2*pow(b1,2)*pow(f1,2) + 
	   4*f2*g1*i1*pow(b1,2)*pow(f1,2) + 
	   8*a1*i1*i2*pow(b1,2)*pow(f1,2) - 
	   32*a1*c1*j2*pow(b1,2)*pow(f1,2) - 
	   2*c1*h1*h2*pow(d1,2)*pow(f1,2) + 
	   2*e2*h1*i1*pow(d1,2)*pow(f1,2) - 
	   2*e1*h2*i1*pow(d1,2)*pow(f1,2) - 
	   2*e1*h1*i2*pow(d1,2)*pow(f1,2) - 
	   2*b1*i1*i2*pow(d1,2)*pow(f1,2) + 
	   8*b1*c1*j2*pow(d1,2)*pow(f1,2) - 
	   2*b1*g1*g2*pow(e1,2)*pow(f1,2) + 
	   2*d2*g1*h1*pow(e1,2)*pow(f1,2) - 
	   2*d1*g2*h1*pow(e1,2)*pow(f1,2) - 
	   2*d1*g1*h2*pow(e1,2)*pow(f1,2) - 
	   2*a1*h1*h2*pow(e1,2)*pow(f1,2) + 
	   8*a1*b1*j2*pow(e1,2)*pow(f1,2) + 
	   4*j2*pow(d1,2)*pow(e1,2)*pow(f1,2) - 
	   2*b2*e1*g1*h1*pow(f1,3) + 2*b1*e2*g1*h1*pow(f1,3) + 
	   2*b1*e1*g2*h1*pow(f1,3) + 2*b1*e1*g1*h2*pow(f1,3) + 
	   2*d1*e1*h1*h2*pow(f1,3) - 2*b2*d1*h1*i1*pow(f1,3) + 
	   2*b1*d2*h1*i1*pow(f1,3) + 2*b1*d1*h2*i1*pow(f1,3) + 
	   2*b1*d1*h1*i2*pow(f1,3) - 8*b1*d1*e1*j2*pow(f1,3) - 
	   4*g2*i1*pow(b1,2)*pow(f1,3) - 
	   4*g1*i2*pow(b1,2)*pow(f1,3) - 2*b1*h1*h2*pow(f1,4) + 
	   4*j2*pow(b1,2)*pow(f1,4) - 4*b2*c1*d1*e1*f1*pow(g1,2) - 
	   4*b1*c2*d1*e1*f1*pow(g1,2) + 
	   4*b1*c1*d2*e1*f1*pow(g1,2) + 
	   4*b1*c1*d1*e2*f1*pow(g1,2) + 
	   4*b1*c1*d1*e1*f2*pow(g1,2) - 
	   8*c1*f1*f2*pow(b1,2)*pow(g1,2) - 
	   8*b1*d1*d2*pow(c1,2)*pow(g1,2) + 
	   16*a2*pow(b1,2)*pow(c1,2)*pow(g1,2) - 
	   2*c1*e1*e2*pow(d1,2)*pow(g1,2) + 
	   4*b2*pow(c1,2)*pow(d1,2)*pow(g1,2) - 
	   8*a2*b1*c1*pow(e1,2)*pow(g1,2) + 
	   2*c1*d1*d2*pow(e1,2)*pow(g1,2) + 
	   d1*e2*f1*pow(e1,2)*pow(g1,2) + 
	   2*b1*f1*f2*pow(e1,2)*pow(g1,2) + 
	   c2*pow(d1,2)*pow(e1,2)*pow(g1,2) - 
	   d2*f1*pow(e1,3)*pow(g1,2) - d1*f2*pow(e1,3)*pow(g1,2) + 
	   a2*pow(e1,4)*pow(g1,2) - 
	   2*b1*e1*e2*pow(f1,2)*pow(g1,2) + 
	   4*c2*pow(b1,2)*pow(f1,2)*pow(g1,2) + 
	   b2*pow(e1,2)*pow(f1,2)*pow(g1,2) - 
	   4*a2*c1*d1*e1*f1*pow(h1,2) - 
	   4*a1*c2*d1*e1*f1*pow(h1,2) + 
	   4*a1*c1*d2*e1*f1*pow(h1,2) + 
	   4*a1*c1*d1*e2*f1*pow(h1,2) + 
	   4*a1*c1*d1*e1*f2*pow(h1,2) - 
	   8*c1*e1*e2*pow(a1,2)*pow(h1,2) - 
	   8*a1*d1*d2*pow(c1,2)*pow(h1,2) + 
	   16*b2*pow(a1,2)*pow(c1,2)*pow(h1,2) - 
	   2*c1*f1*f2*pow(d1,2)*pow(h1,2) + 
	   4*a2*pow(c1,2)*pow(d1,2)*pow(h1,2) - 
	   2*a1*f1*f2*pow(e1,2)*pow(h1,2) + 
	   4*c2*pow(a1,2)*pow(e1,2)*pow(h1,2) - 
	   8*a1*b2*c1*pow(f1,2)*pow(h1,2) + 
	   2*c1*d1*d2*pow(f1,2)*pow(h1,2) + 
	   2*a1*e1*e2*pow(f1,2)*pow(h1,2) + 
	   d1*e1*f2*pow(f1,2)*pow(h1,2) + 
	   c2*pow(d1,2)*pow(f1,2)*pow(h1,2) + 
	   a2*pow(e1,2)*pow(f1,2)*pow(h1,2) - 
	   d2*e1*pow(f1,3)*pow(h1,2) - d1*e2*pow(f1,3)*pow(h1,2) + 
	   b2*pow(f1,4)*pow(h1,2) - 4*a2*b1*d1*e1*f1*pow(i1,2) - 
	   4*a1*b2*d1*e1*f1*pow(i1,2) + 
	   4*a1*b1*d2*e1*f1*pow(i1,2) + 
	   4*a1*b1*d1*e2*f1*pow(i1,2) + 
	   4*a1*b1*d1*e1*f2*pow(i1,2) - 
	   8*b1*e1*e2*pow(a1,2)*pow(i1,2) - 
	   8*a1*f1*f2*pow(b1,2)*pow(i1,2) + 
	   16*c2*pow(a1,2)*pow(b1,2)*pow(i1,2) - 
	   8*a1*b1*c2*pow(d1,2)*pow(i1,2) + 
	   2*a1*e1*e2*pow(d1,2)*pow(i1,2) + 
	   d2*e1*f1*pow(d1,2)*pow(i1,2) + 
	   2*b1*f1*f2*pow(d1,2)*pow(i1,2) - 
	   e2*f1*pow(d1,3)*pow(i1,2) - e1*f2*pow(d1,3)*pow(i1,2) + 
	   c2*pow(d1,4)*pow(i1,2) - 
	   2*a1*d1*d2*pow(e1,2)*pow(i1,2) + 
	   4*b2*pow(a1,2)*pow(e1,2)*pow(i1,2) + 
	   a2*pow(d1,2)*pow(e1,2)*pow(i1,2) - 
	   2*b1*d1*d2*pow(f1,2)*pow(i1,2) + 
	   4*a2*pow(b1,2)*pow(f1,2)*pow(i1,2) + 
	   b2*pow(d1,2)*pow(f1,2)*pow(i1,2);
	rtc[1]=
		-2*c2*d1*e1*g1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2)) + 
   2*c1*d1*e2*g1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2)) + 
   4*b1*c2*f1*g1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2)) - 
   e1*e2*f1*g1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2)) - 
   4*b1*c1*f2*g1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2)) + 
   4*a1*c2*e1*h1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2)) - 
   4*a1*c1*e2*h1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2)) - 
   2*c2*d1*f1*h1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2)) + 
   2*c1*d1*f2*h1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2)) - 
   e1*f1*f2*h1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2)) - 
   8*a1*b1*c2*i1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2)) + 
   2*a1*e1*e2*i1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2)) - 
   d1*e2*f1*i1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2)) - 
   d1*e1*f2*i1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2)) + 
   2*b1*f1*f2*i1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2)) + 
   8*a1*b1*c1*i2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2)) + 
   2*d1*e1*f1*i2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2)) + 
   2*c2*i1*pow(d1,2)*(-(d2*e1*g1) - d1*e2*g1 + 
      2*b2*f1*g1 + 2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 
      2*a2*e1*h1 + 2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 
      2*a1*e1*h2 - d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 
      2*d1*d2*i1 - 4*a1*b1*i2 + i2*pow(d1,2)) - 
   2*c1*i2*pow(d1,2)*(-(d2*e1*g1) - d1*e2*g1 + 
      2*b2*f1*g1 + 2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 
      2*a2*e1*h1 + 2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 
      2*a1*e1*h2 - d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 
      2*d1*d2*i1 - 4*a1*b1*i2 + i2*pow(d1,2)) + 
   f2*g1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2))*pow(e1,2) - 
   2*a1*i2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2))*pow(e1,2) - 
   8*a2*b1*c1*g1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2)) + 
   2*c1*d1*d2*g1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2)) - 
   d2*e1*f1*g1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2)) - 
   d1*e1*f2*g1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2)) + 
   2*b1*f1*f2*g1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2)) + 
   8*a1*b1*c1*g2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2)) + 
   2*d1*e1*f1*g2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2)) + 
   4*a2*c1*d1*h1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2)) - 
   4*a1*c1*d2*h1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2)) - 
   2*a2*e1*f1*h1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2)) + 
   2*a1*e1*f2*h1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2)) - 
   d1*f1*f2*h1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2)) - 
   2*a2*d1*e1*i1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2)) + 
   2*a1*d2*e1*i1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2)) + 
   4*a2*b1*f1*i1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2)) - 
   d1*d2*f1*i1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2)) - 
   4*a1*b1*f2*i1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2)) - 
   2*c1*g2*pow(d1,2)*(-4*b2*c1*g1 - 4*b1*c2*g1 + 
      2*e1*e2*g1 - 4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - 
      e2*f1*h1 - e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - 
      d2*e1*i1 - d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - 
      d1*e1*i2 + 2*b1*f1*i2 + g2*pow(e1,2)) + 
   f2*i1*pow(d1,2)*(-4*b2*c1*g1 - 4*b1*c2*g1 + 
      2*e1*e2*g1 - 4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - 
      e2*f1*h1 - e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - 
      d2*e1*i1 - d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - 
      d1*e1*i2 + 2*b1*f1*i2 + g2*pow(e1,2)) + 
   2*a2*g1*pow(e1,2)*(-4*b2*c1*g1 - 4*b1*c2*g1 + 
      2*e1*e2*g1 - 4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - 
      e2*f1*h1 - e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - 
      d2*e1*i1 - d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - 
      d1*e1*i2 + 2*b1*f1*i2 + g2*pow(e1,2)) - 
   2*a1*g2*pow(e1,2)*(-4*b2*c1*g1 - 4*b1*c2*g1 + 
      2*e1*e2*g1 - 4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - 
      e2*f1*h1 - e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - 
      d2*e1*i1 - d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - 
      d1*e1*i2 + 2*b1*f1*i2 + g2*pow(e1,2)) + 
   e2*h1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2))*pow(f1,2) - 
   2*b1*i2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2))*pow(f1,2) - 
   2*b1*g2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2))*pow(f1,2) + 
   d2*h1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2))*pow(f1,2) - 
   4*b1*c1*g1*g2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2)) + 
   2*c1*d1*g2*h1*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2)) - 
   e1*f1*g2*h1*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2)) + 
   2*c1*d1*g1*h2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2)) - 
   e1*f1*g1*h2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2)) - 
   4*a1*c1*h1*h2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2)) - 
   d1*e1*g2*i1*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2)) + 
   2*b1*f1*g2*i1*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2)) + 
   2*a1*e1*h2*i1*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2)) - 
   d1*f1*h2*i1*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2)) - 
   d1*e1*g1*i2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2)) + 
   2*b1*f1*g1*i2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2)) + 
   2*a1*e1*h1*i2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2)) - 
   d1*f1*h1*i2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2)) - 
   4*a1*b1*i1*i2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2)) + 
   16*a1*b1*c1*j2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2)) + 
   4*d1*e1*f1*j2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2)) + 
   i1*i2*pow(d1,2)*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2)) - 
   4*c1*j2*pow(d1,2)*(8*a2*b1*c1 + 8*a1*b2*c1 + 
      8*a1*b1*c2 - 4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 
      2*d1*e2*f1 + 2*d1*e1*f2 - 4*b1*f1*f2 - 
      2*c2*pow(d1,2) - 2*a2*pow(e1,2) - 2*b2*pow(f1,2)) + 
   g1*g2*pow(e1,2)*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2)) - 
   4*a1*j2*pow(e1,2)*(8*a2*b1*c1 + 8*a1*b2*c1 + 
      8*a1*b1*c2 - 4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 
      2*d1*e2*f1 + 2*d1*e1*f2 - 4*b1*f1*f2 - 
      2*c2*pow(d1,2) - 2*a2*pow(e1,2) - 2*b2*pow(f1,2)) + 
   h1*h2*pow(f1,2)*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2)) - 
   4*b1*j2*pow(f1,2)*(8*a2*b1*c1 + 8*a1*b2*c1 + 
      8*a1*b1*c2 - 4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 
      2*d1*e2*f1 + 2*d1*e1*f2 - 4*b1*f1*f2 - 
      2*c2*pow(d1,2) - 2*a2*pow(e1,2) - 2*b2*pow(f1,2)) + 
   4*b2*c1*d1*g1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) - 
   4*b1*c1*d2*g1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) - 
   d1*e1*e2*g1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) - 
   2*b2*e1*f1*g1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) + 
   2*b1*e2*f1*g1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) - 
   8*a1*b2*c1*h1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) + 
   2*c1*d1*d2*h1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) + 
   2*a1*e1*e2*h1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) - 
   d2*e1*f1*h1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) - 
   d1*e2*f1*h1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) + 
   8*a1*b1*c1*h2*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) + 
   2*d1*e1*f1*h2*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) + 
   4*a1*b2*e1*i1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) - 
   d1*d2*e1*i1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) - 
   4*a1*b1*e2*i1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) - 
   2*b2*d1*f1*i1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) + 
   2*b1*d2*f1*i1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) - 
   2*c1*h2*pow(d1,2)*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) + 
   e2*i1*pow(d1,2)*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) + 
   d2*g1*pow(e1,2)*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) - 
   2*a1*h2*pow(e1,2)*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) + 
   2*b2*h1*pow(f1,2)*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) - 
   2*b1*h2*pow(f1,2)*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2));
	rtc[2]=
		-2*c2*d1*e1*g1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*pow(d2,2)) + 
   2*c1*d1*e2*g1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*pow(d2,2)) + 
   4*b1*c2*f1*g1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*pow(d2,2)) - 
   e1*e2*f1*g1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*pow(d2,2)) - 
   4*b1*c1*f2*g1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*pow(d2,2)) + 
   4*a1*c2*e1*h1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*pow(d2,2)) - 
   4*a1*c1*e2*h1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*pow(d2,2)) - 
   2*c2*d1*f1*h1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*pow(d2,2)) + 
   2*c1*d1*f2*h1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*pow(d2,2)) - 
   e1*f1*f2*h1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*pow(d2,2)) - 
   8*a1*b1*c2*i1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*pow(d2,2)) + 
   2*a1*e1*e2*i1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*pow(d2,2)) - 
   d1*e2*f1*i1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*pow(d2,2)) - 
   d1*e1*f2*i1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*pow(d2,2)) + 
   2*b1*f1*f2*i1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*pow(d2,2)) + 
   8*a1*b1*c1*i2*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*pow(d2,2)) + 
   2*d1*e1*f1*i2*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*pow(d2,2)) + 
   2*c2*i1*pow(d1,2)*(-(d2*e2*g1) + 2*b2*f2*g1 - 
      d2*e1*g2 - d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 
      2*a2*e2*h1 - d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - 
      d2*f1*h2 - d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 
      4*a1*b2*i2 + 2*d1*d2*i2 + i1*pow(d2,2)) - 
   2*c1*i2*pow(d1,2)*(-(d2*e2*g1) + 2*b2*f2*g1 - 
      d2*e1*g2 - d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 
      2*a2*e2*h1 - d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - 
      d2*f1*h2 - d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 
      4*a1*b2*i2 + 2*d1*d2*i2 + i1*pow(d2,2)) + 
   f2*g1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*pow(d2,2))*pow(e1,2) - 
   2*a1*i2*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*pow(d2,2))*pow(e1,2) + 
   f2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*pow(d1,2))*(-4*b2*c1*g1 - 4*b1*c2*g1 + 
      2*e1*e2*g1 - 4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - 
      e2*f1*h1 - e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - 
      d2*e1*i1 - d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - 
      d1*e1*i2 + 2*b1*f1*i2 + g2*pow(e1,2)) - 
   8*a2*b1*c1*g1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2)) + 
   2*c1*d1*d2*g1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2)) - 
   d2*e1*f1*g1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2)) - 
   d1*e1*f2*g1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2)) + 
   2*b1*f1*f2*g1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2)) + 
   8*a1*b1*c1*g2*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2)) + 
   2*d1*e1*f1*g2*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2)) + 
   4*a2*c1*d1*h1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2)) - 
   4*a1*c1*d2*h1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2)) - 
   2*a2*e1*f1*h1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2)) + 
   2*a1*e1*f2*h1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2)) - 
   d1*f1*f2*h1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2)) - 
   2*a2*d1*e1*i1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2)) + 
   2*a1*d2*e1*i1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2)) + 
   4*a2*b1*f1*i1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2)) - 
   d1*d2*f1*i1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2)) - 
   4*a1*b1*f2*i1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2)) - 
   2*c1*g2*pow(d1,2)*(-4*b2*c2*g1 - 4*b2*c1*g2 - 
      4*b1*c2*g2 + 2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 
      2*c2*d1*h2 + 2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - 
      d2*e2*i1 + 2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 
      2*b2*f1*i2 + 2*b1*f2*i2 + g1*pow(e2,2)) + 
   f2*i1*pow(d1,2)*(-4*b2*c2*g1 - 4*b2*c1*g2 - 
      4*b1*c2*g2 + 2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 
      2*c2*d1*h2 + 2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - 
      d2*e2*i1 + 2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 
      2*b2*f1*i2 + 2*b1*f2*i2 + g1*pow(e2,2)) + 
   2*a2*g1*pow(e1,2)*(-4*b2*c2*g1 - 4*b2*c1*g2 - 
      4*b1*c2*g2 + 2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 
      2*c2*d1*h2 + 2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - 
      d2*e2*i1 + 2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 
      2*b2*f1*i2 + 2*b1*f2*i2 + g1*pow(e2,2)) - 
   2*a1*g2*pow(e1,2)*(-4*b2*c2*g1 - 4*b2*c1*g2 - 
      4*b1*c2*g2 + 2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 
      2*c2*d1*h2 + 2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - 
      d2*e2*i1 + 2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 
      2*b2*f1*i2 + 2*b1*f2*i2 + g1*pow(e2,2)) + 
   e2*h1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*pow(d2,2))*pow(f1,2) - 
   2*b1*i2*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*pow(d2,2))*pow(f1,2) - 
   2*b1*g2*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2))*pow(f1,2) + 
   d2*h1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2))*pow(f1,2) + 
   i2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*pow(d1,2))*(8*a2*b1*c1 + 8*a1*b2*c1 + 
      8*a1*b1*c2 - 4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 
      2*d1*e2*f1 + 2*d1*e1*f2 - 4*b1*f1*f2 - 
      2*c2*pow(d1,2) - 2*a2*pow(e1,2) - 2*b2*pow(f1,2)) + 
   g2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2))*
    (8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 4*c1*d1*d2 - 
      4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 2*d1*e1*f2 - 
      4*b1*f1*f2 - 2*c2*pow(d1,2) - 2*a2*pow(e1,2) - 
      2*b2*pow(f1,2)) + 
   e2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*pow(d1,2))*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) + 
   d2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2))*
    (2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - e1*f2*g1 + 
      2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 4*a1*c2*h1 + 
      2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 2*a1*e2*i1 - 
      d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - d1*f1*i2 + 
      h2*pow(f1,2)) + h2*
    (8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 4*c1*d1*d2 - 
      4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 2*d1*e1*f2 - 
      4*b1*f1*f2 - 2*c2*pow(d1,2) - 2*a2*pow(e1,2) - 
      2*b2*pow(f1,2))*(2*c2*d1*g1 + 2*c1*d2*g1 - 
      e2*f1*g1 - e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 
      4*a2*c1*h1 - 4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 
      2*a2*e1*i1 + 2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 
      2*a1*e1*i2 - d1*f1*i2 + h2*pow(f1,2)) - 
   4*b1*c1*g1*g2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2)) + 
   2*c1*d1*g2*h1*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2)) - 
   e1*f1*g2*h1*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2)) + 
   2*c1*d1*g1*h2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2)) - 
   e1*f1*g1*h2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2)) - 
   4*a1*c1*h1*h2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2)) - 
   d1*e1*g2*i1*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2)) + 
   2*b1*f1*g2*i1*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2)) + 
   2*a1*e1*h2*i1*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2)) - 
   d1*f1*h2*i1*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2)) - 
   d1*e1*g1*i2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2)) + 
   2*b1*f1*g1*i2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2)) + 
   2*a1*e1*h1*i2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2)) - 
   d1*f1*h1*i2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2)) - 
   4*a1*b1*i1*i2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2)) + 
   16*a1*b1*c1*j2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2)) + 
   4*d1*e1*f1*j2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2)) + 
   i1*i2*pow(d1,2)*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2)) - 
   4*c1*j2*pow(d1,2)*(8*a2*b2*c1 + 8*a2*b1*c2 + 
      8*a1*b2*c2 - 4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 
      2*d2*e1*f2 + 2*d1*e2*f2 - 4*b2*f1*f2 - 
      2*c1*pow(d2,2) - 2*a1*pow(e2,2) - 2*b1*pow(f2,2)) + 
   g1*g2*pow(e1,2)*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2)) - 
   4*a1*j2*pow(e1,2)*(8*a2*b2*c1 + 8*a2*b1*c2 + 
      8*a1*b2*c2 - 4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 
      2*d2*e1*f2 + 2*d1*e2*f2 - 4*b2*f1*f2 - 
      2*c1*pow(d2,2) - 2*a1*pow(e2,2) - 2*b1*pow(f2,2)) + 
   h1*h2*pow(f1,2)*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2)) - 
   4*b1*j2*pow(f1,2)*(8*a2*b2*c1 + 8*a2*b1*c2 + 
      8*a1*b2*c2 - 4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 
      2*d2*e1*f2 + 2*d1*e2*f2 - 4*b2*f1*f2 - 
      2*c1*pow(d2,2) - 2*a1*pow(e2,2) - 2*b1*pow(f2,2)) + 
   4*b2*c1*d1*g1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) - 
   4*b1*c1*d2*g1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) - 
   d1*e1*e2*g1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) - 
   2*b2*e1*f1*g1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) + 
   2*b1*e2*f1*g1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) - 
   8*a1*b2*c1*h1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) + 
   2*c1*d1*d2*h1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) + 
   2*a1*e1*e2*h1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) - 
   d2*e1*f1*h1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) - 
   d1*e2*f1*h1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) + 
   8*a1*b1*c1*h2*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) + 
   2*d1*e1*f1*h2*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) + 
   4*a1*b2*e1*i1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) - 
   d1*d2*e1*i1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) - 
   4*a1*b1*e2*i1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) - 
   2*b2*d1*f1*i1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) + 
   2*b1*d2*f1*i1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) - 
   2*c1*h2*pow(d1,2)*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) + 
   e2*i1*pow(d1,2)*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) + 
   d2*g1*pow(e1,2)*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) - 
   2*a1*h2*pow(e1,2)*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) + 
   2*b2*h1*pow(f1,2)*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) - 
   2*b1*h2*pow(f1,2)*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) + 
   c2*pow(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*pow(d1,2),2) + 
   a2*pow(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2),2) + 
   j2*pow(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2),2) + 
   b2*pow(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - e1*f2*g1 + 
      2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 4*a1*c2*h1 + 
      2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 2*a1*e2*i1 - 
      d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - d1*f1*i2 + 
      h2*pow(f1,2),2);
	rtc[3]=
		2*c2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*pow(d1,2))*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*pow(d2,2)) - 
   2*c2*d1*e1*g1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) + 
   2*c1*d1*e2*g1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) + 
   4*b1*c2*f1*g1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) - 
   e1*e2*f1*g1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) - 
   4*b1*c1*f2*g1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) + 
   4*a1*c2*e1*h1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) - 
   4*a1*c1*e2*h1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) - 
   2*c2*d1*f1*h1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) + 
   2*c1*d1*f2*h1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) - 
   e1*f1*f2*h1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) - 
   8*a1*b1*c2*i1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) + 
   2*a1*e1*e2*i1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) - 
   d1*e2*f1*i1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) - 
   d1*e1*f2*i1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) + 
   2*b1*f1*f2*i1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) + 
   8*a1*b1*c1*i2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) + 
   2*d1*e1*f1*i2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) + 
   2*c2*i1*pow(d1,2)*(-(d2*e2*g2) + 2*b2*f2*g2 + 
      2*a2*e2*h2 - d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) - 
   2*c1*i2*pow(d1,2)*(-(d2*e2*g2) + 2*b2*f2*g2 + 
      2*a2*e2*h2 - d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) + 
   f2*g1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2))*pow(e1,2) - 
   2*a1*i2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2))*pow(e1,2) + 
   f2*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*pow(d2,2))*(-4*b2*c1*g1 - 4*b1*c2*g1 + 
      2*e1*e2*g1 - 4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - 
      e2*f1*h1 - e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - 
      d2*e1*i1 - d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - 
      d1*e1*i2 + 2*b1*f1*i2 + g2*pow(e1,2)) + 
   f2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*pow(d1,2))*(-4*b2*c2*g1 - 4*b2*c1*g2 - 
      4*b1*c2*g2 + 2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 
      2*c2*d1*h2 + 2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - 
      d2*e2*i1 + 2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 
      2*b2*f1*i2 + 2*b1*f2*i2 + g1*pow(e2,2)) + 
   2*a2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2))*
    (-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 2*e1*e2*g2 + 
      2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 2*c1*d2*h2 - 
      e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 2*b2*f2*i1 - 
      d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 2*b1*f2*i2 + 
      g1*pow(e2,2)) - 8*a2*b1*c1*g1*
    (-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*pow(e2,2)) + 
   2*c1*d1*d2*g1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) - 
   d2*e1*f1*g1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) - 
   d1*e1*f2*g1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) + 
   2*b1*f1*f2*g1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) + 
   8*a1*b1*c1*g2*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) + 
   2*d1*e1*f1*g2*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) + 
   4*a2*c1*d1*h1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) - 
   4*a1*c1*d2*h1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) - 
   2*a2*e1*f1*h1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) + 
   2*a1*e1*f2*h1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) - 
   d1*f1*f2*h1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) - 
   2*a2*d1*e1*i1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) + 
   2*a1*d2*e1*i1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) + 
   4*a2*b1*f1*i1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) - 
   d1*d2*f1*i1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) - 
   4*a1*b1*f2*i1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) - 
   2*c1*g2*pow(d1,2)*(-4*b2*c2*g2 + 2*c2*d2*h2 - 
      e2*f2*h2 - d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) + 
   f2*i1*pow(d1,2)*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) + 
   2*a2*g1*pow(e1,2)*(-4*b2*c2*g2 + 2*c2*d2*h2 - 
      e2*f2*h2 - d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) - 
   2*a1*g2*pow(e1,2)*(-4*b2*c2*g2 + 2*c2*d2*h2 - 
      e2*f2*h2 - d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) + 
   e2*h1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2))*pow(f1,2) - 
   2*b1*i2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2))*pow(f1,2) - 
   2*b1*g2*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2))*pow(f1,2) + 
   d2*h1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*pow(e2,2))*pow(f1,2) + 
   i2*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*pow(d2,2))*(8*a2*b1*c1 + 8*a1*b2*c1 + 
      8*a1*b1*c2 - 4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 
      2*d1*e2*f1 + 2*d1*e1*f2 - 4*b1*f1*f2 - 
      2*c2*pow(d1,2) - 2*a2*pow(e1,2) - 2*b2*pow(f1,2)) + 
   g2*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2))*
    (8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 4*c1*d1*d2 - 
      4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 2*d1*e1*f2 - 
      4*b1*f1*f2 - 2*c2*pow(d1,2) - 2*a2*pow(e1,2) - 
      2*b2*pow(f1,2)) + 
   e2*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*pow(d2,2))*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*pow(f1,2)) + 
   d2*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2))*
    (2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - e1*f2*g1 + 
      2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 4*a1*c2*h1 + 
      2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 2*a1*e2*i1 - 
      d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - d1*f1*i2 + 
      h2*pow(f1,2)) + i2*
    (-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*pow(d1,2))*(8*a2*b2*c1 + 8*a2*b1*c2 + 
      8*a1*b2*c2 - 4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 
      2*d2*e1*f2 + 2*d1*e2*f2 - 4*b2*f1*f2 - 
      2*c1*pow(d2,2) - 2*a1*pow(e2,2) - 2*b1*pow(f2,2)) + 
   g2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2))*
    (8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 4*c2*d1*d2 - 
      4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 2*d1*e2*f2 - 
      4*b2*f1*f2 - 2*c1*pow(d2,2) - 2*a1*pow(e2,2) - 
      2*b1*pow(f2,2)) + 
   2*j2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2))*
    (8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 4*c2*d1*d2 - 
      4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 2*d1*e2*f2 - 
      4*b2*f1*f2 - 2*c1*pow(d2,2) - 2*a1*pow(e2,2) - 
      2*b1*pow(f2,2)) + 
   h2*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - e1*f2*g1 + 
      2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 4*a1*c2*h1 + 
      2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 2*a1*e2*i1 - 
      d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - d1*f1*i2 + 
      h2*pow(f1,2))*(8*a2*b2*c1 + 8*a2*b1*c2 + 
      8*a1*b2*c2 - 4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 
      2*d2*e1*f2 + 2*d1*e2*f2 - 4*b2*f1*f2 - 
      2*c1*pow(d2,2) - 2*a1*pow(e2,2) - 2*b1*pow(f2,2)) - 
   4*b1*c1*g1*g2*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   2*c1*d1*g2*h1*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2)) - 
   e1*f1*g2*h1*(8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*pow(d2,2) - 
      2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   2*c1*d1*g1*h2*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2)) - 
   e1*f1*g1*h2*(8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*pow(d2,2) - 
      2*a2*pow(e2,2) - 2*b2*pow(f2,2)) - 
   4*a1*c1*h1*h2*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2)) - 
   d1*e1*g2*i1*(8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*pow(d2,2) - 
      2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   2*b1*f1*g2*i1*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   2*a1*e1*h2*i1*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2)) - 
   d1*f1*h2*i1*(8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*pow(d2,2) - 
      2*a2*pow(e2,2) - 2*b2*pow(f2,2)) - 
   d1*e1*g1*i2*(8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*pow(d2,2) - 
      2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   2*b1*f1*g1*i2*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   2*a1*e1*h1*i2*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2)) - 
   d1*f1*h1*i2*(8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*pow(d2,2) - 
      2*a2*pow(e2,2) - 2*b2*pow(f2,2)) - 
   4*a1*b1*i1*i2*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   16*a1*b1*c1*j2*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   4*d1*e1*f1*j2*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   i1*i2*pow(d1,2)*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2)) - 
   4*c1*j2*pow(d1,2)*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   g1*g2*pow(e1,2)*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2)) - 
   4*a1*j2*pow(e1,2)*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   h1*h2*pow(f1,2)*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2)) - 
   4*b1*j2*pow(f1,2)*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   e2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*pow(d1,2))*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) + 
   d2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2))*
    (2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 2*c1*d2*g2 - 
      e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 4*a2*c1*h2 - 
      4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - d2*f2*i1 + 
      2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - d1*f2*i2 + 
      h1*pow(f2,2)) + h2*
    (8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 4*c1*d1*d2 - 
      4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 2*d1*e1*f2 - 
      4*b1*f1*f2 - 2*c2*pow(d1,2) - 2*a2*pow(e1,2) - 
      2*b2*pow(f1,2))*(2*c2*d2*g1 - e2*f2*g1 + 
      2*c2*d1*g2 + 2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 
      4*a2*c2*h1 - 4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 
      2*a2*e2*i1 - d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - 
      d2*f1*i2 - d1*f2*i2 + h1*pow(f2,2)) + 
   2*b2*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - e1*f2*g1 + 
      2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 4*a1*c2*h1 + 
      2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 2*a1*e2*i1 - 
      d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - d1*f1*i2 + 
      h2*pow(f1,2))*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) + 
   4*b2*c1*d1*g1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) - 
   4*b1*c1*d2*g1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) - 
   d1*e1*e2*g1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) - 
   2*b2*e1*f1*g1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) + 
   2*b1*e2*f1*g1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) - 
   8*a1*b2*c1*h1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) + 
   2*c1*d1*d2*h1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) + 
   2*a1*e1*e2*h1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) - 
   d2*e1*f1*h1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) - 
   d1*e2*f1*h1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) + 
   8*a1*b1*c1*h2*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) + 
   2*d1*e1*f1*h2*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) + 
   4*a1*b2*e1*i1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) - 
   d1*d2*e1*i1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) - 
   4*a1*b1*e2*i1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) - 
   2*b2*d1*f1*i1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) + 
   2*b1*d2*f1*i1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) - 
   2*c1*h2*pow(d1,2)*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) + 
   e2*i1*pow(d1,2)*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) + 
   d2*g1*pow(e1,2)*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) - 
   2*a1*h2*pow(e1,2)*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) + 
   2*b2*h1*pow(f1,2)*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) - 
   2*b1*h2*pow(f1,2)*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2));
	rtc[4]=
		2*c2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*pow(d1,2))*(-(d2*e2*g2) + 2*b2*f2*g2 + 
      2*a2*e2*h2 - d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) + 
   f2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - d2*f2*h2 - 
      4*a2*b2*i2 + i2*pow(d2,2))*
    (-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 4*b1*c1*g2 + 
      2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - e1*f2*h1 + 
      2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - d1*e2*i1 + 
      2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 2*b1*f1*i2 + 
      g2*pow(e1,2)) + f2*
    (-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*pow(d2,2))*(-4*b2*c2*g1 - 4*b2*c1*g2 - 
      4*b1*c2*g2 + 2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 
      2*c2*d1*h2 + 2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - 
      d2*e2*i1 + 2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 
      2*b2*f1*i2 + 2*b1*f2*i2 + g1*pow(e2,2)) + 
   f2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*pow(d1,2))*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) + 
   2*a2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2))*
    (-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*pow(e2,2)) + 
   i2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - d2*f2*h2 - 
      4*a2*b2*i2 + i2*pow(d2,2))*
    (8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 4*c1*d1*d2 - 
      4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 2*d1*e1*f2 - 
      4*b1*f1*f2 - 2*c2*pow(d1,2) - 2*a2*pow(e1,2) - 
      2*b2*pow(f1,2)) + 
   g2*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*pow(e2,2))*
    (8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 4*c1*d1*d2 - 
      4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 2*d1*e1*f2 - 
      4*b1*f1*f2 - 2*c2*pow(d1,2) - 2*a2*pow(e1,2) - 
      2*b2*pow(f1,2)) + 
   e2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - d2*f2*h2 - 
      4*a2*b2*i2 + i2*pow(d2,2))*
    (2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - e1*f2*g1 + 
      2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 4*a1*c2*h1 + 
      2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 2*a1*e2*i1 - 
      d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - d1*f1*i2 + 
      h2*pow(f1,2)) + d2*
    (-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*pow(e2,2))*
    (2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - e1*f2*g1 + 
      2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 4*a1*c2*h1 + 
      2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 2*a1*e2*i1 - 
      d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - d1*f1*i2 + 
      h2*pow(f1,2)) + i2*
    (-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*pow(d2,2))*(8*a2*b2*c1 + 8*a2*b1*c2 + 
      8*a1*b2*c2 - 4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 
      2*d2*e1*f2 + 2*d1*e2*f2 - 4*b2*f1*f2 - 
      2*c1*pow(d2,2) - 2*a1*pow(e2,2) - 2*b1*pow(f2,2)) + 
   g2*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2))*
    (8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 4*c2*d1*d2 - 
      4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 2*d1*e2*f2 - 
      4*b2*f1*f2 - 2*c1*pow(d2,2) - 2*a1*pow(e2,2) - 
      2*b1*pow(f2,2)) + 
   i2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*pow(d1,2))*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   g2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2))*
    (8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*pow(d2,2) - 
      2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   2*j2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*pow(d1,2) - 
      2*a2*pow(e1,2) - 2*b2*pow(f1,2))*
    (8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*pow(d2,2) - 
      2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   h2*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - e1*f2*g1 + 
      2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 4*a1*c2*h1 + 
      2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 2*a1*e2*i1 - 
      d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - d1*f1*i2 + 
      h2*pow(f1,2))*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   e2*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*pow(d2,2))*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2)) + 
   d2*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2))*
    (2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 2*c1*d2*g2 - 
      e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 4*a2*c1*h2 - 
      4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - d2*f2*i1 + 
      2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - d1*f2*i2 + 
      h1*pow(f2,2)) + h2*
    (8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 4*c2*d1*d2 - 
      4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 2*d1*e2*f2 - 
      4*b2*f1*f2 - 2*c1*pow(d2,2) - 2*a1*pow(e2,2) - 
      2*b1*pow(f2,2))*(2*c2*d2*g1 - e2*f2*g1 + 
      2*c2*d1*g2 + 2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 
      4*a2*c2*h1 - 4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 
      2*a2*e2*i1 - d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - 
      d2*f1*i2 - d1*f2*i2 + h1*pow(f2,2)) + 
   e2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*pow(d1,2))*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) + 
   d2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*pow(e1,2))*
    (2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 2*a2*e2*i2 - 
      d2*f2*i2 + h2*pow(f2,2)) + 
   h2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 4*c1*d1*d2 - 
      4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 2*d1*e1*f2 - 
      4*b1*f1*f2 - 2*c2*pow(d1,2) - 2*a2*pow(e1,2) - 
      2*b2*pow(f1,2))*(2*c2*d2*g2 - e2*f2*g2 - 
      4*a2*c2*h2 + 2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) + 
   2*b2*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - e1*f2*g1 + 
      2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 4*a1*c2*h1 + 
      2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 2*a1*e2*i1 - 
      d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - d1*f1*i2 + 
      h2*pow(f1,2))*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) + 
   c2*pow(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*pow(d2,2),2) + 
   a2*pow(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2),2) + 
   j2*pow(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2),2) + 
   b2*pow(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*pow(f2,2),2);
	rtc[5]=
		2*c2*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*pow(d2,2))*(-(d2*e2*g2) + 2*b2*f2*g2 + 
      2*a2*e2*h2 - d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2)) + 
   f2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - d2*f2*h2 - 
      4*a2*b2*i2 + i2*pow(d2,2))*
    (-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 2*e1*e2*g2 + 
      2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 2*c1*d2*h2 - 
      e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 2*b2*f2*i1 - 
      d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 2*b1*f2*i2 + 
      g1*pow(e2,2)) + f2*
    (-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*pow(d2,2))*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*pow(e2,2)) + 
   2*a2*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2))*
    (-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*pow(e2,2)) + 
   i2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - d2*f2*h2 - 
      4*a2*b2*i2 + i2*pow(d2,2))*
    (8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 4*c2*d1*d2 - 
      4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 2*d1*e2*f2 - 
      4*b2*f1*f2 - 2*c1*pow(d2,2) - 2*a1*pow(e2,2) - 
      2*b1*pow(f2,2)) + 
   g2*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*pow(e2,2))*
    (8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 4*c2*d1*d2 - 
      4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 2*d1*e2*f2 - 
      4*b2*f1*f2 - 2*c1*pow(d2,2) - 2*a1*pow(e2,2) - 
      2*b1*pow(f2,2)) + 
   i2*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*pow(d2,2))*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   g2*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2))*
    (8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*pow(d2,2) - 
      2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   2*j2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*pow(d2,2) - 
      2*a1*pow(e2,2) - 2*b1*pow(f2,2))*
    (8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*pow(d2,2) - 
      2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   e2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - d2*f2*h2 - 
      4*a2*b2*i2 + i2*pow(d2,2))*
    (2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 2*c1*d2*g2 - 
      e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 4*a2*c1*h2 - 
      4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - d2*f2*i1 + 
      2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - d1*f2*i2 + 
      h1*pow(f2,2)) + d2*
    (-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*pow(e2,2))*
    (2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 2*c1*d2*g2 - 
      e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 4*a2*c1*h2 - 
      4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - d2*f2*i1 + 
      2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - d1*f2*i2 + 
      h1*pow(f2,2)) + h2*
    (8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*pow(d2,2) - 
      2*a2*pow(e2,2) - 2*b2*pow(f2,2))*
    (2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 2*c1*d2*g2 - 
      e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 4*a2*c1*h2 - 
      4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - d2*f2*i1 + 
      2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - d1*f2*i2 + 
      h1*pow(f2,2)) + e2*
    (-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*pow(d2,2))*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) + 
   d2*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*pow(e2,2))*
    (2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 2*a2*e2*i2 - 
      d2*f2*i2 + h2*pow(f2,2)) + 
   h2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 4*c2*d1*d2 - 
      4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 2*d1*e2*f2 - 
      4*b2*f1*f2 - 2*c1*pow(d2,2) - 2*a1*pow(e2,2) - 
      2*b1*pow(f2,2))*(2*c2*d2*g2 - e2*f2*g2 - 
      4*a2*c2*h2 + 2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2)) + 
   2*b2*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 2*c1*d2*g2 - 
      e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 4*a2*c1*h2 - 
      4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - d2*f2*i1 + 
      2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - d1*f2*i2 + 
      h1*pow(f2,2))*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2));
	rtc[6]=
		f2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - d2*f2*h2 - 
      4*a2*b2*i2 + i2*pow(d2,2))*
    (-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*pow(e2,2)) + 
   i2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - d2*f2*h2 - 
      4*a2*b2*i2 + i2*pow(d2,2))*
    (8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*pow(d2,2) - 
      2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   g2*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*pow(e2,2))*
    (8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*pow(d2,2) - 
      2*a2*pow(e2,2) - 2*b2*pow(f2,2)) + 
   e2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - d2*f2*h2 - 
      4*a2*b2*i2 + i2*pow(d2,2))*
    (2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 2*a2*e2*i2 - 
      d2*f2*i2 + h2*pow(f2,2)) + 
   d2*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*pow(e2,2))*
    (2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 2*a2*e2*i2 - 
      d2*f2*i2 + h2*pow(f2,2)) + 
   h2*(8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*pow(d2,2) - 
      2*a2*pow(e2,2) - 2*b2*pow(f2,2))*
    (2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 2*a2*e2*i2 - 
      d2*f2*i2 + h2*pow(f2,2)) + 
   c2*pow(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*pow(d2,2),2) + 
   a2*pow(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*pow(e2,2),2) + 
   j2*pow(8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*pow(d2,2) - 
      2*a2*pow(e2,2) - 2*b2*pow(f2,2),2) + 
   b2*pow(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*pow(f2,2),2);

	int order=0;
	for (int k=0;k<=6;k++){
	  if (fabs(rtc[k]) > EPS){
		order=k;
	    }
	}
	if (order==0)
	    return false;

	for (int k=0;k<=order;k++)
	    rtc[k]/=rtc[order];

#ifndef NDEBUG  // tested: order==6 whether or not in contact
	g_debuginf<<"root6.cpp: g_iteration="<<g_iteration
		  <<" order="<<order<<std::endl;
#endif

	if (!zrhqr(rtc, order, rtr, rti)) // find roots for a polynomial using eigenvalues method
	    return false;

#ifndef NDEBUG
	for (int k=1;k<=order;k++){
	  g_debuginf<<"root6.cpp: rtr="<<rtr[k]<<" rti="<<rti[k]<<std::endl;
	}
#endif

	REAL lamda[6];
	int jj=0;
	for (int k=1;k<=order;k++){
	    if(fabs(rti[k]) < EPS)
		lamda[jj++]=rtr[k];
	}

#ifndef NDEBUG
	g_debuginf<<"root6.cpp: jj="<<jj<<std::endl;
#endif

	REAL x, y, z, det, within, itself;
	bool found=false;
	REAL deepest = 0;
	for (int k=0;k<jj;k++){
	    x=0;
	    y=0;
	    z=0;
	    within=0;
            det=
	                8*a1*b1*c1 + 2*d1*e1*f1 - 2*c1*pow(d1,2) - 
			2*a1*pow(e1,2) - 2*b1*pow(f1,2) + 
			lamda[k]*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
			4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 
			2*d1*e2*f1 + 2*d1*e1*f2 - 4*b1*f1*f2 - 
			2*c2*pow(d1,2) - 2*a2*pow(e1,2) - 2*b2*pow(f1,2))\
			+ (8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
			4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 
			2*d2*e1*f2 + 2*d1*e2*f2 - 4*b2*f1*f2 - 
			2*c1*pow(d2,2) - 2*a1*pow(e2,2) - 2*b1*pow(f2,2))*
			pow(lamda[k],2) + (8*a2*b2*c2 + 2*d2*e2*f2 - 
			2*c2*pow(d2,2) - 2*a2*pow(e2,2) - 2*b2*pow(f2,2))*
	                pow(lamda[k],3);
	    if(fabs(det) > EPS){ // if determinant is zero, there are infinite solutions
		                    // it is necessary to skip this case, otherwise computatation will fail.
		x=
			((-4*b1*c1*g1 + 2*c1*d1*h1 - e1*f1*h1 - d1*e1*i1 + 
			2*b1*f1*i1 + g1*pow(e1,2) + 
			lamda[k]*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
			4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
			e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
			d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
			2*b1*f1*i2 + g2*pow(e1,2)) + 
			(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 2*e1*e2*g2 + 
			2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 2*c1*d2*h2 - 
			e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 2*b2*f2*i1 - 
			d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 2*b1*f2*i2 + 
			g1*pow(e2,2))*pow(lamda[k],2) + 
			(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
			 2*b2*f2*i2 + g2*pow(e2,2))*pow(lamda[k],3)) )/det;
		y=
			((2*c1*d1*g1 - e1*f1*g1 - 4*a1*c1*h1 + 2*a1*e1*i1 - 
			d1*f1*i1 + h1*pow(f1,2) + 
			lamda[k]*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
			e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
			4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 
			2*a2*e1*i1 + 2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 
			2*a1*e1*i2 - d1*f1*i2 + h2*pow(f1,2)) + 
			(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 2*c1*d2*g2 - 
			e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 4*a2*c1*h2 - 
			4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - d2*f2*i1 + 
			2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - d1*f2*i2 + 
			h1*pow(f2,2))*pow(lamda[k],2) + 
			(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 2*a2*e2*i2 - 
			 d2*f2*i2 + h2*pow(f2,2))*pow(lamda[k],3)) )/det;
		z=
			((-(d1*e1*g1) + 2*b1*f1*g1 + 2*a1*e1*h1 - d1*f1*h1 - 
			4*a1*b1*i1 + i1*pow(d1,2) + 
			lamda[k]*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
			2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
			2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
			d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
			4*a1*b1*i2 + i2*pow(d1,2)) + 
			(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
			2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
			2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
			4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
			2*d1*d2*i2 + i1*pow(d2,2))*pow(lamda[k],2) + 
			(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - d2*f2*h2 - 
			 4*a2*b2*i2 + i2*pow(d2,2))*pow(lamda[k],3)) )/det;
	    
		within = a1*x*x+b1*y*y+c1*z*z+d1*x*y+e1*y*z+f1*z*x+g1*x+h1*y+i1*z+j1;
		itself = a2*x*x+b2*y*y+c2*z*z+d2*x*y+e2*y*z+f2*z*x+g2*x+h2*y+i2*z+j2;
		
#ifndef NDEBUG
		g_debuginf<<"root6.cpp: k="<<k
			  <<" det="<<det
			  <<" x="<< x <<" y="<< y <<" z="<< z
			  <<" within="<<within
			  <<" itself="<<itself<<std::endl;
#endif

		// first, the found point must be on the surface of particle 2, i.e., itself < EPS
		// second, the point must penetrate deepest into particle 1, i.e., smallest negative within 
		if(fabs(itself) < EPS && within < deepest){
		  deepest = within;
		  point=vec(x,y,z);
		  found=true;
#ifndef NDEBUG
		  g_debuginf<<"root6.cpp: selected within="<<deepest<<std::endl;
#endif
		}
	    }
	    else {
#ifndef NDEBUG
	      g_debuginf<<"root6: k="<<k
			<<" det="<<det
			<<" determinant is 0!" <<std::endl;
#endif
	    }
	}
	return found;
}

} // namespace dem ends
