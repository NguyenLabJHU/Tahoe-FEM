//    <<Numerical Recipes in C>>
//    void zrhqr(float a[], int m, float rtr[], float rti[])
//    Find all the roots of a polynomial with real coefficients, i=0~m, a(i)xi, given the degree m
//    and the coefficients a[0..m]. The method is to construct an upper Hessenberg matrix whose
//    eigenvalues are the desired roots, and then use the routines balanc and hqr. The real and
//    imaginary parts of the roots are returned in rtr[1..m] and rti[1..m], respectively

//    The purpose of this function is to find the point on the surface of particle 2, which is innermost to
//    the surface of particle 1. However, the algorithm in the thesis can't guarantee that the found point
//    is inside particle 1, which means, even if the two particles are completely separated (not-in-touch),
//    the algorithm could still find a point!
//
//    return values:
//      true  - non-overlapped
//      false - overlapped and vector point returns the innermost point.

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "const.h"
#include "root6.h"
#include "nr.h"
//#define DEBUG
using namespace std;

namespace dem {

extern ofstream g_exceptioninf;
extern int g_iteration;

bool root6(long double coef1[],long double coef2[],vec& point){
	if((coef1[0]==1.0&&coef1[1]==1.0&&coef1[2]==1.0&&
		coef1[3]==0.0&&coef1[4]==0.0&&coef1[5]==0.0)&&
		(coef2[0]==1.0&&coef2[1]==1.0&&coef2[2]==1.0&&
		coef2[3]==0.0&&coef2[4]==0.0&&coef2[5]==0.0)){
		long double x1=-coef1[6]/2;
		long double y1=-coef1[7]/2;
		long double z1=-coef1[8]/2;
		long double R1=sqrtl(x1*x1+y1*y1+z1*z1-coef1[9]);
		long double x2=-coef2[6]/2;
		long double y2=-coef2[7]/2;
		long double z2=-coef2[8]/2;
		long double R2=sqrtl(x2*x2+y2*y2+z2*z2-coef2[9]);
		vec dirc=vec(x1-x2,y1-y2,z1-z2);
		dirc=normalize(dirc);
		point=vec(x2,y2,z2)+R2*dirc;
		vec judge=point-vec(x1,y1,z1);
		if(vfabsl(judge)<R1 && judge%dirc<0) // ptcl_1 is outside ptcl_2 and is in touch with ptcl_2.
		    return true;  // overlapped!    // if fabsl(judge)>R1 and judge%dirc>0, then ptcl_1 is inside ptcl_2.
		else
		    return false; // non-overlapped!
	}

	long double a1=coef1[0];
	long double b1=coef1[1];
	long double c1=coef1[2];
	long double d1=coef1[3];
	long double e1=coef1[4];
	long double f1=coef1[5];
	long double g1=coef1[6];
	long double h1=coef1[7];
	long double i1=coef1[8];
	long double j1=coef1[9];

	long double a2=coef2[0];
	long double b2=coef2[1];
	long double c2=coef2[2];
	long double d2=coef2[3];
	long double e2=coef2[4];
	long double f2=coef2[5];
	long double g2=coef2[6];
	long double h2=coef2[7];
	long double i2=coef2[8];
	long double j2=coef2[9];
	long double rtc[7], rti[7], rtr[7];
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
	   16*b2*c1*e1*h1*i1*powl(a1,2) - 
	   16*b1*c2*e1*h1*i1*powl(a1,2) + 
	   16*b1*c1*e2*h1*i1*powl(a1,2) + 
	   16*b1*c1*e1*h2*i1*powl(a1,2) + 
	   16*b1*c1*e1*h1*i2*powl(a1,2) - 
	   16*a2*c1*f1*g1*i1*powl(b1,2) - 
	   16*a1*c2*f1*g1*i1*powl(b1,2) + 
	   16*a1*c1*f2*g1*i1*powl(b1,2) + 
	   16*a1*c1*f1*g2*i1*powl(b1,2) + 
	   16*a1*c1*f1*g1*i2*powl(b1,2) - 
	   32*c1*i1*i2*powl(a1,2)*powl(b1,2) - 
	   16*a2*b1*d1*g1*h1*powl(c1,2) - 
	   16*a1*b2*d1*g1*h1*powl(c1,2) + 
	   16*a1*b1*d2*g1*h1*powl(c1,2) + 
	   16*a1*b1*d1*g2*h1*powl(c1,2) + 
	   16*a1*b1*d1*g1*h2*powl(c1,2) - 
	   32*b1*h1*h2*powl(a1,2)*powl(c1,2) - 
	   32*a1*g1*g2*powl(b1,2)*powl(c1,2) + 
	   64*j2*powl(a1,2)*powl(b1,2)*powl(c1,2) + 
	   2*c2*e1*f1*g1*h1*powl(d1,2) - 
	   2*c1*e2*f1*g1*h1*powl(d1,2) - 
	   2*c1*e1*f2*g1*h1*powl(d1,2) + 
	   6*c1*e1*f1*g2*h1*powl(d1,2) + 
	   6*c1*e1*f1*g1*h2*powl(d1,2) - 
	   2*c1*d2*e1*g1*i1*powl(d1,2) - 
	   4*b2*c1*f1*g1*i1*powl(d1,2) + 
	   4*b1*c2*f1*g1*i1*powl(d1,2) - 
	   4*b1*c1*f2*g1*i1*powl(d1,2) - 
	   4*b1*c1*f1*g2*i1*powl(d1,2) - 
	   4*a2*c1*e1*h1*i1*powl(d1,2) + 
	   4*a1*c2*e1*h1*i1*powl(d1,2) - 
	   4*a1*c1*e2*h1*i1*powl(d1,2) - 
	   2*c1*d2*f1*h1*i1*powl(d1,2) - 
	   4*a1*c1*e1*h2*i1*powl(d1,2) - 
	   4*b1*c1*f1*g1*i2*powl(d1,2) - 
	   4*a1*c1*e1*h1*i2*powl(d1,2) + 
	   16*a1*b1*c1*i1*i2*powl(d1,2) + 
	   8*b1*g1*g2*powl(c1,2)*powl(d1,2) + 
	   4*d2*g1*h1*powl(c1,2)*powl(d1,2) + 
	   8*a1*h1*h2*powl(c1,2)*powl(d1,2) - 
	   32*a1*b1*j2*powl(c1,2)*powl(d1,2) - 
	   2*c2*e1*g1*i1*powl(d1,3) + 2*c1*e2*g1*i1*powl(d1,3) + 
	   2*c1*e1*g2*i1*powl(d1,3) - 2*c2*f1*h1*i1*powl(d1,3) + 
	   2*c1*f2*h1*i1*powl(d1,3) + 2*c1*f1*h2*i1*powl(d1,3) + 
	   2*c1*e1*g1*i2*powl(d1,3) + 2*c1*f1*h1*i2*powl(d1,3) + 
	   2*e1*f1*i1*i2*powl(d1,3) - 8*c1*e1*f1*j2*powl(d1,3) - 
	   4*g2*h1*powl(c1,2)*powl(d1,3) - 
	   4*g1*h2*powl(c1,2)*powl(d1,3) - 2*c1*i1*i2*powl(d1,4) + 
	   4*j2*powl(c1,2)*powl(d1,4) + 
	   16*a1*b1*c1*g1*g2*powl(e1,2) + 
	   4*a2*c1*d1*g1*h1*powl(e1,2) - 
	   4*a1*c2*d1*g1*h1*powl(e1,2) - 
	   4*a1*c1*d2*g1*h1*powl(e1,2) - 
	   2*a1*e2*f1*g1*h1*powl(e1,2) - 
	   4*a1*c1*d1*g2*h1*powl(e1,2) - 
	   4*a1*c1*d1*g1*h2*powl(e1,2) - 
	   2*a1*d1*e2*g1*i1*powl(e1,2) + 
	   4*a2*b1*f1*g1*i1*powl(e1,2) - 
	   4*a1*b2*f1*g1*i1*powl(e1,2) - 
	   4*a1*b1*f2*g1*i1*powl(e1,2) - 
	   4*a1*b1*f1*g2*i1*powl(e1,2) + 
	   2*a2*d1*f1*h1*i1*powl(e1,2) - 
	   2*a1*d2*f1*h1*i1*powl(e1,2) - 
	   2*a1*d1*f2*h1*i1*powl(e1,2) + 
	   6*a1*d1*f1*h2*i1*powl(e1,2) - 
	   4*a1*b1*f1*g1*i2*powl(e1,2) + 
	   6*a1*d1*f1*h1*i2*powl(e1,2) + 
	   8*c1*h1*h2*powl(a1,2)*powl(e1,2) + 
	   4*e2*h1*i1*powl(a1,2)*powl(e1,2) + 
	   8*b1*i1*i2*powl(a1,2)*powl(e1,2) - 
	   32*b1*c1*j2*powl(a1,2)*powl(e1,2) - 
	   2*c1*g1*g2*powl(d1,2)*powl(e1,2) + 
	   2*f2*g1*i1*powl(d1,2)*powl(e1,2) - 
	   2*f1*g2*i1*powl(d1,2)*powl(e1,2) - 
	   2*f1*g1*i2*powl(d1,2)*powl(e1,2) - 
	   2*a1*i1*i2*powl(d1,2)*powl(e1,2) + 
	   8*a1*c1*j2*powl(d1,2)*powl(e1,2) + 
	   2*d1*f1*g1*g2*powl(e1,3) - 2*a2*f1*g1*h1*powl(e1,3) + 
	   2*a1*f2*g1*h1*powl(e1,3) + 2*a1*f1*g2*h1*powl(e1,3) + 
	   2*a1*f1*g1*h2*powl(e1,3) - 2*a2*d1*g1*i1*powl(e1,3) + 
	   2*a1*d2*g1*i1*powl(e1,3) + 2*a1*d1*g2*i1*powl(e1,3) + 
	   2*a1*d1*g1*i2*powl(e1,3) - 8*a1*d1*f1*j2*powl(e1,3) - 
	   4*h2*i1*powl(a1,2)*powl(e1,3) - 
	   4*h1*i2*powl(a1,2)*powl(e1,3) - 2*a1*g1*g2*powl(e1,4) + 
	   4*j2*powl(a1,2)*powl(e1,4) + 4*b2*c1*d1*g1*h1*powl(f1,2) - 
	   4*b1*c2*d1*g1*h1*powl(f1,2) - 
	   4*b1*c1*d2*g1*h1*powl(f1,2) - 
	   2*b1*e1*f2*g1*h1*powl(f1,2) - 
	   4*b1*c1*d1*g2*h1*powl(f1,2) - 
	   4*b1*c1*d1*g1*h2*powl(f1,2) + 
	   16*a1*b1*c1*h1*h2*powl(f1,2) + 
	   2*b2*d1*e1*g1*i1*powl(f1,2) - 
	   2*b1*d2*e1*g1*i1*powl(f1,2) - 
	   2*b1*d1*e2*g1*i1*powl(f1,2) + 
	   6*b1*d1*e1*g2*i1*powl(f1,2) - 
	   4*a2*b1*e1*h1*i1*powl(f1,2) + 
	   4*a1*b2*e1*h1*i1*powl(f1,2) - 
	   4*a1*b1*e2*h1*i1*powl(f1,2) - 
	   2*b1*d1*f2*h1*i1*powl(f1,2) - 
	   4*a1*b1*e1*h2*i1*powl(f1,2) + 
	   6*b1*d1*e1*g1*i2*powl(f1,2) - 
	   4*a1*b1*e1*h1*i2*powl(f1,2) + 
	   8*c1*g1*g2*powl(b1,2)*powl(f1,2) + 
	   4*f2*g1*i1*powl(b1,2)*powl(f1,2) + 
	   8*a1*i1*i2*powl(b1,2)*powl(f1,2) - 
	   32*a1*c1*j2*powl(b1,2)*powl(f1,2) - 
	   2*c1*h1*h2*powl(d1,2)*powl(f1,2) + 
	   2*e2*h1*i1*powl(d1,2)*powl(f1,2) - 
	   2*e1*h2*i1*powl(d1,2)*powl(f1,2) - 
	   2*e1*h1*i2*powl(d1,2)*powl(f1,2) - 
	   2*b1*i1*i2*powl(d1,2)*powl(f1,2) + 
	   8*b1*c1*j2*powl(d1,2)*powl(f1,2) - 
	   2*b1*g1*g2*powl(e1,2)*powl(f1,2) + 
	   2*d2*g1*h1*powl(e1,2)*powl(f1,2) - 
	   2*d1*g2*h1*powl(e1,2)*powl(f1,2) - 
	   2*d1*g1*h2*powl(e1,2)*powl(f1,2) - 
	   2*a1*h1*h2*powl(e1,2)*powl(f1,2) + 
	   8*a1*b1*j2*powl(e1,2)*powl(f1,2) + 
	   4*j2*powl(d1,2)*powl(e1,2)*powl(f1,2) - 
	   2*b2*e1*g1*h1*powl(f1,3) + 2*b1*e2*g1*h1*powl(f1,3) + 
	   2*b1*e1*g2*h1*powl(f1,3) + 2*b1*e1*g1*h2*powl(f1,3) + 
	   2*d1*e1*h1*h2*powl(f1,3) - 2*b2*d1*h1*i1*powl(f1,3) + 
	   2*b1*d2*h1*i1*powl(f1,3) + 2*b1*d1*h2*i1*powl(f1,3) + 
	   2*b1*d1*h1*i2*powl(f1,3) - 8*b1*d1*e1*j2*powl(f1,3) - 
	   4*g2*i1*powl(b1,2)*powl(f1,3) - 
	   4*g1*i2*powl(b1,2)*powl(f1,3) - 2*b1*h1*h2*powl(f1,4) + 
	   4*j2*powl(b1,2)*powl(f1,4) - 4*b2*c1*d1*e1*f1*powl(g1,2) - 
	   4*b1*c2*d1*e1*f1*powl(g1,2) + 
	   4*b1*c1*d2*e1*f1*powl(g1,2) + 
	   4*b1*c1*d1*e2*f1*powl(g1,2) + 
	   4*b1*c1*d1*e1*f2*powl(g1,2) - 
	   8*c1*f1*f2*powl(b1,2)*powl(g1,2) - 
	   8*b1*d1*d2*powl(c1,2)*powl(g1,2) + 
	   16*a2*powl(b1,2)*powl(c1,2)*powl(g1,2) - 
	   2*c1*e1*e2*powl(d1,2)*powl(g1,2) + 
	   4*b2*powl(c1,2)*powl(d1,2)*powl(g1,2) - 
	   8*a2*b1*c1*powl(e1,2)*powl(g1,2) + 
	   2*c1*d1*d2*powl(e1,2)*powl(g1,2) + 
	   d1*e2*f1*powl(e1,2)*powl(g1,2) + 
	   2*b1*f1*f2*powl(e1,2)*powl(g1,2) + 
	   c2*powl(d1,2)*powl(e1,2)*powl(g1,2) - 
	   d2*f1*powl(e1,3)*powl(g1,2) - d1*f2*powl(e1,3)*powl(g1,2) + 
	   a2*powl(e1,4)*powl(g1,2) - 
	   2*b1*e1*e2*powl(f1,2)*powl(g1,2) + 
	   4*c2*powl(b1,2)*powl(f1,2)*powl(g1,2) + 
	   b2*powl(e1,2)*powl(f1,2)*powl(g1,2) - 
	   4*a2*c1*d1*e1*f1*powl(h1,2) - 
	   4*a1*c2*d1*e1*f1*powl(h1,2) + 
	   4*a1*c1*d2*e1*f1*powl(h1,2) + 
	   4*a1*c1*d1*e2*f1*powl(h1,2) + 
	   4*a1*c1*d1*e1*f2*powl(h1,2) - 
	   8*c1*e1*e2*powl(a1,2)*powl(h1,2) - 
	   8*a1*d1*d2*powl(c1,2)*powl(h1,2) + 
	   16*b2*powl(a1,2)*powl(c1,2)*powl(h1,2) - 
	   2*c1*f1*f2*powl(d1,2)*powl(h1,2) + 
	   4*a2*powl(c1,2)*powl(d1,2)*powl(h1,2) - 
	   2*a1*f1*f2*powl(e1,2)*powl(h1,2) + 
	   4*c2*powl(a1,2)*powl(e1,2)*powl(h1,2) - 
	   8*a1*b2*c1*powl(f1,2)*powl(h1,2) + 
	   2*c1*d1*d2*powl(f1,2)*powl(h1,2) + 
	   2*a1*e1*e2*powl(f1,2)*powl(h1,2) + 
	   d1*e1*f2*powl(f1,2)*powl(h1,2) + 
	   c2*powl(d1,2)*powl(f1,2)*powl(h1,2) + 
	   a2*powl(e1,2)*powl(f1,2)*powl(h1,2) - 
	   d2*e1*powl(f1,3)*powl(h1,2) - d1*e2*powl(f1,3)*powl(h1,2) + 
	   b2*powl(f1,4)*powl(h1,2) - 4*a2*b1*d1*e1*f1*powl(i1,2) - 
	   4*a1*b2*d1*e1*f1*powl(i1,2) + 
	   4*a1*b1*d2*e1*f1*powl(i1,2) + 
	   4*a1*b1*d1*e2*f1*powl(i1,2) + 
	   4*a1*b1*d1*e1*f2*powl(i1,2) - 
	   8*b1*e1*e2*powl(a1,2)*powl(i1,2) - 
	   8*a1*f1*f2*powl(b1,2)*powl(i1,2) + 
	   16*c2*powl(a1,2)*powl(b1,2)*powl(i1,2) - 
	   8*a1*b1*c2*powl(d1,2)*powl(i1,2) + 
	   2*a1*e1*e2*powl(d1,2)*powl(i1,2) + 
	   d2*e1*f1*powl(d1,2)*powl(i1,2) + 
	   2*b1*f1*f2*powl(d1,2)*powl(i1,2) - 
	   e2*f1*powl(d1,3)*powl(i1,2) - e1*f2*powl(d1,3)*powl(i1,2) + 
	   c2*powl(d1,4)*powl(i1,2) - 
	   2*a1*d1*d2*powl(e1,2)*powl(i1,2) + 
	   4*b2*powl(a1,2)*powl(e1,2)*powl(i1,2) + 
	   a2*powl(d1,2)*powl(e1,2)*powl(i1,2) - 
	   2*b1*d1*d2*powl(f1,2)*powl(i1,2) + 
	   4*a2*powl(b1,2)*powl(f1,2)*powl(i1,2) + 
	   b2*powl(d1,2)*powl(f1,2)*powl(i1,2);
	rtc[1]=
		-2*c2*d1*e1*g1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2)) + 
   2*c1*d1*e2*g1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2)) + 
   4*b1*c2*f1*g1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2)) - 
   e1*e2*f1*g1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2)) - 
   4*b1*c1*f2*g1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2)) + 
   4*a1*c2*e1*h1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2)) - 
   4*a1*c1*e2*h1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2)) - 
   2*c2*d1*f1*h1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2)) + 
   2*c1*d1*f2*h1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2)) - 
   e1*f1*f2*h1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2)) - 
   8*a1*b1*c2*i1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2)) + 
   2*a1*e1*e2*i1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2)) - 
   d1*e2*f1*i1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2)) - 
   d1*e1*f2*i1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2)) + 
   2*b1*f1*f2*i1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2)) + 
   8*a1*b1*c1*i2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2)) + 
   2*d1*e1*f1*i2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2)) + 
   2*c2*i1*powl(d1,2)*(-(d2*e1*g1) - d1*e2*g1 + 
      2*b2*f1*g1 + 2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 
      2*a2*e1*h1 + 2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 
      2*a1*e1*h2 - d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 
      2*d1*d2*i1 - 4*a1*b1*i2 + i2*powl(d1,2)) - 
   2*c1*i2*powl(d1,2)*(-(d2*e1*g1) - d1*e2*g1 + 
      2*b2*f1*g1 + 2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 
      2*a2*e1*h1 + 2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 
      2*a1*e1*h2 - d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 
      2*d1*d2*i1 - 4*a1*b1*i2 + i2*powl(d1,2)) + 
   f2*g1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2))*powl(e1,2) - 
   2*a1*i2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2))*powl(e1,2) - 
   8*a2*b1*c1*g1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2)) + 
   2*c1*d1*d2*g1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2)) - 
   d2*e1*f1*g1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2)) - 
   d1*e1*f2*g1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2)) + 
   2*b1*f1*f2*g1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2)) + 
   8*a1*b1*c1*g2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2)) + 
   2*d1*e1*f1*g2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2)) + 
   4*a2*c1*d1*h1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2)) - 
   4*a1*c1*d2*h1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2)) - 
   2*a2*e1*f1*h1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2)) + 
   2*a1*e1*f2*h1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2)) - 
   d1*f1*f2*h1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2)) - 
   2*a2*d1*e1*i1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2)) + 
   2*a1*d2*e1*i1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2)) + 
   4*a2*b1*f1*i1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2)) - 
   d1*d2*f1*i1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2)) - 
   4*a1*b1*f2*i1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2)) - 
   2*c1*g2*powl(d1,2)*(-4*b2*c1*g1 - 4*b1*c2*g1 + 
      2*e1*e2*g1 - 4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - 
      e2*f1*h1 - e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - 
      d2*e1*i1 - d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - 
      d1*e1*i2 + 2*b1*f1*i2 + g2*powl(e1,2)) + 
   f2*i1*powl(d1,2)*(-4*b2*c1*g1 - 4*b1*c2*g1 + 
      2*e1*e2*g1 - 4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - 
      e2*f1*h1 - e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - 
      d2*e1*i1 - d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - 
      d1*e1*i2 + 2*b1*f1*i2 + g2*powl(e1,2)) + 
   2*a2*g1*powl(e1,2)*(-4*b2*c1*g1 - 4*b1*c2*g1 + 
      2*e1*e2*g1 - 4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - 
      e2*f1*h1 - e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - 
      d2*e1*i1 - d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - 
      d1*e1*i2 + 2*b1*f1*i2 + g2*powl(e1,2)) - 
   2*a1*g2*powl(e1,2)*(-4*b2*c1*g1 - 4*b1*c2*g1 + 
      2*e1*e2*g1 - 4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - 
      e2*f1*h1 - e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - 
      d2*e1*i1 - d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - 
      d1*e1*i2 + 2*b1*f1*i2 + g2*powl(e1,2)) + 
   e2*h1*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2))*powl(f1,2) - 
   2*b1*i2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2))*powl(f1,2) - 
   2*b1*g2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2))*powl(f1,2) + 
   d2*h1*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2))*powl(f1,2) - 
   4*b1*c1*g1*g2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2)) + 
   2*c1*d1*g2*h1*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2)) - 
   e1*f1*g2*h1*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2)) + 
   2*c1*d1*g1*h2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2)) - 
   e1*f1*g1*h2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2)) - 
   4*a1*c1*h1*h2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2)) - 
   d1*e1*g2*i1*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2)) + 
   2*b1*f1*g2*i1*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2)) + 
   2*a1*e1*h2*i1*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2)) - 
   d1*f1*h2*i1*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2)) - 
   d1*e1*g1*i2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2)) + 
   2*b1*f1*g1*i2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2)) + 
   2*a1*e1*h1*i2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2)) - 
   d1*f1*h1*i2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2)) - 
   4*a1*b1*i1*i2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2)) + 
   16*a1*b1*c1*j2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2)) + 
   4*d1*e1*f1*j2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2)) + 
   i1*i2*powl(d1,2)*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2)) - 
   4*c1*j2*powl(d1,2)*(8*a2*b1*c1 + 8*a1*b2*c1 + 
      8*a1*b1*c2 - 4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 
      2*d1*e2*f1 + 2*d1*e1*f2 - 4*b1*f1*f2 - 
      2*c2*powl(d1,2) - 2*a2*powl(e1,2) - 2*b2*powl(f1,2)) + 
   g1*g2*powl(e1,2)*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2)) - 
   4*a1*j2*powl(e1,2)*(8*a2*b1*c1 + 8*a1*b2*c1 + 
      8*a1*b1*c2 - 4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 
      2*d1*e2*f1 + 2*d1*e1*f2 - 4*b1*f1*f2 - 
      2*c2*powl(d1,2) - 2*a2*powl(e1,2) - 2*b2*powl(f1,2)) + 
   h1*h2*powl(f1,2)*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2)) - 
   4*b1*j2*powl(f1,2)*(8*a2*b1*c1 + 8*a1*b2*c1 + 
      8*a1*b1*c2 - 4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 
      2*d1*e2*f1 + 2*d1*e1*f2 - 4*b1*f1*f2 - 
      2*c2*powl(d1,2) - 2*a2*powl(e1,2) - 2*b2*powl(f1,2)) + 
   4*b2*c1*d1*g1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) - 
   4*b1*c1*d2*g1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) - 
   d1*e1*e2*g1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) - 
   2*b2*e1*f1*g1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) + 
   2*b1*e2*f1*g1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) - 
   8*a1*b2*c1*h1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) + 
   2*c1*d1*d2*h1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) + 
   2*a1*e1*e2*h1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) - 
   d2*e1*f1*h1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) - 
   d1*e2*f1*h1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) + 
   8*a1*b1*c1*h2*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) + 
   2*d1*e1*f1*h2*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) + 
   4*a1*b2*e1*i1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) - 
   d1*d2*e1*i1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) - 
   4*a1*b1*e2*i1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) - 
   2*b2*d1*f1*i1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) + 
   2*b1*d2*f1*i1*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) - 
   2*c1*h2*powl(d1,2)*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) + 
   e2*i1*powl(d1,2)*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) + 
   d2*g1*powl(e1,2)*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) - 
   2*a1*h2*powl(e1,2)*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) + 
   2*b2*h1*powl(f1,2)*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) - 
   2*b1*h2*powl(f1,2)*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2));
	rtc[2]=
		-2*c2*d1*e1*g1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*powl(d2,2)) + 
   2*c1*d1*e2*g1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*powl(d2,2)) + 
   4*b1*c2*f1*g1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*powl(d2,2)) - 
   e1*e2*f1*g1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*powl(d2,2)) - 
   4*b1*c1*f2*g1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*powl(d2,2)) + 
   4*a1*c2*e1*h1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*powl(d2,2)) - 
   4*a1*c1*e2*h1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*powl(d2,2)) - 
   2*c2*d1*f1*h1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*powl(d2,2)) + 
   2*c1*d1*f2*h1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*powl(d2,2)) - 
   e1*f1*f2*h1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*powl(d2,2)) - 
   8*a1*b1*c2*i1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*powl(d2,2)) + 
   2*a1*e1*e2*i1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*powl(d2,2)) - 
   d1*e2*f1*i1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*powl(d2,2)) - 
   d1*e1*f2*i1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*powl(d2,2)) + 
   2*b1*f1*f2*i1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*powl(d2,2)) + 
   8*a1*b1*c1*i2*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*powl(d2,2)) + 
   2*d1*e1*f1*i2*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*powl(d2,2)) + 
   2*c2*i1*powl(d1,2)*(-(d2*e2*g1) + 2*b2*f2*g1 - 
      d2*e1*g2 - d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 
      2*a2*e2*h1 - d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - 
      d2*f1*h2 - d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 
      4*a1*b2*i2 + 2*d1*d2*i2 + i1*powl(d2,2)) - 
   2*c1*i2*powl(d1,2)*(-(d2*e2*g1) + 2*b2*f2*g1 - 
      d2*e1*g2 - d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 
      2*a2*e2*h1 - d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - 
      d2*f1*h2 - d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 
      4*a1*b2*i2 + 2*d1*d2*i2 + i1*powl(d2,2)) + 
   f2*g1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*powl(d2,2))*powl(e1,2) - 
   2*a1*i2*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*powl(d2,2))*powl(e1,2) + 
   f2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*powl(d1,2))*(-4*b2*c1*g1 - 4*b1*c2*g1 + 
      2*e1*e2*g1 - 4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - 
      e2*f1*h1 - e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - 
      d2*e1*i1 - d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - 
      d1*e1*i2 + 2*b1*f1*i2 + g2*powl(e1,2)) - 
   8*a2*b1*c1*g1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2)) + 
   2*c1*d1*d2*g1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2)) - 
   d2*e1*f1*g1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2)) - 
   d1*e1*f2*g1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2)) + 
   2*b1*f1*f2*g1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2)) + 
   8*a1*b1*c1*g2*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2)) + 
   2*d1*e1*f1*g2*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2)) + 
   4*a2*c1*d1*h1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2)) - 
   4*a1*c1*d2*h1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2)) - 
   2*a2*e1*f1*h1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2)) + 
   2*a1*e1*f2*h1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2)) - 
   d1*f1*f2*h1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2)) - 
   2*a2*d1*e1*i1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2)) + 
   2*a1*d2*e1*i1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2)) + 
   4*a2*b1*f1*i1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2)) - 
   d1*d2*f1*i1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2)) - 
   4*a1*b1*f2*i1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2)) - 
   2*c1*g2*powl(d1,2)*(-4*b2*c2*g1 - 4*b2*c1*g2 - 
      4*b1*c2*g2 + 2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 
      2*c2*d1*h2 + 2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - 
      d2*e2*i1 + 2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 
      2*b2*f1*i2 + 2*b1*f2*i2 + g1*powl(e2,2)) + 
   f2*i1*powl(d1,2)*(-4*b2*c2*g1 - 4*b2*c1*g2 - 
      4*b1*c2*g2 + 2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 
      2*c2*d1*h2 + 2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - 
      d2*e2*i1 + 2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 
      2*b2*f1*i2 + 2*b1*f2*i2 + g1*powl(e2,2)) + 
   2*a2*g1*powl(e1,2)*(-4*b2*c2*g1 - 4*b2*c1*g2 - 
      4*b1*c2*g2 + 2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 
      2*c2*d1*h2 + 2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - 
      d2*e2*i1 + 2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 
      2*b2*f1*i2 + 2*b1*f2*i2 + g1*powl(e2,2)) - 
   2*a1*g2*powl(e1,2)*(-4*b2*c2*g1 - 4*b2*c1*g2 - 
      4*b1*c2*g2 + 2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 
      2*c2*d1*h2 + 2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - 
      d2*e2*i1 + 2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 
      2*b2*f1*i2 + 2*b1*f2*i2 + g1*powl(e2,2)) + 
   e2*h1*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*powl(d2,2))*powl(f1,2) - 
   2*b1*i2*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*powl(d2,2))*powl(f1,2) - 
   2*b1*g2*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2))*powl(f1,2) + 
   d2*h1*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2))*powl(f1,2) + 
   i2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*powl(d1,2))*(8*a2*b1*c1 + 8*a1*b2*c1 + 
      8*a1*b1*c2 - 4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 
      2*d1*e2*f1 + 2*d1*e1*f2 - 4*b1*f1*f2 - 
      2*c2*powl(d1,2) - 2*a2*powl(e1,2) - 2*b2*powl(f1,2)) + 
   g2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2))*
    (8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 4*c1*d1*d2 - 
      4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 2*d1*e1*f2 - 
      4*b1*f1*f2 - 2*c2*powl(d1,2) - 2*a2*powl(e1,2) - 
      2*b2*powl(f1,2)) + 
   e2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*powl(d1,2))*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) + 
   d2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2))*
    (2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - e1*f2*g1 + 
      2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 4*a1*c2*h1 + 
      2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 2*a1*e2*i1 - 
      d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - d1*f1*i2 + 
      h2*powl(f1,2)) + h2*
    (8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 4*c1*d1*d2 - 
      4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 2*d1*e1*f2 - 
      4*b1*f1*f2 - 2*c2*powl(d1,2) - 2*a2*powl(e1,2) - 
      2*b2*powl(f1,2))*(2*c2*d1*g1 + 2*c1*d2*g1 - 
      e2*f1*g1 - e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 
      4*a2*c1*h1 - 4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 
      2*a2*e1*i1 + 2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 
      2*a1*e1*i2 - d1*f1*i2 + h2*powl(f1,2)) - 
   4*b1*c1*g1*g2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2)) + 
   2*c1*d1*g2*h1*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2)) - 
   e1*f1*g2*h1*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2)) + 
   2*c1*d1*g1*h2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2)) - 
   e1*f1*g1*h2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2)) - 
   4*a1*c1*h1*h2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2)) - 
   d1*e1*g2*i1*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2)) + 
   2*b1*f1*g2*i1*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2)) + 
   2*a1*e1*h2*i1*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2)) - 
   d1*f1*h2*i1*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2)) - 
   d1*e1*g1*i2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2)) + 
   2*b1*f1*g1*i2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2)) + 
   2*a1*e1*h1*i2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2)) - 
   d1*f1*h1*i2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2)) - 
   4*a1*b1*i1*i2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2)) + 
   16*a1*b1*c1*j2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2)) + 
   4*d1*e1*f1*j2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2)) + 
   i1*i2*powl(d1,2)*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2)) - 
   4*c1*j2*powl(d1,2)*(8*a2*b2*c1 + 8*a2*b1*c2 + 
      8*a1*b2*c2 - 4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 
      2*d2*e1*f2 + 2*d1*e2*f2 - 4*b2*f1*f2 - 
      2*c1*powl(d2,2) - 2*a1*powl(e2,2) - 2*b1*powl(f2,2)) + 
   g1*g2*powl(e1,2)*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2)) - 
   4*a1*j2*powl(e1,2)*(8*a2*b2*c1 + 8*a2*b1*c2 + 
      8*a1*b2*c2 - 4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 
      2*d2*e1*f2 + 2*d1*e2*f2 - 4*b2*f1*f2 - 
      2*c1*powl(d2,2) - 2*a1*powl(e2,2) - 2*b1*powl(f2,2)) + 
   h1*h2*powl(f1,2)*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2)) - 
   4*b1*j2*powl(f1,2)*(8*a2*b2*c1 + 8*a2*b1*c2 + 
      8*a1*b2*c2 - 4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 
      2*d2*e1*f2 + 2*d1*e2*f2 - 4*b2*f1*f2 - 
      2*c1*powl(d2,2) - 2*a1*powl(e2,2) - 2*b1*powl(f2,2)) + 
   4*b2*c1*d1*g1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) - 
   4*b1*c1*d2*g1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) - 
   d1*e1*e2*g1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) - 
   2*b2*e1*f1*g1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) + 
   2*b1*e2*f1*g1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) - 
   8*a1*b2*c1*h1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) + 
   2*c1*d1*d2*h1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) + 
   2*a1*e1*e2*h1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) - 
   d2*e1*f1*h1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) - 
   d1*e2*f1*h1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) + 
   8*a1*b1*c1*h2*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) + 
   2*d1*e1*f1*h2*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) + 
   4*a1*b2*e1*i1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) - 
   d1*d2*e1*i1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) - 
   4*a1*b1*e2*i1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) - 
   2*b2*d1*f1*i1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) + 
   2*b1*d2*f1*i1*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) - 
   2*c1*h2*powl(d1,2)*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) + 
   e2*i1*powl(d1,2)*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) + 
   d2*g1*powl(e1,2)*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) - 
   2*a1*h2*powl(e1,2)*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) + 
   2*b2*h1*powl(f1,2)*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) - 
   2*b1*h2*powl(f1,2)*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) + 
   c2*powl(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
      2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
      2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
      d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
      4*a1*b1*i2 + i2*powl(d1,2),2) + 
   a2*powl(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2),2) + 
   j2*powl(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2),2) + 
   b2*powl(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - e1*f2*g1 + 
      2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 4*a1*c2*h1 + 
      2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 2*a1*e2*i1 - 
      d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - d1*f1*i2 + 
      h2*powl(f1,2),2);
	rtc[3]=
		2*c2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*powl(d1,2))*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - 
      d1*e2*g2 + 2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - 
      d2*f2*h1 + 2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - 
      d1*f2*h2 - 4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
      2*d1*d2*i2 + i1*powl(d2,2)) - 
   2*c2*d1*e1*g1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) + 
   2*c1*d1*e2*g1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) + 
   4*b1*c2*f1*g1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) - 
   e1*e2*f1*g1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) - 
   4*b1*c1*f2*g1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) + 
   4*a1*c2*e1*h1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) - 
   4*a1*c1*e2*h1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) - 
   2*c2*d1*f1*h1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) + 
   2*c1*d1*f2*h1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) - 
   e1*f1*f2*h1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) - 
   8*a1*b1*c2*i1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) + 
   2*a1*e1*e2*i1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) - 
   d1*e2*f1*i1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) - 
   d1*e1*f2*i1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) + 
   2*b1*f1*f2*i1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) + 
   8*a1*b1*c1*i2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) + 
   2*d1*e1*f1*i2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) + 
   2*c2*i1*powl(d1,2)*(-(d2*e2*g2) + 2*b2*f2*g2 + 
      2*a2*e2*h2 - d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) - 
   2*c1*i2*powl(d1,2)*(-(d2*e2*g2) + 2*b2*f2*g2 + 
      2*a2*e2*h2 - d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) + 
   f2*g1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2))*powl(e1,2) - 
   2*a1*i2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2))*powl(e1,2) + 
   f2*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*powl(d2,2))*(-4*b2*c1*g1 - 4*b1*c2*g1 + 
      2*e1*e2*g1 - 4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - 
      e2*f1*h1 - e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - 
      d2*e1*i1 - d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - 
      d1*e1*i2 + 2*b1*f1*i2 + g2*powl(e1,2)) + 
   f2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*powl(d1,2))*(-4*b2*c2*g1 - 4*b2*c1*g2 - 
      4*b1*c2*g2 + 2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 
      2*c2*d1*h2 + 2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - 
      d2*e2*i1 + 2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 
      2*b2*f1*i2 + 2*b1*f2*i2 + g1*powl(e2,2)) + 
   2*a2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2))*
    (-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 2*e1*e2*g2 + 
      2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 2*c1*d2*h2 - 
      e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 2*b2*f2*i1 - 
      d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 2*b1*f2*i2 + 
      g1*powl(e2,2)) - 8*a2*b1*c1*g1*
    (-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*powl(e2,2)) + 
   2*c1*d1*d2*g1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) - 
   d2*e1*f1*g1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) - 
   d1*e1*f2*g1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) + 
   2*b1*f1*f2*g1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) + 
   8*a1*b1*c1*g2*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) + 
   2*d1*e1*f1*g2*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) + 
   4*a2*c1*d1*h1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) - 
   4*a1*c1*d2*h1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) - 
   2*a2*e1*f1*h1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) + 
   2*a1*e1*f2*h1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) - 
   d1*f1*f2*h1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) - 
   2*a2*d1*e1*i1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) + 
   2*a1*d2*e1*i1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) + 
   4*a2*b1*f1*i1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) - 
   d1*d2*f1*i1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) - 
   4*a1*b1*f2*i1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) - 
   2*c1*g2*powl(d1,2)*(-4*b2*c2*g2 + 2*c2*d2*h2 - 
      e2*f2*h2 - d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) + 
   f2*i1*powl(d1,2)*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) + 
   2*a2*g1*powl(e1,2)*(-4*b2*c2*g2 + 2*c2*d2*h2 - 
      e2*f2*h2 - d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) - 
   2*a1*g2*powl(e1,2)*(-4*b2*c2*g2 + 2*c2*d2*h2 - 
      e2*f2*h2 - d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) + 
   e2*h1*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2))*powl(f1,2) - 
   2*b1*i2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2))*powl(f1,2) - 
   2*b1*g2*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2))*powl(f1,2) + 
   d2*h1*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*powl(e2,2))*powl(f1,2) + 
   i2*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*powl(d2,2))*(8*a2*b1*c1 + 8*a1*b2*c1 + 
      8*a1*b1*c2 - 4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 
      2*d1*e2*f1 + 2*d1*e1*f2 - 4*b1*f1*f2 - 
      2*c2*powl(d1,2) - 2*a2*powl(e1,2) - 2*b2*powl(f1,2)) + 
   g2*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2))*
    (8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 4*c1*d1*d2 - 
      4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 2*d1*e1*f2 - 
      4*b1*f1*f2 - 2*c2*powl(d1,2) - 2*a2*powl(e1,2) - 
      2*b2*powl(f1,2)) + 
   e2*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*powl(d2,2))*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
      e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
      4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 
      2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - 
      d1*f1*i2 + h2*powl(f1,2)) + 
   d2*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2))*
    (2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - e1*f2*g1 + 
      2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 4*a1*c2*h1 + 
      2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 2*a1*e2*i1 - 
      d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - d1*f1*i2 + 
      h2*powl(f1,2)) + i2*
    (-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*powl(d1,2))*(8*a2*b2*c1 + 8*a2*b1*c2 + 
      8*a1*b2*c2 - 4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 
      2*d2*e1*f2 + 2*d1*e2*f2 - 4*b2*f1*f2 - 
      2*c1*powl(d2,2) - 2*a1*powl(e2,2) - 2*b1*powl(f2,2)) + 
   g2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2))*
    (8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 4*c2*d1*d2 - 
      4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 2*d1*e2*f2 - 
      4*b2*f1*f2 - 2*c1*powl(d2,2) - 2*a1*powl(e2,2) - 
      2*b1*powl(f2,2)) + 
   2*j2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2))*
    (8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 4*c2*d1*d2 - 
      4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 2*d1*e2*f2 - 
      4*b2*f1*f2 - 2*c1*powl(d2,2) - 2*a1*powl(e2,2) - 
      2*b1*powl(f2,2)) + 
   h2*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - e1*f2*g1 + 
      2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 4*a1*c2*h1 + 
      2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 2*a1*e2*i1 - 
      d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - d1*f1*i2 + 
      h2*powl(f1,2))*(8*a2*b2*c1 + 8*a2*b1*c2 + 
      8*a1*b2*c2 - 4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 
      2*d2*e1*f2 + 2*d1*e2*f2 - 4*b2*f1*f2 - 
      2*c1*powl(d2,2) - 2*a1*powl(e2,2) - 2*b1*powl(f2,2)) - 
   4*b1*c1*g1*g2*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   2*c1*d1*g2*h1*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2)) - 
   e1*f1*g2*h1*(8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*powl(d2,2) - 
      2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   2*c1*d1*g1*h2*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2)) - 
   e1*f1*g1*h2*(8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*powl(d2,2) - 
      2*a2*powl(e2,2) - 2*b2*powl(f2,2)) - 
   4*a1*c1*h1*h2*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2)) - 
   d1*e1*g2*i1*(8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*powl(d2,2) - 
      2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   2*b1*f1*g2*i1*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   2*a1*e1*h2*i1*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2)) - 
   d1*f1*h2*i1*(8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*powl(d2,2) - 
      2*a2*powl(e2,2) - 2*b2*powl(f2,2)) - 
   d1*e1*g1*i2*(8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*powl(d2,2) - 
      2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   2*b1*f1*g1*i2*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   2*a1*e1*h1*i2*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2)) - 
   d1*f1*h1*i2*(8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*powl(d2,2) - 
      2*a2*powl(e2,2) - 2*b2*powl(f2,2)) - 
   4*a1*b1*i1*i2*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   16*a1*b1*c1*j2*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   4*d1*e1*f1*j2*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   i1*i2*powl(d1,2)*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2)) - 
   4*c1*j2*powl(d1,2)*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   g1*g2*powl(e1,2)*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2)) - 
   4*a1*j2*powl(e1,2)*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   h1*h2*powl(f1,2)*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2)) - 
   4*b1*j2*powl(f1,2)*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   e2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*powl(d1,2))*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) + 
   d2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2))*
    (2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 2*c1*d2*g2 - 
      e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 4*a2*c1*h2 - 
      4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - d2*f2*i1 + 
      2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - d1*f2*i2 + 
      h1*powl(f2,2)) + h2*
    (8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 4*c1*d1*d2 - 
      4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 2*d1*e1*f2 - 
      4*b1*f1*f2 - 2*c2*powl(d1,2) - 2*a2*powl(e1,2) - 
      2*b2*powl(f1,2))*(2*c2*d2*g1 - e2*f2*g1 + 
      2*c2*d1*g2 + 2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 
      4*a2*c2*h1 - 4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 
      2*a2*e2*i1 - d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - 
      d2*f1*i2 - d1*f2*i2 + h1*powl(f2,2)) + 
   2*b2*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - e1*f2*g1 + 
      2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 4*a1*c2*h1 + 
      2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 2*a1*e2*i1 - 
      d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - d1*f1*i2 + 
      h2*powl(f1,2))*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) + 
   4*b2*c1*d1*g1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) - 
   4*b1*c1*d2*g1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) - 
   d1*e1*e2*g1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) - 
   2*b2*e1*f1*g1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) + 
   2*b1*e2*f1*g1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) - 
   8*a1*b2*c1*h1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) + 
   2*c1*d1*d2*h1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) + 
   2*a1*e1*e2*h1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) - 
   d2*e1*f1*h1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) - 
   d1*e2*f1*h1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) + 
   8*a1*b1*c1*h2*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) + 
   2*d1*e1*f1*h2*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) + 
   4*a1*b2*e1*i1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) - 
   d1*d2*e1*i1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) - 
   4*a1*b1*e2*i1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) - 
   2*b2*d1*f1*i1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) + 
   2*b1*d2*f1*i1*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) - 
   2*c1*h2*powl(d1,2)*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) + 
   e2*i1*powl(d1,2)*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) + 
   d2*g1*powl(e1,2)*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) - 
   2*a1*h2*powl(e1,2)*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) + 
   2*b2*h1*powl(f1,2)*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) - 
   2*b1*h2*powl(f1,2)*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2));
	rtc[4]=
		2*c2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*powl(d1,2))*(-(d2*e2*g2) + 2*b2*f2*g2 + 
      2*a2*e2*h2 - d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) + 
   f2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - d2*f2*h2 - 
      4*a2*b2*i2 + i2*powl(d2,2))*
    (-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 4*b1*c1*g2 + 
      2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - e1*f2*h1 + 
      2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - d1*e2*i1 + 
      2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 2*b1*f1*i2 + 
      g2*powl(e1,2)) + f2*
    (-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*powl(d2,2))*(-4*b2*c2*g1 - 4*b2*c1*g2 - 
      4*b1*c2*g2 + 2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 
      2*c2*d1*h2 + 2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - 
      d2*e2*i1 + 2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 
      2*b2*f1*i2 + 2*b1*f2*i2 + g1*powl(e2,2)) + 
   f2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*powl(d1,2))*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) + 
   2*a2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2))*
    (-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*powl(e2,2)) + 
   i2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - d2*f2*h2 - 
      4*a2*b2*i2 + i2*powl(d2,2))*
    (8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 4*c1*d1*d2 - 
      4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 2*d1*e1*f2 - 
      4*b1*f1*f2 - 2*c2*powl(d1,2) - 2*a2*powl(e1,2) - 
      2*b2*powl(f1,2)) + 
   g2*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*powl(e2,2))*
    (8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 4*c1*d1*d2 - 
      4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 2*d1*e1*f2 - 
      4*b1*f1*f2 - 2*c2*powl(d1,2) - 2*a2*powl(e1,2) - 
      2*b2*powl(f1,2)) + 
   e2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - d2*f2*h2 - 
      4*a2*b2*i2 + i2*powl(d2,2))*
    (2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - e1*f2*g1 + 
      2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 4*a1*c2*h1 + 
      2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 2*a1*e2*i1 - 
      d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - d1*f1*i2 + 
      h2*powl(f1,2)) + d2*
    (-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*powl(e2,2))*
    (2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - e1*f2*g1 + 
      2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 4*a1*c2*h1 + 
      2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 2*a1*e2*i1 - 
      d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - d1*f1*i2 + 
      h2*powl(f1,2)) + i2*
    (-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*powl(d2,2))*(8*a2*b2*c1 + 8*a2*b1*c2 + 
      8*a1*b2*c2 - 4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 
      2*d2*e1*f2 + 2*d1*e2*f2 - 4*b2*f1*f2 - 
      2*c1*powl(d2,2) - 2*a1*powl(e2,2) - 2*b1*powl(f2,2)) + 
   g2*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2))*
    (8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 4*c2*d1*d2 - 
      4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 2*d1*e2*f2 - 
      4*b2*f1*f2 - 2*c1*powl(d2,2) - 2*a1*powl(e2,2) - 
      2*b1*powl(f2,2)) + 
   i2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*powl(d1,2))*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   g2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2))*
    (8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*powl(d2,2) - 
      2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   2*j2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
      4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 
      2*d1*e1*f2 - 4*b1*f1*f2 - 2*c2*powl(d1,2) - 
      2*a2*powl(e1,2) - 2*b2*powl(f1,2))*
    (8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*powl(d2,2) - 
      2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   h2*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - e1*f2*g1 + 
      2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 4*a1*c2*h1 + 
      2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 2*a1*e2*i1 - 
      d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - d1*f1*i2 + 
      h2*powl(f1,2))*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   e2*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*powl(d2,2))*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2)) + 
   d2*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2))*
    (2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 2*c1*d2*g2 - 
      e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 4*a2*c1*h2 - 
      4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - d2*f2*i1 + 
      2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - d1*f2*i2 + 
      h1*powl(f2,2)) + h2*
    (8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 4*c2*d1*d2 - 
      4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 2*d1*e2*f2 - 
      4*b2*f1*f2 - 2*c1*powl(d2,2) - 2*a1*powl(e2,2) - 
      2*b1*powl(f2,2))*(2*c2*d2*g1 - e2*f2*g1 + 
      2*c2*d1*g2 + 2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 
      4*a2*c2*h1 - 4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 
      2*a2*e2*i1 - d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - 
      d2*f1*i2 - d1*f2*i2 + h1*powl(f2,2)) + 
   e2*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 2*b1*f2*g1 - 
      d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 2*a1*e2*h1 - 
      d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - d1*f1*h2 - 
      4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 4*a1*b1*i2 + 
      i2*powl(d1,2))*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) + 
   d2*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
      4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
      e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
      d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
      2*b1*f1*i2 + g2*powl(e1,2))*
    (2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 2*a2*e2*i2 - 
      d2*f2*i2 + h2*powl(f2,2)) + 
   h2*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 4*c1*d1*d2 - 
      4*a1*e1*e2 + 2*d2*e1*f1 + 2*d1*e2*f1 + 2*d1*e1*f2 - 
      4*b1*f1*f2 - 2*c2*powl(d1,2) - 2*a2*powl(e1,2) - 
      2*b2*powl(f1,2))*(2*c2*d2*g2 - e2*f2*g2 - 
      4*a2*c2*h2 + 2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) + 
   2*b2*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - e1*f2*g1 + 
      2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 4*a1*c2*h1 + 
      2*f1*f2*h1 - 4*a1*c1*h2 + 2*a2*e1*i1 + 2*a1*e2*i1 - 
      d2*f1*i1 - d1*f2*i1 + 2*a1*e1*i2 - d1*f1*i2 + 
      h2*powl(f1,2))*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) + 
   c2*powl(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*powl(d2,2),2) + 
   a2*powl(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2),2) + 
   j2*powl(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2),2) + 
   b2*powl(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 
      2*c1*d2*g2 - e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 
      4*a2*c1*h2 - 4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - 
      d2*f2*i1 + 2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - 
      d1*f2*i2 + h1*powl(f2,2),2);
	rtc[5]=
		2*c2*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*powl(d2,2))*(-(d2*e2*g2) + 2*b2*f2*g2 + 
      2*a2*e2*h2 - d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2)) + 
   f2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - d2*f2*h2 - 
      4*a2*b2*i2 + i2*powl(d2,2))*
    (-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 2*e1*e2*g2 + 
      2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 2*c1*d2*h2 - 
      e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 2*b2*f2*i1 - 
      d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 2*b1*f2*i2 + 
      g1*powl(e2,2)) + f2*
    (-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*powl(d2,2))*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - 
      d2*e2*i2 + 2*b2*f2*i2 + g2*powl(e2,2)) + 
   2*a2*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2))*
    (-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*powl(e2,2)) + 
   i2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - d2*f2*h2 - 
      4*a2*b2*i2 + i2*powl(d2,2))*
    (8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 4*c2*d1*d2 - 
      4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 2*d1*e2*f2 - 
      4*b2*f1*f2 - 2*c1*powl(d2,2) - 2*a1*powl(e2,2) - 
      2*b1*powl(f2,2)) + 
   g2*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*powl(e2,2))*
    (8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 4*c2*d1*d2 - 
      4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 2*d1*e2*f2 - 
      4*b2*f1*f2 - 2*c1*powl(d2,2) - 2*a1*powl(e2,2) - 
      2*b1*powl(f2,2)) + 
   i2*(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*powl(d2,2))*(8*a2*b2*c2 + 2*d2*e2*f2 - 
      2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   g2*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2))*
    (8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*powl(d2,2) - 
      2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   2*j2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
      4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 
      2*d1*e2*f2 - 4*b2*f1*f2 - 2*c1*powl(d2,2) - 
      2*a1*powl(e2,2) - 2*b1*powl(f2,2))*
    (8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*powl(d2,2) - 
      2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   e2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - d2*f2*h2 - 
      4*a2*b2*i2 + i2*powl(d2,2))*
    (2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 2*c1*d2*g2 - 
      e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 4*a2*c1*h2 - 
      4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - d2*f2*i1 + 
      2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - d1*f2*i2 + 
      h1*powl(f2,2)) + d2*
    (-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*powl(e2,2))*
    (2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 2*c1*d2*g2 - 
      e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 4*a2*c1*h2 - 
      4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - d2*f2*i1 + 
      2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - d1*f2*i2 + 
      h1*powl(f2,2)) + h2*
    (8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*powl(d2,2) - 
      2*a2*powl(e2,2) - 2*b2*powl(f2,2))*
    (2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 2*c1*d2*g2 - 
      e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 4*a2*c1*h2 - 
      4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - d2*f2*i1 + 
      2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - d1*f2*i2 + 
      h1*powl(f2,2)) + e2*
    (-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
      2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
      2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
      4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 2*d1*d2*i2 + 
      i1*powl(d2,2))*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) + 
   d2*(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 
      2*e1*e2*g2 + 2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 
      2*c1*d2*h2 - e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 
      2*b2*f2*i1 - d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 
      2*b1*f2*i2 + g1*powl(e2,2))*
    (2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 2*a2*e2*i2 - 
      d2*f2*i2 + h2*powl(f2,2)) + 
   h2*(8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 4*c2*d1*d2 - 
      4*a2*e1*e2 + 2*d2*e2*f1 + 2*d2*e1*f2 + 2*d1*e2*f2 - 
      4*b2*f1*f2 - 2*c1*powl(d2,2) - 2*a1*powl(e2,2) - 
      2*b1*powl(f2,2))*(2*c2*d2*g2 - e2*f2*g2 - 
      4*a2*c2*h2 + 2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2)) + 
   2*b2*(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 2*c1*d2*g2 - 
      e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 4*a2*c1*h2 - 
      4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - d2*f2*i1 + 
      2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - d1*f2*i2 + 
      h1*powl(f2,2))*(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2));
	rtc[6]=
		f2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - d2*f2*h2 - 
      4*a2*b2*i2 + i2*powl(d2,2))*
    (-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*powl(e2,2)) + 
   i2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - d2*f2*h2 - 
      4*a2*b2*i2 + i2*powl(d2,2))*
    (8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*powl(d2,2) - 
      2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   g2*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*powl(e2,2))*
    (8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*powl(d2,2) - 
      2*a2*powl(e2,2) - 2*b2*powl(f2,2)) + 
   e2*(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - d2*f2*h2 - 
      4*a2*b2*i2 + i2*powl(d2,2))*
    (2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 2*a2*e2*i2 - 
      d2*f2*i2 + h2*powl(f2,2)) + 
   d2*(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*powl(e2,2))*
    (2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 2*a2*e2*i2 - 
      d2*f2*i2 + h2*powl(f2,2)) + 
   h2*(8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*powl(d2,2) - 
      2*a2*powl(e2,2) - 2*b2*powl(f2,2))*
    (2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 2*a2*e2*i2 - 
      d2*f2*i2 + h2*powl(f2,2)) + 
   c2*powl(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - 
      d2*f2*h2 - 4*a2*b2*i2 + i2*powl(d2,2),2) + 
   a2*powl(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
      2*b2*f2*i2 + g2*powl(e2,2),2) + 
   j2*powl(8*a2*b2*c2 + 2*d2*e2*f2 - 2*c2*powl(d2,2) - 
      2*a2*powl(e2,2) - 2*b2*powl(f2,2),2) + 
   b2*powl(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 
      2*a2*e2*i2 - d2*f2*i2 + h2*powl(f2,2),2);

	int order=0;
	for (int k=0;k<=6;k++){
	    if (rtc[k]!=0){
		order=k;
	    }
	}
	for (int k=0;k<=order;k++)
	    rtc[k]/=rtc[order];
	if (order==0)
	    return false; // usually order==6

	#ifdef DEBUG
	g_exceptioninf<<endl<<"g_iteration= "<<setw(16)<<g_iteration<<endl
		    <<"order="<<setw(10)<<order<<endl;
	#endif

	if (!zrhqr(rtc, order, rtr, rti)) // find roots for a polynomial using eigenvalues method
	    return false;

	#ifdef DEBUG
	for (int k=1;k<=order;k++){
	    g_exceptioninf<<setw(16)<<rtr[k]<<setw(16)<<rti[k]<<endl;
	}
	#endif

	long double lamda[6];
	int jj=0;
	for (int k=1;k<=order;k++){
	    if(fabsl(rti[k])<PREC)
		lamda[jj++]=rtr[k];
	}

	#ifdef DEBUG
	g_exceptioninf<<"jj="<<setw(10)<<jj<<endl;
	#endif

	long double x, y, z, det, within;
	bool found=false;
	for (int k=0;k<jj;k++){
	    x=0;
	    y=0;
	    z=0;
	    within=0;
            det=
	                8*a1*b1*c1 + 2*d1*e1*f1 - 2*c1*powl(d1,2) - 
			2*a1*powl(e1,2) - 2*b1*powl(f1,2) + 
			lamda[k]*(8*a2*b1*c1 + 8*a1*b2*c1 + 8*a1*b1*c2 - 
			4*c1*d1*d2 - 4*a1*e1*e2 + 2*d2*e1*f1 + 
			2*d1*e2*f1 + 2*d1*e1*f2 - 4*b1*f1*f2 - 
			2*c2*powl(d1,2) - 2*a2*powl(e1,2) - 2*b2*powl(f1,2))\
			+ (8*a2*b2*c1 + 8*a2*b1*c2 + 8*a1*b2*c2 - 
			4*c2*d1*d2 - 4*a2*e1*e2 + 2*d2*e2*f1 + 
			2*d2*e1*f2 + 2*d1*e2*f2 - 4*b2*f1*f2 - 
			2*c1*powl(d2,2) - 2*a1*powl(e2,2) - 2*b1*powl(f2,2))*
			powl(lamda[k],2) + (8*a2*b2*c2 + 2*d2*e2*f2 - 
			2*c2*powl(d2,2) - 2*a2*powl(e2,2) - 2*b2*powl(f2,2))*
	                powl(lamda[k],3);
	    if(det!=0){ // if determinant is zero, there are infinite solutions
		        // it is essential to skip this case, otherwise computatation will fail.
		x=
			((-4*b1*c1*g1 + 2*c1*d1*h1 - e1*f1*h1 - d1*e1*i1 + 
			2*b1*f1*i1 + g1*powl(e1,2) + 
			lamda[k]*(-4*b2*c1*g1 - 4*b1*c2*g1 + 2*e1*e2*g1 - 
			4*b1*c1*g2 + 2*c2*d1*h1 + 2*c1*d2*h1 - e2*f1*h1 - 
			e1*f2*h1 + 2*c1*d1*h2 - e1*f1*h2 - d2*e1*i1 - 
			d1*e2*i1 + 2*b2*f1*i1 + 2*b1*f2*i1 - d1*e1*i2 + 
			2*b1*f1*i2 + g2*powl(e1,2)) + 
			(-4*b2*c2*g1 - 4*b2*c1*g2 - 4*b1*c2*g2 + 2*e1*e2*g2 + 
			2*c2*d2*h1 - e2*f2*h1 + 2*c2*d1*h2 + 2*c1*d2*h2 - 
			e2*f1*h2 - e1*f2*h2 - d2*e2*i1 + 2*b2*f2*i1 - 
			d2*e1*i2 - d1*e2*i2 + 2*b2*f1*i2 + 2*b1*f2*i2 + 
			g1*powl(e2,2))*powl(lamda[k],2) + 
			(-4*b2*c2*g2 + 2*c2*d2*h2 - e2*f2*h2 - d2*e2*i2 + 
			 2*b2*f2*i2 + g2*powl(e2,2))*powl(lamda[k],3)) )/det;
		y=
			((2*c1*d1*g1 - e1*f1*g1 - 4*a1*c1*h1 + 2*a1*e1*i1 - 
			d1*f1*i1 + h1*powl(f1,2) + 
			lamda[k]*(2*c2*d1*g1 + 2*c1*d2*g1 - e2*f1*g1 - 
			e1*f2*g1 + 2*c1*d1*g2 - e1*f1*g2 - 4*a2*c1*h1 - 
			4*a1*c2*h1 + 2*f1*f2*h1 - 4*a1*c1*h2 + 
			2*a2*e1*i1 + 2*a1*e2*i1 - d2*f1*i1 - d1*f2*i1 + 
			2*a1*e1*i2 - d1*f1*i2 + h2*powl(f1,2)) + 
			(2*c2*d2*g1 - e2*f2*g1 + 2*c2*d1*g2 + 2*c1*d2*g2 - 
			e2*f1*g2 - e1*f2*g2 - 4*a2*c2*h1 - 4*a2*c1*h2 - 
			4*a1*c2*h2 + 2*f1*f2*h2 + 2*a2*e2*i1 - d2*f2*i1 + 
			2*a2*e1*i2 + 2*a1*e2*i2 - d2*f1*i2 - d1*f2*i2 + 
			h1*powl(f2,2))*powl(lamda[k],2) + 
			(2*c2*d2*g2 - e2*f2*g2 - 4*a2*c2*h2 + 2*a2*e2*i2 - 
			 d2*f2*i2 + h2*powl(f2,2))*powl(lamda[k],3)) )/det;
		z=
			((-(d1*e1*g1) + 2*b1*f1*g1 + 2*a1*e1*h1 - d1*f1*h1 - 
			4*a1*b1*i1 + i1*powl(d1,2) + 
			lamda[k]*(-(d2*e1*g1) - d1*e2*g1 + 2*b2*f1*g1 + 
			2*b1*f2*g1 - d1*e1*g2 + 2*b1*f1*g2 + 2*a2*e1*h1 + 
			2*a1*e2*h1 - d2*f1*h1 - d1*f2*h1 + 2*a1*e1*h2 - 
			d1*f1*h2 - 4*a2*b1*i1 - 4*a1*b2*i1 + 2*d1*d2*i1 - 
			4*a1*b1*i2 + i2*powl(d1,2)) + 
			(-(d2*e2*g1) + 2*b2*f2*g1 - d2*e1*g2 - d1*e2*g2 + 
			2*b2*f1*g2 + 2*b1*f2*g2 + 2*a2*e2*h1 - d2*f2*h1 + 
			2*a2*e1*h2 + 2*a1*e2*h2 - d2*f1*h2 - d1*f2*h2 - 
			4*a2*b2*i1 - 4*a2*b1*i2 - 4*a1*b2*i2 + 
			2*d1*d2*i2 + i1*powl(d2,2))*powl(lamda[k],2) + 
			(-(d2*e2*g2) + 2*b2*f2*g2 + 2*a2*e2*h2 - d2*f2*h2 - 
			 4*a2*b2*i2 + i2*powl(d2,2))*powl(lamda[k],3)) )/det;
	    
		within = a1*x*x+b1*y*y+c1*z*z+d1*x*y+e1*y*z+f1*z*x+g1*x+h1*y+i1*z+j1;
		
                #ifdef DEBUG
		g_exceptioninf<<setw(10)<<k<<setw(10)<<" det= "<<setw(16)<<det
			    <<" x y z= "<<setw(16)<<x<<setw(16)<<y<<setw(16)<<z
			    <<" within= "<<setw(16)<<within<<endl;
                #endif
		
		if(within<0){
		    point=vec(x,y,z);
		    found=true;
		    break;
		}
	    }
	    else {
		#ifdef DEBUG
		g_exceptioninf<<setw(10)<<k<<setw(10)<<" det= "<<setw(16)<<det
			    <<" ????????????????????????????????????????????????"
			    <<"????????????????????????????????"<<endl;
                #endif
	    }
	}
	return found;
}

} // namespace dem ends
