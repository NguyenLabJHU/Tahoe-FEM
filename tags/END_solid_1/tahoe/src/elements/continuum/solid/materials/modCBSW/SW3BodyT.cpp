/* $Id: SW3BodyT.cpp,v 1.1.1.1 2001-01-29 08:20:26 paklein Exp $ */
/* created: paklein (05/22/1997)                                          */

#include "SW3BodyT.h"
#include <math.h>
#include "SWDataT.h"
#include "dMatrixT.h"
#include "iArray2DT.h"
#include "ThermalDilatationT.h"

/* parameters */
const int kNumVars = 3; //number of arguments in Phi

/* constructor */
SW3BodyT::SW3BodyT(const dArrayT& lengths,
	const dArrayT& angles, const iArray2DT& bondpairs,
	const ThermalDilatationT* thermal, const SWDataT& SW):
	ThreeBodyT(lengths, angles, bondpairs, thermal),
	fSW(SW)
{

}

/* triggers recomputation */
void SW3BodyT::Set(void)
{
	/* shallow workspace */
	dMatrixT sh_mat;
	dArrayT  sh_vec;

	/* expansion factor */
	double a = fSW.fa;
	if (fThermal) a *= (1.0 + fThermal->PercentElongation());
		
	for (int i = 0; i < fPairs.MajorDim(); i++)
	{
		double r1  = fLengths[fPairs(i,0)];
		double r2  = fLengths[fPairs(i,1)];
		double c12 = fAngles[i];
	
		/* phi */
		fPhi[i] = ComputePhi(r1, r2, c12, a);
		
		/* gradient */		
		sh_vec.Set(kNumVars, fdPhi(i));
		ComputeGradient(sh_vec, r1, r2, c12, a);

		/* hessian */		
		sh_mat.Set(kNumVars, kNumVars, fddPhi(i));
		ComputeHessian(sh_mat, r1, r2, c12, a);
	}
}

/**********************************************************************
* Private
**********************************************************************/

/* return the 3-body potential */	
double SW3BodyT::ComputePhi(double r1, double r2, double c12, double a) const
{
	return( pow(1./3 + c12,2)*fSW.feps*fSW.flambda*
	         ( (r1 < fSW.frcut*a) ?
	         	exp(fSW.fgamma/(-fSW.frcut + r1/a)) : 0.0 )*
( (r2 < fSW.frcut*a) ?
	exp(fSW.fgamma/(-fSW.frcut + r2/a)) : 0.0 ) );
}

/* compute Gradient of the 3-body potential */	
void SW3BodyT::ComputeGradient(dArrayT& grad, double r1, double r2,
	double c12, double a)
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9;

	z1 = 0.3333333333333333 + c12;
	z2 = pow(a,-1.);
	z3 = -fSW.frcut;
	z4 = pow(z1,2.);
	z5 = r1*z2;
	z6 = r2*z2;
	z5 = z3 + z5;
	z3 = z3 + z6;
	z6 = pow(z5,-2.);
	z5 = pow(z5,-1.);
	z7 = pow(z3,-2.);
	z3 = pow(z3,-1.);
	z5 = fSW.fgamma*z5;
	z3 = fSW.fgamma*z3;
	z5 = exp(z5);
	z3 = exp(z3);
	z8 = (r1 < fSW.frcut) ? z5 : 0;
	z9 = (r2 < fSW.frcut) ? z3 : 0;
	z1 = 2.*fSW.feps*fSW.flambda*z1*z8*z9;
	z5 = (r1 < fSW.frcut) ? -fSW.fgamma*z2*z5*z6 : 0;
	z2 = (r2 < fSW.frcut) ? -fSW.fgamma*z2*z3*z7 : 0;
	z3 = fSW.feps*fSW.flambda*z4;
	z4 = z3*z5*z9;
	z2 = z2*z3*z8;

	//returned	
	//List(z4,z2,z1);
	
	grad[0] = z4;
	grad[1] = z2;
	grad[2] = z1;
}

/* compute Hessian of the 3-body potential */
void SW3BodyT::ComputeHessian(dMatrixT& hessian, double r1, double r2,
	double c12, double a)
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18;

	z1 = 0.3333333333333333 + c12;
	z2 = pow(a,-2.);
	z3 = pow(a,-1.);
	z4 = pow(fSW.fgamma,2.);
	z5 = -fSW.frcut;
	z6 = pow(z1,2.);
	z7 = r1*z3;
	z8 = r2*z3;
	z7 = z5 + z7;
	z5 = z5 + z8;
	z8 = pow(z5,-4.);
	z9 = pow(z5,-3.);
	z10 = pow(z5,-2.);
	z5 = pow(z5,-1.);
	z11 = pow(z7,-4.);
	z12 = pow(z7,-3.);
	z13 = pow(z7,-2.);
	z7 = pow(z7,-1.);
	z5 = fSW.fgamma*z5;
	z7 = fSW.fgamma*z7;
	z5 = exp(z5);
	z7 = exp(z7);
	z14 = (r1 < fSW.frcut) ? z7 : 0;
	z15 = (r2 < fSW.frcut) ? z5 : 0;
	z16 = 2.*fSW.fgamma*z2;
	z17 = z2*z5;
	z18 = 2.*fSW.feps*fSW.flambda*z14*z15;
	z9 = z16*z5*z9;
	z12 = z12*z16*z7;
	z2 = z11*z2*z4*z7;
	z7 = (r1 < fSW.frcut) ? -fSW.fgamma*z13*z3*z7 : 0;
	z3 = (r2 < fSW.frcut) ? -fSW.fgamma*z10*z3*z5 : 0;
	z1 = 2.*fSW.feps*fSW.flambda*z1;
	z5 = z1*z15*z7;
	z1 = z1*z14*z3;
	z4 = z17*z4*z8;
	z3 = fSW.feps*fSW.flambda*z3*z6*z7;
	z2 = (r1 < fSW.frcut) ? z12 + z2 : 0;
	z2 = fSW.feps*fSW.flambda*z15*z2*z6;
	z4 = (r2 < fSW.frcut) ? z4 + z9 : 0;

	//returned
	//{{z2, z3, z5},
	// {z3, fSW.feps fSW.flambda z14 z4 z6, z1},
	// {z5, z1, z18}}

	hessian(0,0) = z2;
	hessian(1,1) = fSW.feps*fSW.flambda*z14*z4*z6;
	hessian(2,2) = z18;
	hessian(2,1) = hessian(1,2) = z1;
	hessian(2,0) = hessian(0,2) = z5;
	hessian(1,0) = hessian(0,1) = z3;
}
