/* $Id: DetCheckT.cpp,v 1.2 2001-07-03 01:35:41 paklein Exp $ */
/* created: paklein (09/11/1997)                                          */

#include "DetCheckT.h"
#include <math.h>
#include "ExceptionCodes.h"
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "dArrayT.h"

/* constants */
const double Pi = acos(-1.0);

/* constructor */
DetCheckT::DetCheckT(const dSymMatrixT& s_jl, const dMatrixT& c_ijkl):
	fs_jl(s_jl),
	fc_ijkl(c_ijkl)
{

}

/* inline functions */

/* determinant function and derivatives */
inline double DetCheckT::det(double t) const
{
	return A0 + A2*sin(2.0*(phi2+t)) + A4*sin(4.0*(phi4+t));
}

inline double DetCheckT::ddet(double t) const
{
	return 2.0*A2*cos(2.0*(phi2+t)) + 4.0*A4*cos(4.0*(phi4+t));
}

inline double DetCheckT::dddet(double t) const
{
	return -4.0*A2*sin(2.0*(phi2+t)) - 16.0*A4*sin(4.0*(phi4+t));
}

/*
* Returns 1 if acoustic tensor isn't positive definite,
* and returns the normal to the surface of localization.
* Returns 0, otherwise.
*/
int DetCheckT::IsLocalized(dArrayT& normal)
{
	if (fs_jl.Rows() == 2)
		return DetCheck2D(normal);
	else
	//TEMP - not implemented
		return 0;
}

/**********************************************************************
* Private
**********************************************************************/

/* 2D determinant check function */
int DetCheckT::DetCheck2D(dArrayT& normal)
{
	/* set det function */
	ComputeCoefficients();

	/* quick exit */
	if ( A0 > 0 && A0 > (fabs(A2) + fabs(A4)) ) return 0;

	/* local minima */
	double min2  = ((A2 > 0.0) ? 3.0*Pi/4.0 : Pi/4.0) - phi2;

	double dat[4] = {-Pi/8.0, 3.0*Pi/8.0, 7.0*Pi/8.0, 11.0*Pi/8.0};

	dArrayT mins4(4,dat);
	if (A4 > 0) mins4 -= phi4;
	else mins4 -= (phi4 + Pi/4.0);

	/* closest local minima in sin2t and sin4t*/
	int    dex;
	double min = Pi;
	for (int i = 0; i < 4; i++)
	{	
		double diff = fabs(mins4[i] - min2);

		if (diff < min)
		{	
			dex = i;
			min = diff;
		}
	}
	
	/* search starting value */
	double angle = 0.5*(min2 + mins4[dex]);
	double maxstep = Pi/18.0; //10 degrees max per step
	
	/* Newton search */
	double res  = ddet(angle);
	double res0 = res;
	int   count = 0;
	while ( fabs(res/res0) > 1.0e-10 && ++count < 15)
	{
		double dangle = res/dddet(angle);
		dangle = (fabs(dangle) > maxstep) ? maxstep*((dangle < 0) ? -1:1) : dangle;
	
		angle -= dangle;
		res  = ddet(angle);
	}

	/* check minima */
	if (dddet(angle) < 0.0)
	{
		cout << "\n DetCheckT::IsLocalized: ERROR:\n";
		cout << "f =  " << A0 << " + " << A2 << "*sin(2.0*(";
		cout << phi2 << "+t)) + " << A4 << "*sin(4.0*(";
		cout << phi4 << "+t))" << '\n';
		cout << "root found at " << angle << " rad" << endl;
		cout << "starting at   " << 0.5*(min2 + mins4[dex]) << " rad" << endl;


		/* search starting value */
		double angle = 0.5*(min2 + mins4[dex]);
		
		/* Newton search */
		double res  = ddet(angle);
		double res0 = res;
		int   count = 0;
		while( fabs(res/res0) > 1.0e-10 && ++count < 10)
		{
			double dangle = res/dddet(angle);

			cout << count << '\t' << res << '\t' << dangle;
		
			angle -= dangle;
			res  = ddet(angle);
		}
			throw eGeneralFail;
	}
	
	if (det(angle) > 0.0)
		return 0;
	else
	{
		normal.Allocate(2);

		/* angle in [0, Pi] */
		if (angle > Pi) angle -= Pi;
		
		/* compute normal */
		normal[0] = cos(angle);
		normal[1] = sin(angle);
		
		return 1;
	}
}

/* compute coefficients of det(theta) function */
void DetCheckT::ComputeCoefficients(void)
{
	/* moduli components */
	double c11 = fc_ijkl(0,0);
	double c22 = fc_ijkl(1,1);
	double c33 = fc_ijkl(2,2);
	double c23 = fc_ijkl(1,2);
	double c13 = fc_ijkl(0,2);
	double c12 = fc_ijkl(0,1);
	
	/* stress components */
	double s11 = fs_jl[0];
	double s22 = fs_jl[1];
	double s12 = fs_jl[2];

	/* intermediate values */
	double s2t = (-(c12*c13) + c13*c22 + c11*c23 - c12*c23 + c13*s11 +
			c23*s11 + c11*s12 + c22*s12 + 2*c33*s12 + 2*s11*s12 + c13*s22 +
			c23*s22 + 2*s12*s22)/2;
	double s4t = (-(c12*c13) - c13*c22 + c11*c23 + c12*c23 + c13*s11 +
			c23*s11 + c11*s12 - c22*s12 + 2*s11*s12 - c13*s22 - c23*s22 -
			2*s12*s22)/4;
			
	double c2t = (-c13*c13 + c23*c23 + c11*c33 - c22*c33 + c11*s11 + c33*s11
			+ s11*s11 - c22*s22 - c33*s22 - s22*s22)/2;
	double c4t = (c12*c12 - c13*c13 - c11*c22 - 2*c13*c23 - c23*c23 +
			c11*c33 + 2*c12*c33 + c22*c33 + c11*s11 - c22*s11 + s11*s11 -
			4*c13*s12 - 4*c23*s12 - 4*s12*s12 - c11*s22 + c22*s22 - 2*s11*s22
			+ s22*s22)/8;

	/* phase shifts */
	phi2 = atan2(c2t,s2t)/2.0;
	phi4 = atan2(c4t,s4t)/4.0;
	
	/* amplitudes */
	A0 = (-c12*c12 - 3*c13*c13 + c11*c22 + 2*c13*c23 - 3*c23*c23 +
			3*c11*c33 - 2*c12*c33 + 3*c22*c33 + 3*c11*s11 + c22*s11 +
			4*c33*s11 + 3*s11*s11 + 4*c13*s12 + 4*c23*s12 + 4*s12*s12 +
			c11*s22 + 3*c22*s22 + 4*c33*s22 + 2*s11*s22 + 3*s22*s22)/8;

	A2 = c2t/sin(2.0*phi2);
	A4 = c4t/sin(4.0*phi4);
}
