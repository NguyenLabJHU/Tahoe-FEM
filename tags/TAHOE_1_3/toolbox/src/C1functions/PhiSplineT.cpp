/* $Id: PhiSplineT.cpp,v 1.4 2003-11-21 22:41:27 paklein Exp $ */
/* created: paklein (01/30/2000)                                          */
/* PhiSplineT.h                                                           */

#include "PhiSplineT.h"

/* constructor */

using namespace Tahoe;

PhiSplineT::PhiSplineT(const dArray2DT& points, FixityT fixity, double r_cut):
	CubicSplineT(points, fixity),
	fr_cut(r_cut)
{

}

/* returning values */
double PhiSplineT::Function(double r) const
{
	if (r > fr_cut)
		return 0.0;
	else
	{
		double z = function(r);
		return z*z/r;
	}
}

double PhiSplineT::DFunction(double r) const
{
	if (r > fr_cut)
		return 0.0;
	else
	{
		double z, Dz, DDz;
		all_functions(r, z, Dz, DDz);
		return z*(2.0*Dz - z/r)/r;	
	}
}

double PhiSplineT::DDFunction(double r) const
{
	if (r > fr_cut)
		return 0.0;
	else
	{
		double z, Dz, DDz;
		all_functions(r, z, Dz, DDz);

		double dphi = z*(2.0*Dz - z/r)/r;
		return 2.0*(Dz*Dz + z*DDz - dphi)/r;
	}
}

/* returning values in groups - returns refence to out to allow:
*
*	dArrayT& goodname = pfunc->MapFunction(in, tempspace);
*/
dArrayT& PhiSplineT::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension check */
	if ( in.Length() != out.Length() ) throw ExceptionT::kGeneralFail;
	
	const double *pin   = in.Pointer();
	double *pout  = out.Pointer();
	int    length = in.Length();
	
	/* fast mapping */
	for (int i = 0; i < length; i++)
		*pout++ = PhiSplineT::Function(*pin++);

	return out;
}

dArrayT& PhiSplineT::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension check */
	if ( in.Length() != out.Length() ) throw ExceptionT::kGeneralFail;
	
	const double *pin   =  in.Pointer();
	double *pout  = out.Pointer();
	int    length = in.Length();
	
	/* fast mapping */
	for (int i = 0; i < length; i++)
		*pout++ = PhiSplineT::DFunction(*pin++);

	return out;
}

dArrayT& PhiSplineT::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension check */
	if ( in.Length() != out.Length() ) throw ExceptionT::kGeneralFail;
	
	const double *pin   =  in.Pointer();
	double *pout  = out.Pointer();
	int    length = in.Length();
	
	/* fast mapping */
	for (int i = 0; i < length; i++)
		*pout++ = PhiSplineT::DDFunction(*pin++);

	return out;
}

/*
* Return 0th, 1st, and 2nd derivative in the respective
* fields of the dArrayT.
*/  	
void PhiSplineT::SetAll(double r, dArrayT& data) const
{
	if (r > fr_cut)
	{
		data[0] = 0.0;
		data[1] = 0.0;
		data[2] = 0.0;
	}
	else
	{
		double z, Dz, DDz;
		all_functions(r, z, Dz, DDz);

		double   phi = z*z/r;
		double  Dphi = (2.0*z*Dz - phi)/r;
		double DDphi = 2.0*(Dz*Dz + z*DDz - Dphi)/r;
	
		data[0] = phi;
		data[1] = Dphi;
		data[2] = DDphi;
	}
}
