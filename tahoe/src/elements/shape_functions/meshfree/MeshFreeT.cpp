/* $Id: MeshFreeT.cpp,v 1.1.1.1.4.1 2001-06-18 20:13:42 paklein Exp $ */
/* created: paklein (12/08/1999)                                          */

#include "MeshFreeT.h"

#include <iostream.h>
#include "ExceptionCodes.h"

istream& operator>>(istream& in, MeshFreeT::FormulationT& code)
{
	int i_code;
	in >> i_code;
	switch (i_code)
	{
		case MeshFreeT::kEFG:
			code = MeshFreeT::kEFG;
			break;
		case MeshFreeT::kRKPM:
			code = MeshFreeT::kRKPM;
			break;
		default:
			cout << "\n operator>>MeshFreeT::FormulationT: unknown code: "
			<< i_code<< endl;
			throw eBadInputValue;	
	}
	return in;
}

istream& operator>>(istream& in, MeshFreeT::WindowTypeT& code)
{
	int i_code;
	in >> i_code;
	switch (i_code)
	{
		case MeshFreeT::kGaussian:
			code = MeshFreeT::kGaussian;
			break;
		case MeshFreeT::kCubicSpline:
			code = MeshFreeT::kCubicSpline;
			break;
		case MeshFreeT::kBrick:
			code = MeshFreeT::kBrick;
			break;
		default:
			cout << "\n operator>>MeshFreeT::WindowTypeT: unknown code: "
			<< i_code<< endl;
			throw eBadInputValue;	
	}
	return in;
}
