/* $Id: MeshFreeT.cpp,v 1.6 2004-03-04 08:54:29 paklein Exp $ */
/* created: paklein (12/08/1999)                                          */

#include "MeshFreeT.h"

#include <iostream.h>
#include "ExceptionT.h"

namespace Tahoe {

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
			throw ExceptionT::kBadInputValue;	
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
		case MeshFreeT::kRectCubicSpline:
			code = MeshFreeT::kRectCubicSpline;
			break;
		case MeshFreeT::kBrick:
			code = MeshFreeT::kBrick;
			break;
		case MeshFreeT::kCubicSpline:
			code = MeshFreeT::kCubicSpline;
			break;
		default:
			cout << "\n operator>>MeshFreeT::WindowTypeT: unknown code: "
			<< i_code<< endl;
			throw ExceptionT::kBadInputValue;	
	}
	return in;
}

}
