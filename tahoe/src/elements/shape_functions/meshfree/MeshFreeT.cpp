/* $Id: MeshFreeT.cpp,v 1.1.1.1 2001-01-29 08:20:31 paklein Exp $ */
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
